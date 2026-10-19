# -*- coding: utf-8 -*-

"""Non-graphical part of the Thermal Conductivity step in a SEAMM flowchart.

The thermal conductivity is calculated from the fluctuations of the heat flux in
equilibrium MD, by the Green-Kubo integral of its autocorrelation function and by the
Helfand moments. The MD is run before this step -- in a Loop for independent runs --
or read from files; the step analyzes it.
"""

import importlib
import logging
from math import floor, log10
from pathlib import Path
import warnings

import numpy as np
from seamm_md_analysis import (
    Estimate,
    collective_helfand_moments,
    component_check,
    fit_green_kubo_integral,
    fit_helfand_moment,
    helfand_slope,
    isotropic_average,
    isotropy,
    effective_correlation_functions,
    flux_correlation_matrices,
    flux_helfand_matrices,
    mean_and_stderr,
    open_md,
    schur_complement,
    tensor_correlation_functions,
)
from seamm_md_analysis.requirements import Need, Requirements
from seamm_md_analysis.step_support import (
    check_md_flowchart,
    md_inputs,
    print_details,
    read_files,
    report_state,
    size_warnings,
    state_data,
)
from seamm_md_analysis.plots import (
    plot_correlation_functions,
    plot_GK_integrals,
    plot_helfand_moments,
    plot_helfand_slopes,
)
from tabulate import tabulate

from .thermal_conductivity_parameters import ThermalConductivityParameters
import thermal_conductivity_step
import molsystem
import seamm
import seamm_util
from seamm_util import Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)

#: The shortest run, or block of a run, that can be analyzed: the correlation
#: functions are kept to a tenth of it, and fitting them needs a little over 0.8 ps.
MIN_RUN_TIME = Q_(10.0, "ps")
#: The methods, in the order of the tables
METHODS = ("Green-Kubo", "Helfand moments", "Helfand derivative")
#: The names of the diagonal and off-diagonal components
DIAGONAL = ("κ_xx", "κ_yy", "κ_zz")
OFF_DIAGONAL = ("κ_xy", "κ_xz", "κ_yz")
#: The results stored for each method: the key of κ and its standard error
RESULT_KEYS = {
    "Green-Kubo": "kappa",
    "Helfand moments": "kappa (Helfand)",
    "Helfand derivative": "kappa (Helfand derivative)",
}

job = printing.getPrinter()
printer = printing.getPrinter("Thermal Conductivity")

# Add this module's properties to the standard properties
path = importlib.resources.files("thermal_conductivity_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


def _fit(function, *args, **kwargs):
    """The value of a fit, or NaN if it fails."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = function(*args, **kwargs)
    except Exception:
        return float("nan")
    value = result.value if hasattr(result, "value") else result.slope
    return float(value)


class ThermalConductivity(seamm.Node):
    """
    The non-graphical part of a Thermal Conductivity step in a flowchart.

    Attributes
    ----------
    parameters : ThermalConductivityParameters
        The control parameters for Thermal Conductivity.

    See Also
    --------
    TkThermalConductivity,
    ThermalConductivity, ThermalConductivityParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Thermal Conductivity",
        extension=None,
        logger=logger,
    ):
        """A step for Thermal Conductivity in a SEAMM flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.
        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes
        """
        logger.debug(f"Creating Thermal Conductivity {self}")
        super().__init__(
            flowchart=flowchart,
            title="Thermal Conductivity",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = thermal_conductivity_step.metadata
        self.parameters = ThermalConductivityParameters()
        self.tensor_labels = [
            ("xx", "red", "rgba(255,0,0,0.1)"),
            ("yy", "green", "rgba(0,255,0,0.1)"),
            ("zz", "blue", "rgba(0,0,255,0.1)"),
            ("xy", "rgb(127,127,0)", "rgba(127,127,0,0.1)"),
            ("xz", "rgb(255,0,255)", "rgba(255,0,255,0.1)"),
            ("yz", "rgb(0,255,255)", "rgba(0,255,255,0.1)"),
        ]
        self._reset()

        # Set our citation level to 1!
        self.citation_level = 1

    def _reset(self):
        """Clear the data of earlier runs, e.g. from an earlier loop iteration."""
        self.paths = []  # The md.h5md file of each run
        self.V = []  # The volume of each run or block
        self.T = []  # The temperature of each run or block
        self.timestep = []  # The interval between samples of the heat flux, fs
        self.Ms = []  # The Helfand moments of each run or block, 6 x m
        self.M = None  # The mean Helfand moments, 6 x m
        self.M_err = None  # Their spread over the runs or blocks, 6 x m
        self.Jcfs = []  # The correlation functions of J, 6 x m, of each
        self.Jcf = None  # Their mean, 6 x m
        self.Jcf_err = None  # Their spread, 6 x m
        self.GK_integrals = []  # The running integrals of the Jcf's, of each
        self.GK_integral = None  # Their mean
        self.GK_integral_err = None  # Their spread
        # For a mixture, the Green-Kubo integrals without the correction for the
        # interdiffusion of the species, of each run or block
        self.raw_GK_integrals = []
        self._species = None  # The species of a mixture whose fluxes were sampled
        # For a mixture, the matrices of the heat flux and species fluxes of each run
        # or block: the correlation functions, their integrals, and Helfand moments
        self._C_matrices = []
        self._L_matrices = []
        self._M_matrices = []
        self._total_time = Q_(0.0, "ps")  # The length of all the runs
        self._isotropy = None  # Whether the runs are isotropic, and why

    @property
    def version(self):
        """The semantic version of this module."""
        return thermal_conductivity_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return thermal_conductivity_step.__git_revision__

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "thermal-conductivity-step"
        parser = seamm_util.getParser()

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        super().create_parser(name=parser_name)

        if not parser_exists:
            # Any options for thermal conductivity itself
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )

        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        short : bool
            Unused: the description is always short.

        Returns
        -------
        str
            A description of the current step.
        """
        if P is None:
            P = self.parameters.values_to_dict()

        source = P["md source"]
        text = (
            "Calculate the thermal conductivity from the heat flux, by the Green-Kubo "
            "and Helfand methods, from "
        )
        if source == read_files:
            text += f"the MD in {P['md files']}."
        else:
            text += f"{source} earlier in the flowchart."
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Run a Thermal Conductivity step.

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)

        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Print what we are doing
        printer.important(__(self.description_text(P), indent=self.indent))

        _, configuration = self.get_system_configuration()

        self._reset()
        self.paths = md_inputs(self, P, name="heat_flux")
        if len(self.paths) == 0:
            raise RuntimeError("No MD run finished, so there is nothing to analyze.")
        runs = [self._read_md_file(path) for path in self.paths]
        if self._species is not None:
            self.references.cite(
                raw=self._bibliography["Bertossa2019"],
                alias="Bertossa2019",
                module="thermal_conductivity_step",
                level=1,
                note="The thermal conductivity of a mixture.",
            )

        text = f"Analyzing {len(runs)} MD run{'s' if len(runs) > 1 else ''}"
        if len(runs) == 1:
            source = "blocks"
            samples = self._split_into_blocks(runs[0], P["number of blocks"])
            text += f", divided into {len(samples)} blocks to estimate the errors"
        else:
            source = "runs"
            samples = runs
        text += ": " + ", ".join(str(p) for p in self.paths)
        printer.job(__(text, indent=self.indent + 4 * " "))
        printer.job("")

        for i, data in enumerate(samples, start=1):
            self.process_run(i, data, what="Block" if source == "blocks" else "Run")
        self._combine()

        self.analyze(
            P=P,
            configuration=configuration,
            source=source,
            uncorrected=self._uncorrected_mixture(configuration, runs),
        )

        return next_node

    def analyze(
        self,
        indent="",
        P=None,
        configuration=None,
        source="runs",
        uncorrected=None,
        **kwargs,
    ):
        """Fit, print, plot and store the thermal conductivity.

        The value of each method is the average over the runs (or blocks) of the
        average of the three diagonal components fitted to each run; its error is
        the standard error of these averages. The curves averaged over the runs are
        fitted too, for the graphs.

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        P : dict
            The values of the parameters.
        configuration : molsystem._Configuration
            The configuration for the results.
        source : str
            What the samples are: "runs" or "blocks".
        uncorrected : [str]
            The species of a mixture whose mass fluxes were not sampled, so the
            interdiffusion could not be removed; None for a pure substance.
        """
        n_samples = len(self.Ms)
        dt = Q_(self.timestep[0], "fs")
        ts = (np.arange(self.M.shape[1]) * dt.m_as("ps")).tolist()
        notes = []

        # The derivative of the averaged Helfand moments, fitted like the Green-Kubo
        # integral. Its time constants set where the Helfand moments are fitted.
        slopes = []
        xs = []
        errs = []
        fit0 = []
        for i in range(6):
            slope, x, err = helfand_slope(self.M[i], ts, sigma=self.M_err[i])
            slopes.append(slope)
            xs.append(x)
            errs.append(err)
            if i < 3:
                kappa, fit_err, a, a_err, tau, tau_err, tf, yf = (
                    fit_green_kubo_integral(slope, x, sigma=err)
                )
                fit0.append(
                    {
                        "kappa": kappa,
                        "stderr": fit_err,
                        "xs": tf,
                        "ys": yf,
                        "a": a,
                        "a_err": a_err,
                        "tau": tau,
                        "tau_err": tau_err,
                    }
                )
        self._plot_helfand_derivatives(slopes, xs, errs, fit0)

        starts = [max(1.0, 3 * max(f["tau"])) for f in fit0]
        starts.extend([max(starts)] * 3)

        # The slopes of the averaged Helfand moments
        fit = []
        for i in range(6):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    slope, err, fxs, fys = fit_helfand_moment(
                        self.M[i],
                        ts,
                        sigma=self.M_err[i],
                        start=starts[i],
                        logger=self.logger,
                    )
                fit.append({"kappa": slope, "stderr": err, "xs": fxs, "ys": fys})
            except Exception as e:
                logger.warning(f"The fit of the Helfand moments failed: {e}")
                if i == 0:
                    notes.append(
                        "The Helfand moments could not be fitted. The fit starts at "
                        f"{starts[i]:.2f} ps, but the moments extend only to "
                        f"{ts[-1]:.2f} ps, a tenth of the length of each run (or "
                        "block). Use longer runs, or fewer blocks."
                    )
                break
        self.plot_helfand_moment(self.M, ts, M_err=self.M_err, fit=fit)
        helfand_fitted = len(fit) > 0

        # The Green-Kubo integrals averaged over the runs
        cf_ts = (np.arange(self.Jcf.shape[1]) * dt.m_as("ps")).tolist()
        end = len(cf_ts) // 4
        self.plot_JJ_correlation(
            self.Jcf[:, :end], cf_ts[:end], err=self.Jcf_err[:, :end]
        )
        gk_ts = (np.arange(self.GK_integral.shape[1]) * dt.m_as("ps")).tolist()
        fit = []
        for i in range(3):
            kappa, fit_err, a, a_err, tau, tau_err, tf, yf = fit_green_kubo_integral(
                self.GK_integral[i], gk_ts, sigma=self.GK_integral_err[i]
            )
            fit.append(
                {
                    "kappa": kappa,
                    "stderr": fit_err,
                    "xs": tf,
                    "ys": yf,
                    "a": a,
                    "a_err": a_err,
                    "tau": tau,
                    "tau_err": tau_err,
                }
            )
        self.plot_GK_integral(
            self.GK_integral, gk_ts, err=self.GK_integral_err, fit=fit
        )

        # The values fitted to each run or block
        per_sample, off_diagonal = self._per_sample(ts, gk_ts, starts)
        per_sample = self._pseudo_values(per_sample)

        estimates = {}
        components = {}
        what = "run" if source == "runs" else "block"
        for method in METHODS:
            values = per_sample[method]
            estimate = isotropic_average(values)
            estimate.source = source
            estimate.units = "W/m/K"
            estimates[method] = estimate
            components[method] = [
                Estimate.from_samples(values[:, i], units="W/m/K", source=source)
                for i in range(3)
            ]
            failed = int(np.sum(~np.all(np.isfinite(values), axis=1)))
            if failed == n_samples:
                # The failure of the Helfand moments is explained above
                if method != "Helfand moments" or helfand_fitted:
                    notes.append(f"The {method} fit failed for every {what}.")
            elif failed > 0:
                notes.append(
                    f"The {method} fit failed for {failed} of the {n_samples} "
                    f"{source}, which are left out."
                )

        if self._species is not None:
            notes.append(self._mixture_note(gk_ts, source))
        elif uncorrected is not None:
            notes.append(
                f"The system is a mixture of {len(uncorrected)} species, but the MD "
                "did not sample their mass fluxes, so the energy carried by their "
                "interdiffusion could not be removed: κ is too large, by as much as "
                "40% for 1 m NaCl in water. Run the MD again with a current version "
                "of the LAMMPS step, which samples them with the heat flux."
            )
            printer.important(
                __(
                    "Warning: the thermal conductivity of this mixture is too large: "
                    "the MD did not sample the mass fluxes of the species.",
                    indent=self.indent + 4 * " ",
                )
            )

        isotropic, why = self._isotropy
        table = {"Method": [], "κ": [], "± (95%)": [], "Units": []}

        def row(label, estimate):
            if estimate.n == 0 or not np.isfinite(estimate.value):
                v, e = "--", ""
            else:
                v, e = estimate.formatted()
            table["Method"].append(label)
            table["κ"].append(v)
            table["± (95%)"].append(e)
            table["Units"].append("W/m/K")

        for method in METHODS:
            row(method, estimates[method])
            if not isotropic:
                for name, estimate in zip(DIAGONAL, components[method]):
                    row(f"{method}, {name}", estimate)

        text = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("left", "decimal", "decimal", "left"),
        )
        length = len(text.splitlines()[0])
        printer.normal(
            __(
                "\n" + "Thermal Conductivity".center(length) + "\n" + text,
                indent=8 * " ",
                wrap=False,
                dedent=False,
            )
        )

        # What the table is from, and the checks
        nruns = len(self.paths)
        total = self._total_time.m_as("ps")
        line = (
            f"From {nruns} run{'s' if nruns > 1 else ''} of MD, {total:.0f} ps in "
            f"all, with the heat flux every {dt.m_as('fs'):.3g} fs"
        )
        if source == "blocks":
            line += f", divided into {n_samples} blocks"
        lines = [line + ". " + estimates["Green-Kubo"].describe_error()]
        if isotropic:
            check = component_check(per_sample["Green-Kubo"], names=DIAGONAL)
            line = check.text
            off = self._off_diagonal_check(off_diagonal)
            if off != "":
                line += " " + off
            lines.append(line)
        else:
            lines.append(
                f"The {'runs are' if nruns > 1 else 'run is'} not treated as "
                f"isotropic ({why}), so the diagonal components are given; κ is a "
                "third of their sum."
            )
        for line in lines + notes:
            printer.normal(__(line, indent=8 * " "))
        printer.normal("")

        self._print_samples(per_sample, off_diagonal, source)

        state = report_state(self, printer, self.paths)

        # Store the results
        self.model = self._find_model()
        data = {"nruns": nruns}
        for method in METHODS:
            estimate = estimates[method]
            if estimate.n == 0 or not np.isfinite(estimate.value):
                continue
            key = RESULT_KEYS[method]
            data[key] = estimate.value
            if estimate.stderr is not None:
                data[f"{key},stderr"] = estimate.stderr
        for name, estimate in zip(("xx", "yy", "zz"), components["Green-Kubo"]):
            if estimate.n == 0:
                continue
            data[f"kappa_{name}"] = estimate.value
            if estimate.stderr is not None:
                data[f"kappa_{name},stderr"] = estimate.stderr
        if state is not None:
            data.update(state_data(state))
        self.store_results(configuration=configuration, data=data)

    def _mixture_note(self, gk_ts, source):
        """What the correction for the interdiffusion did, for the notes."""
        values = np.full((len(self.raw_GK_integrals), 3), np.nan)
        for s, integral in enumerate(self.raw_GK_integrals):
            for i in range(3):
                values[s, i] = _fit(fit_green_kubo_integral, integral[i], gk_ts)
        raw = isotropic_average(values)
        raw.source = source
        raw.units = "W/m/K"
        n = len(self._species)
        text = (
            f"The system is a mixture of {n} species, so the energy carried by "
            "their interdiffusion, which is not conduction of heat, is removed from "
            "the heat flux using the mass fluxes of the species (Bertossa et al., "
            "Phys. Rev. Lett. 122, 255901 (2019))."
        )
        if len(self._L_matrices) > 1:
            text += (
                f" The correction uses the averages over the {source}, since that "
                f"for a single one is too noisy, so the errors are from the {source} "
                "by the jackknife (leaving out each in turn), and the values by "
                f"{source[:-1]} in the detailed output are jackknife pseudo-values."
            )
        if raw.n > 0 and np.isfinite(raw.value):
            value, error = raw.formatted()
            if error != "":
                value += f" ± {error}"
            text += f" Without this correction Green-Kubo would give {value} W/m/K."
        return text

    def _uncorrected_mixture(self, configuration, runs):
        """The species of a mixture whose MD has no species fluxes, or None.

        Parameters
        ----------
        configuration : molsystem._Configuration
            The system.
        runs : [tuple]
            The data of the runs.

        Returns
        -------
        [str] or None
        """
        if self._species is not None or configuration is None:
            return None
        try:
            species = sorted(set(configuration.get_molecule_smiles()))
        except Exception:
            return None
        return species if len(species) > 1 else None

    def _per_sample(self, ts, gk_ts, starts):
        """κ fitted to each run or block, by each method.

        Parameters
        ----------
        ts : [float]
            The times of the Helfand moments, in ps.
        gk_ts : [float]
            The times of the Green-Kubo integrals, in ps.
        starts : [float]
            Where to start the fit of each component of the Helfand moments, in ps.

        Returns
        -------
        {str: numpy.ndarray(n, 3)}
            The diagonal components by each method; NaN where a fit failed.
        numpy.ndarray(n, 3)
            The off-diagonal components from the Helfand moments.
        """
        n = len(self.Ms)
        result = {method: np.full((n, 3), np.nan) for method in METHODS}
        off_diagonal = np.full((n, 3), np.nan)
        zeros = np.zeros(len(ts))
        for s, (M, integral) in enumerate(zip(self.Ms, self.GK_integrals)):
            for i in range(3):
                result["Green-Kubo"][s, i] = _fit(
                    fit_green_kubo_integral, integral[i], gk_ts
                )
                result["Helfand moments"][s, i] = _fit(
                    fit_helfand_moment, M[i], ts, start=starts[i]
                )
                slope, x, _ = helfand_slope(M[i], ts, zeros)
                result["Helfand derivative"][s, i] = _fit(
                    fit_green_kubo_integral, slope, x
                )
            for i in range(3, 6):
                off_diagonal[s, i - 3] = _fit(
                    fit_helfand_moment, M[i], ts, start=starts[i]
                )
        return result, off_diagonal

    @staticmethod
    def _off_diagonal_check(off_diagonal):
        """One sentence: whether the off-diagonal terms are zero within their
        errors, as they should be for an isotropic system; empty if it cannot be
        told.

        Each term is tested against zero with Student's t (n - 1 degrees of freedom,
        since the errors come from only n runs or blocks), and the smallest p-value
        is corrected for the three tests (Bonferroni).
        """
        values = np.asarray(off_diagonal, dtype=float)
        values = values[np.all(np.isfinite(values), axis=1)]
        n = len(values)
        if n < 2:
            return ""
        means, stderrs = mean_and_stderr(values)
        if not np.all(stderrs > 0):
            return ""
        from scipy import stats

        t = np.abs(means / stderrs)
        worst = int(np.argmax(t))
        p = min(1.0, len(t) * float(2 * stats.t.sf(t[worst], n - 1)))
        name = OFF_DIAGONAL[worst]
        how_far = f"{t[worst]:.1f} standard errors from zero; p = {p:.2g}"
        if p >= 0.01:
            return (
                "The off-diagonal terms are zero within their errors (the largest, "
                f"{name}, is {how_far})."
            )
        return (
            f"The off-diagonal terms are not zero within their errors ({name} is "
            f"{how_far}): the system may not be isotropic, or the runs too short."
        )

    def _print_samples(self, per_sample, off_diagonal, source):
        """The values fitted to each run or block, in the detailed output."""
        what = "Run" if source == "runs" else "Block"
        diagonal = np.concatenate([per_sample[m] for m in METHODS])
        largest = np.nanmax(np.abs(diagonal)) if np.any(np.isfinite(diagonal)) else 1
        decimals = max(0, 3 - floor(log10(largest))) if largest > 0 else 3

        def fmt(value):
            return f"{value:.{decimals}f}" if np.isfinite(value) else "--"

        for method in METHODS:
            values = per_sample[method]
            table = {what: [], "κ": []}
            for name in DIAGONAL:
                table[name] = []
            if method == "Helfand moments":
                for name in OFF_DIAGONAL:
                    table[name] = []
            for s, components in enumerate(values):
                table[what].append(str(s + 1))
                table["κ"].append(fmt(components.mean()))
                for name, value in zip(DIAGONAL, components):
                    table[name].append(fmt(value))
                if method == "Helfand moments":
                    for name, value in zip(OFF_DIAGONAL, off_diagonal[s]):
                        table[name].append(fmt(value))
            text = tabulate(
                table,
                headers="keys",
                tablefmt="rounded_outline",
                disable_numparse=True,
                colalign=("center",) + ("decimal",) * (len(table) - 1),
            )
            length = len(text.splitlines()[0])
            title = f"{method} by {what.lower()} (W/m/K)".center(length)
            print_details(
                printer,
                str(__("\n" + title + "\n" + text + "\n", indent=8 * " ", wrap=False)),
            )

    def _find_model(self):
        """The model (forcefield or model chemistry) of the MD, for properties."""
        try:
            if self.variable_exists("_model_chemistry"):
                return self.get_variable("_model_chemistry")["level"]
            ff = self.get_variable("_forcefield")
            return ff.current_forcefield
        except Exception:
            return "unknown"

    def _plot_helfand_derivatives(self, slopes, xs, errs, fit):
        """Plot the derivatives of the Helfand moments and their fits."""
        figure = self.create_figure(
            module_path=("seamm",),
            template="line.graph_template",
            title="Helfand Derivatives",
        )
        plot_helfand_slopes(
            figure, slopes, xs[0], err=errs, fit=fit, labels=self.tensor_labels
        )
        figure.grid_plots("Slope")

        path = Path(self.directory) / "HelfandDerivatives.graph"
        figure.dump(path)
        if "html" in self.options and self.options["html"]:
            path = path.with_suffix(".html")
            figure.template = "line.html_template"
            figure.dump(path)

    def plot_GK_integral(self, x, ts, err=None, fit=None, filename="GK_integral.graph"):
        """Create a plot of the GK integrals and any fit to them.

        Parameters
        ----------
        x : numpy.ndarray(6, m)
            The cumulative integrals

        ts : [float]
            List of times for the points, in ps

        err : numpy.ndarray(6, m)
            The standard errors of the integrals (optional)

        fit :
            The fit parameters for the integrals (optional)

        filename : str
            The filename to write
        """
        figure = self.create_figure(
            module_path=("seamm",),
            template="line.graph_template",
            title="Green-Kubo Integrals",
        )

        max_tau = 0
        for i in range(3):
            tau = fit[i]["tau"][-1]
            if tau > max_tau:
                max_tau = tau

        if max_tau * 4 > ts[-1]:
            rng = None
        else:
            rng = [0, max_tau * 4]
        plot_GK_integrals(
            figure, x, ts, err=err, fit=fit, labels=self.tensor_labels, _range=rng
        )

        figure.grid_plots("GKI")

        # Write to disk
        path = Path(self.directory) / filename
        figure.dump(path)

        if "html" in self.options and self.options["html"]:
            path = path.with_suffix(".html")
            figure.template = "line.html_template"
            figure.dump(path)

    def plot_helfand_moment(
        self, M, ts, M_err=None, fit=None, filename="HelfandMoments.graph"
    ):
        """Create a plot of the Helfand moments and any fit to them.

        Parameters
        ----------
        M : numpy.ndarray(6, m)
            The Helfand moments.

        ts : [float]
            List of times for the points, in ps
        """
        figure = self.create_figure(
            module_path=("seamm",),
            template="line.graph_template",
            title="Helfand Moments",
        )

        plot_helfand_moments(
            figure, M, ts, err=M_err, fit=fit, labels=self.tensor_labels
        )

        figure.grid_plots("M")

        # Write to disk
        path = Path(self.directory) / filename
        figure.dump(path)

        if "html" in self.options and self.options["html"]:
            path = path.with_suffix(".html")
            figure.template = "line.html_template"
            figure.dump(path)

    def plot_JJ_correlation(self, CCF, ts, err=None, fit=None, filename="JJ.graph"):
        """Create a plot of the correlation functions and any fit to them.

        Parameters
        ----------
        CCF : numpy.ndarray(6, m)
            The cross-correlation functions.

        ts : [float]
            List of times for the points, in ps

        err : numpy.ndarray(6, m)
            The standard errors of the CCF (optional)

        fit :
            The fit parameters for any fit of the CFF (optional)

        filename : str
            The filename to write
        """
        figure = self.create_figure(
            module_path=("seamm",),
            template="line.graph_template",
            title="Heat Flux Correlation Functions",
        )

        plot_correlation_functions(
            figure, CCF, ts, err=err, fit=fit, labels=self.tensor_labels
        )

        figure.grid_plots("CCF")

        # Write to disk
        path = Path(self.directory) / filename
        figure.dump(path)

        if "html" in self.options and self.options["html"]:
            path = path.with_suffix(".html")
            figure.template = "line.html_template"
            figure.dump(path)

    def process_run(self, run, data, what="Run"):
        """Calculate the Helfand moments and correlation functions of one run or
        block.

        Parameters
        ----------
        run : int
            The number of the run or block, for messages.
        data : (numpy.ndarray, pint.Quantity, float, float, numpy.ndarray)
            The heat flux (3 × n, W/m^2), the sampling interval, the temperature
            (K), the volume (Å^3), and for a mixture the mass fluxes of all the
            species but one (k × 3 × n), or None.
        what : str
            "Run" or "Block", for messages.
        """
        J, dt, T, V, S = data
        self.V.append(V)
        self.T.append(T)
        self.timestep.append(dt.m_as("fs"))

        T = Q_(T, "K")
        V = Q_(V, "Å^3")
        k_B = Q_("k_B")
        Jsq = Q_("W^2/m^4")

        # Limit the lengths of the data
        n = J.shape[1]
        m = min(n // 10, 10000)
        if n * dt < MIN_RUN_TIME:
            raise RuntimeError(
                f"{what} {run} is only {(n * dt).to('ps'):~P} long, but at least "
                f"{MIN_RUN_TIME:~P} is needed for the correlation functions to be "
                "fitted."
            )

        # The Helfand moments
        constants = (Jsq * V * dt**2 / (2 * k_B * T**2)).m_as("W/m/K*ps")
        M = collective_helfand_moments(J, m=m) * constants
        self.Ms.append(M)
        if S is not None:
            # For a mixture, also the matrices with the species' mass fluxes, from
            # which _combine removes their interdiffusion
            self._M_matrices.append(flux_helfand_matrices(J, S, m=m) * constants)

        # And correlation functions for Green-Kubo method
        constants = (V / (k_B * T**2) * Jsq * dt).m_as("W/m/K")
        Jcf, integral = tensor_correlation_functions(J, m=m)
        self.Jcfs.append(Jcf)
        self.GK_integrals.append(integral * constants)
        if S is not None:
            self.raw_GK_integrals.append(integral * constants)
            C, L = flux_correlation_matrices(J, S, m=m)
            self._C_matrices.append(C)
            self._L_matrices.append(L * constants)

    def _combine(self):
        """The averages of the curves over the runs or blocks, and their spread.

        For a mixture the thermal conductivity is the Schur complement of the matrix
        of the heat flux and the species' mass fluxes. The complement of a single
        run or block is very noisy -- the species' part is small and poorly known --
        so the matrices are averaged first. The diagonal curves of each run or block
        are then those of the others (leave-one-out, or jackknife), whose fits
        :meth:`analyze` turns into pseudo-values; their spread is scaled to match.
        """
        if len(self._L_matrices) > 0:
            self._remove_interdiffusion()
        for name in ("M", "Jcf", "GK_integral"):
            curves = np.stack(getattr(self, name + "s"))
            setattr(self, name, np.average(curves, axis=0))
            setattr(self, name + "_err", np.std(curves, axis=0))
        if len(self._L_matrices) > 0:
            # The averages are the complements of the averaged matrices, and the
            # spread of the leave-one-out curves is (n - 1) times too small
            n = len(self._L_matrices)
            C = np.mean(self._C_matrices, axis=0)
            L = np.mean(self._L_matrices, axis=0)
            self.M[:3] = schur_complement(np.mean(self._M_matrices, axis=0))
            self.GK_integral[:3] = schur_complement(L)
            self.Jcf[:3] = effective_correlation_functions(C, L)
            for name in ("M", "Jcf", "GK_integral"):
                getattr(self, name + "_err")[:3] *= max(1, n - 1)

    def _remove_interdiffusion(self):
        """Replace the diagonal curves of each run or block of a mixture by those
        without the interdiffusion of the species: the Schur complements of the
        matrices averaged over the other runs or blocks, or of the one there is."""
        n = len(self._L_matrices)
        averages = {}
        for name in ("_C_matrices", "_L_matrices", "_M_matrices"):
            averages[name] = np.mean(getattr(self, name), axis=0)

        def others(name, i):
            if n == 1:
                return averages[name]
            return (n * averages[name] - getattr(self, name)[i]) / (n - 1)

        for i in range(n):
            L = others("_L_matrices", i)
            self.GK_integrals[i][:3] = schur_complement(L)
            self.Jcfs[i][:3] = effective_correlation_functions(
                others("_C_matrices", i), L
            )
            self.Ms[i][:3] = schur_complement(others("_M_matrices", i))

    def _pseudo_values(self, per_sample):
        """The jackknife pseudo-values from the fits to the leave-one-out curves of
        a mixture: mean - (n - 1) (value - mean), whose mean and standard error are
        the jackknife estimate and its error. Unchanged for a pure substance or a
        single run or block."""
        n = len(self._L_matrices)
        if n < 2:
            return per_sample
        result = {}
        for method, values in per_sample.items():
            mean = np.nanmean(values, axis=0) if np.any(np.isfinite(values)) else 0.0
            result[method] = mean - (n - 1) * (values - mean)
        return result

    def md_requirements(self, P=None):
        """What the analysis needs from the MD: the heat flux every two steps, in a
        periodic cell. It analyzes the dynamics."""
        return Requirements(
            needs=[
                Need(
                    "heat_flux",
                    max_steps=2,
                    why="the thermal conductivity, to follow its fastest oscillations",
                )
            ],
            periodic=True,
            dynamics=True,
        )

    def _static(self, key):
        """A parameter's value without evaluating variables, or None."""
        parameter = self.parameters[key]
        if parameter.is_expr:
            return None
        return parameter.get()

    def check_flowchart(self, preceding):
        """Check, before the job runs, that the MD will sample what is needed and is
        set up for the analysis.

        Parameters
        ----------
        preceding : [seamm.Node]
            The steps that run before this one.

        Returns
        -------
        [str]
            A description of each problem; empty if none.
        """
        source = self._static("md source")
        if source is None:
            return []
        return check_md_flowchart(self, preceding, self.md_requirements(), source)

    def check_configuration(self, configuration):
        """Called by an MD step before it runs MD that this step will analyze.

        Returns
        -------
        [str]
            Warnings about the structure, e.g. that it is small.

        Raises
        ------
        RuntimeError
            If the system is not periodic.
        """
        if configuration.periodicity != 3:
            raise RuntimeError(
                f"'{self.title}' will analyze the MD, but the system is not "
                "periodic: the heat flux is divided by the volume of the cell."
            )
        return [f"'{self.title}': {w}" for w in size_warnings(configuration)]

    def _read_md_file(self, path):
        """The heat flux, interval, temperature, volume and species fluxes of one
        md.h5md, noting its length, its species and whether it is isotropic."""
        result = self._read_h5md(path)
        if result is None:
            raise RuntimeError(f"There is no heat flux in {path}.")
        J, dt, _, _, S, species = result
        result = result[:5]
        # Every run must be corrected for the interdiffusion of the species, or none
        first = self._total_time.m == 0
        if first:
            self._species = species
        elif (species is None) != (self._species is None):
            raise RuntimeError(
                f"Some of the runs sampled the mass fluxes of the species of the "
                f"mixture, but not all of them ({path} "
                f"{'did not' if species is None else 'did'}). Analyze runs made the "
                "same way."
            )
        elif species != self._species:
            raise RuntimeError(
                f"The species of the MD in {path}, {', '.join(species)}, are not those "
                f"of the other runs, {', '.join(self._species)}."
            )
        self._total_time += J.shape[1] * dt
        with open_md(path) as md:
            isotropic, why = isotropy(md)
        # The runs are isotropic only if all are; otherwise say why not
        if self._isotropy is None or (self._isotropy[0] and not isotropic):
            self._isotropy = (isotropic, why)
        return result

    @staticmethod
    def _split_into_blocks(data, n_blocks):
        """Divide one run into consecutive blocks, analyzed as separate runs."""
        J, dt, T, V, S = data
        n_blocks = max(1, int(n_blocks))
        n = J.shape[1] // n_blocks
        if n * dt < MIN_RUN_TIME:
            raise RuntimeError(
                f"The run is only {(J.shape[1] * dt).to('ps'):~P} long, too short for "
                f"{n_blocks} blocks: each block needs at least {MIN_RUN_TIME:~P} for "
                "the correlation functions to be fitted. Use fewer blocks or a longer "
                "run."
            )
        if n < 2000:
            raise RuntimeError(
                f"The run has only {J.shape[1]} samples of the heat flux, too few for "
                f"{n_blocks} blocks: each block needs at least 2000 for the "
                "correlation functions to be fitted. Use fewer blocks or a longer run."
            )
        return [
            (
                J[:, i * n : (i + 1) * n],
                dt,
                T,
                V,
                None if S is None else S[..., i * n : (i + 1) * n],
            )
            for i in range(n_blocks)
        ]

    @staticmethod
    def _read_h5md(path):
        """The heat flux, temperature, volume and species fluxes in an md.h5md file.

        Returns
        -------
        (numpy.ndarray(3, n), pint.Quantity, float, float, numpy.ndarray, [str])
            The heat flux in W/m^2, the interval between samples, the mean
            temperature in K, the mean volume in Å^3, and for a mixture the mass
            fluxes of all the species but the first (k × 3 × n, in any units) and
            the species, or None and None; or None if the file has no heat flux.
        """
        with open_md(path) as md:
            if not md.has_observable("heat_flux"):
                return None
            if not md.periodic:
                raise RuntimeError(
                    f"The MD in {path} is not periodic, so it has no thermal "
                    "conductivity."
                )
            flux = md.observable("heat_flux", units="W/m^2", time_units="fs")
            T = float(md.observable("temperature", units="K").values.mean())
            if md.has_observable("volume"):
                V = float(md.observable("volume", units="Å^3").values.mean())
            else:
                V = float(abs(np.linalg.det(md.cell(units="Å"))))
            J = flux.values.T
            S = None
            species = None
            if md.has_observable("species_mass_flux"):
                fluxes = md.observable("species_mass_flux")
                species = list(fluxes.attrs.get("species", []))
                if fluxes.values.shape[1] > 1:
                    # The first species, with the most atoms, is fixed by the others
                    # and the total momentum
                    _, i, j = np.intersect1d(
                        flux.steps, fluxes.steps, return_indices=True
                    )
                    J = J[:, i]
                    S = fluxes.values[j, 1:, :].transpose(1, 2, 0)
                else:
                    species = None
        return J, flux.dt, T, V, S, species
