# -*- coding: utf-8 -*-

"""Non-graphical part of the Thermal Conductivity step in a SEAMM flowchart
"""

import logging
from math import log10, ceil
from pathlib import Path
import pkg_resources
import sys
import traceback

import numpy as np
from tabulate import tabulate

from .analysis import (
    create_correlation_functions,
    create_helfand_moments,
    fit_green_kubo_integral,
    fit_helfand_moment,
    get_helfand_slope,
    plot_correlation_functions,
    plot_GK_integrals,
    plot_helfand_moments,
    plot_helfand_slopes,
)
from .thermal_conductivity_parameters import ThermalConductivityParameters
import thermal_conductivity_step
import molsystem
import seamm
import seamm_util
from seamm_util import ureg, Q_  # noqa: F401
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
job = printing.getPrinter()
printer = printing.getPrinter("Thermal Conductivity")

# Add this module's properties to the standard properties
path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


def fmt_err(value, err, precision=2):
    decimals = -ceil(log10(err)) + precision
    if decimals < 0:
        decimals = 0
    fmt = f".{decimals}f"
    e = f"{err:{fmt}}"
    v = f"{value:{fmt}}"
    return v, e


class ThermalConductivity(seamm.Node):
    """
    The non-graphical part of a Thermal Conductivity step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

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
        namespace="org.molssi.seamm",
        extension=None,
        logger=logger,
    ):
        """A step for Thermal Conductivity in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        namespace : str
            The namespace for the plug-ins of the subflowchart
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Thermal Conductivity {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Thermal Conductivity", namespace=namespace
        )

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
        self._file_handler = None

        # Various data
        self.J = None  # The heat fluxes, nruns x 3 x n
        self.V = []  # The volumes of the runs
        self.T = []  # The temperatures of the runs
        self.timestep = []  # The timestep in fs for the run
        self.Ms = []  # Helfand moments per run
        self.M = None  # The Helfand moments, 6 x m
        self.M_err = None  # The stderr on HM, 6 x m
        self.Jcfs = []  # The ccf's of J, 6 x n, per run
        self.Jcf = None  # The ccf's of J, 6 x n
        self.Jcf_err = None  # The stderr of ccf's, 6 x n
        self.GK_integrals = []  # Cumulative integral of Jcf, per run
        self.GK_integral = None  # Cumulative integral of Jcf
        self.GK_integral_err = None  # Error of cumulative integral of Jcf

        # Set our citation level to 1!
        self.citation_level = 1

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

        # Now need to walk through the steps in the subflowchart...
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            node = node.create_parser()

        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        Returns
        -------
        str
            A description of the current step.
        """
        if P is None:
            P = self.parameters.values_to_dict()

        text = (
            f"Calculate the thermal conductivity at {P['T']} using "
            f"the {P['approach']} approach. {P['nruns']} runs will be averaged for "
            "the final results.\n\n"
        )
        if not short:
            # The subflowchart
            self.subflowchart.root_directory = self.flowchart.root_directory

            # Get the first real node
            node = self.subflowchart.get_node("1").next()

            while node is not None:
                try:
                    text += __(node.description_text()).__str__()
                except Exception as e:
                    print(
                        f"Error describing thermal_conductivity flowchart: {e} in "
                        f"{node}"
                    )
                    logger.critical(
                        f"Error describing thermal_conductivity flowchart: {e} in "
                        f"{node}"
                    )
                    raise
                except:  # noqa: E722
                    print(
                        "Unexpected error describing thermal_conductivity flowchart: "
                        f"{sys.exc_info()[0]} in {str(node)}"
                    )
                    logger.critical(
                        "Unexpected error describing thermal_conductivity flowchart: "
                        f"{sys.exc_info()[0]} in {str(node)}"
                    )
                    raise
                text += "\n"
                node = node.next()

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def run(self):
        """Run a Thermal Conductivity step.

        Parameters
        ----------
        None

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
        printer.important(__(self.description_text(P, short=True), indent=self.indent))

        # Find the handler for job.out and set the level up
        job_handler = None
        out_handler = None
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ):
                job_handler = handler
                job_level = job_handler.level
                job_handler.setLevel(printing.JOB)
            elif isinstance(handler, logging.StreamHandler):
                out_handler = handler
                out_level = out_handler.level
                out_handler.setLevel(printing.JOB)

        # Get the first real node
        first_node = self.subflowchart.get_node("1").next()

        # Ensure the nodes have their options
        node = first_node
        while node is not None:
            node.all_options = self.all_options
            node = node.next()

        # Loop over the runs
        nruns = P["nruns"]
        fmt = f"0{len(str(nruns))}"
        for run in range(1, nruns + 1):
            # Direct most output to iteration.out
            run_dir = Path(self.directory) / f"run_{run:{fmt}}"
            run_dir.mkdir(parents=True, exist_ok=True)

            # A handler for the file
            if self._file_handler is not None:
                self._file_handler.close()
                job.removeHandler(self._file_handler)
            path = run_dir / "Run.out"
            path.unlink(missing_ok=True)
            self._file_handler = logging.FileHandler(path)
            self._file_handler.setLevel(printing.NORMAL)
            formatter = logging.Formatter(fmt="{message:s}", style="{")
            self._file_handler.setFormatter(formatter)
            job.addHandler(self._file_handler)

            # Add the run to the ids so the directory structure is reasonable
            self.flowchart.reset_visited()
            self.set_subids((*self._id, f"run_{run:{fmt}}"))

            # Run through the steps in the loop body
            node = first_node
            try:
                while node is not None:
                    node = node.run()
            except DeprecationWarning as e:
                printer.normal("\nDeprecation warning: " + str(e))
                traceback.print_exc(file=sys.stderr)
                traceback.print_exc(file=sys.stdout)
            except Exception as e:
                printer.job(f"Caught exception in run {run}: {str(e)}")
                with open(run_dir / "stderr.out", "a") as fd:
                    traceback.print_exc(file=fd)
                if "continue" in P["errors"]:
                    continue
                elif "exit" in P["errors"]:
                    break
                else:
                    raise
            else:
                self.process_run(run, run_dir)
                if True and len(self.V) > 1:
                    if job_handler is not None:
                        job_handler.setLevel(job_level)
                    if out_handler is not None:
                        out_handler.setLevel(out_level)

                    self.analyze(P=P, style="1-line", run=run)

                    if job_handler is not None:
                        job_handler.setLevel(printing.JOB)
                    if out_handler is not None:
                        out_handler.setLevel(printing.JOB)

            self.logger.debug(f"End of run {run}")

        # Remove any redirection of printing.
        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
            self._file_handler = None
        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)

        # Analyze the results
        self.analyze(P=P, style="full")

        # Add other citations here or in the appropriate place in the code.
        # Add the bibtex to data/references.bib, and add a self.reference.cite
        # similar to the above to actually add the citation to the references.

        return next_node

    def analyze(self, indent="", P=None, style="full", run=None, **kwargs):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        """
        if style == "1-line":
            table = {
                "Run": [],
                "Kxx": [],
                "exx": [],
                "Kyy": [],
                "eyy": [],
                "Kzz": [],
                "ezz": [],
                "Kxy": [],
                "exy": [],
                "Kxz": [],
                "exz": [],
                "Kyz": [],
                "eyz": [],
            }
        elif style == "full":
            table = {
                "Method": [],
                "Dir": [],
                "Kappa": [],
                "±": [],
                "95%": [],
                "Units": [],
            }

        dt = Q_(self.timestep[0], "fs")
        ts = np.arange(self.M.shape[1]) * dt.m_as("ps")  # Scale to ps
        ts = ts.tolist()

        # Get and plot the slopes of the Helfand moments
        slopes = []
        xs = []
        errs = []
        fit0 = []
        for i in range(6):
            slope, x, err = get_helfand_slope(self.M[i], ts, sigma=self.M_err[i])
            slopes.append(slope)
            xs.append(x)
            errs.append(err)

            if i < 3:
                (
                    kappa,
                    kappa_err,
                    a,
                    a_err,
                    tau,
                    tau_err,
                    tf,
                    yf,
                ) = fit_green_kubo_integral(slope, x, sigma=err)
                # print(f"{kappa=} {tau=}")
                fit0.append(
                    {
                        "kappa": kappa,
                        "stderr": kappa_err,
                        "xs": tf,
                        "ys": yf,
                        "a": a,
                        "a_err": a_err,
                        "tau": tau,
                        "tau_err": tau_err,
                    }
                )

        figure = self.create_figure(
            module_path=("seamm",),
            template="line.graph_template",
            title="Helfand Slopes",
        )

        plot_helfand_slopes(
            figure, slopes, xs[0], err=errs, fit=fit0, labels=self.tensor_labels
        )

        figure.grid_plots("Slope")

        # Write to disk
        filename = "Helfand_slopes.graph"
        path = Path(self.directory) / filename
        figure.dump(path)

        if "html" in self.options and self.options["html"]:
            path = path.with_suffix(".html")
            figure.template = "line.html_template"
            figure.dump(path)

        # Fit the slopes
        fit = []
        for i in range(6):
            if i < len(fit0):
                start = max(fit0[i]["tau"]) * 3
            else:
                start = 1
            slope, err, xs, ys = fit_helfand_moment(
                self.M[i], ts, sigma=self.M_err[i], start=start
            )
            fit.append(
                {
                    "kappa": slope,
                    "stderr": err,
                    "xs": xs,
                    "ys": ys,
                }
            )
            v, e = fmt_err(slope, 2 * err)
            if style == "1-line":
                if i == 0:
                    table["Run"].append(len(self.V))
                alpha = self.tensor_labels[i][0]
                table["K" + alpha].append(v)
                table["e" + alpha].append(e)
            elif style == "full":
                table["Method"].append("Helfand" if i == 0 else "")
                table["Dir"].append(self.tensor_labels[i][0])
                table["Kappa"].append(v)
                table["±"].append("±")
                table["95%"].append(e)
                table["Units"].append("W/m/K" if i == 0 else "")

        self.plot_helfand_moment(self.M, ts, M_err=self.M_err, fit=fit)

        # Now for the Green-Kubo method

        dt = Q_(self.timestep[0], "fs")

        # Plot the correlation functions
        ts = (np.arange(self.Jcf.shape[1]) * dt.m_as("ps")).tolist()
        end = len(ts) // 4
        self.plot_JJ_correlation(self.Jcf[:, :end], ts[:end], err=self.Jcf_err[:, :end])

        # Fit the integrals to 1 - e(-tau/t0)
        ts = np.arange(self.GK_integral.shape[1]) * dt.m_as("ps")
        ts = ts.tolist()

        fit = []
        for i in range(3):
            kappa, kappa_err, a, a_err, tau, tau_err, tf, yf = fit_green_kubo_integral(
                self.GK_integral[i], ts, sigma=self.GK_integral_err[i]
            )
            fit.append(
                {
                    "kappa": kappa,
                    "stderr": kappa_err,
                    "xs": tf,
                    "ys": yf,
                    "a": a,
                    "a_err": a_err,
                    "tau": tau,
                    "tau_err": tau_err,
                }
            )
            try:
                v, e = fmt_err(kappa, 2 * kappa_err)
            except Exception:
                v = f"{kappa:.2f}"
                e = "--"

            if style == "1-line":
                if i == 0:
                    table["Run"].append("")
                alpha = self.tensor_labels[i][0]
                table["K" + alpha].append(v)
                table["e" + alpha].append(e)
            elif style == "full":
                table["Method"].append("Green-Kubo" if i == 0 else "")
                table["Dir"].append(self.tensor_labels[i][0])
                table["Kappa"].append(v)
                table["±"].append("±")
                table["95%"].append(e)
                table["Units"].append("W/m/K" if i == 0 else "")

        for i in range(3, 6):
            alpha = self.tensor_labels[i][0]
            if style == "1-line":
                # blank the off-diagonals
                table["K" + alpha].append("")
                table["e" + alpha].append("")

        # Plot the Green-Kubo integrals
        self.plot_GK_integral(self.GK_integral, ts, err=self.GK_integral_err, fit=fit)

        # Print the table of results
        if style == "1-line":
            text = ""
            tmp = tabulate(
                table,
                headers="keys",
                tablefmt="simple_outline",
                disable_numparse=True,
            )
            if len(self.V) == 2:
                length = len(tmp.splitlines()[0])
                text += "\n"
                text += "Thermal Conductivity".center(length)
                text += "\n"
                text += "--------------------".center(length)
                text += "\n"
                text += "First line is Helfand moments; second, Green-Kubo".center(
                    length
                )
                text += "\n"
                text += "\n".join(tmp.splitlines()[0:-1])
            else:
                text = tmp.splitlines()[-3]
                text += "\n"
                text += tmp.splitlines()[-2]
                if run is not None and run == P["nruns"]:
                    text += "\n"
                    text += tmp.splitlines()[-1]

            printer.normal(__(text, indent=8 * " ", wrap=False, dedent=False))
        else:
            text = ""
            tmp = tabulate(
                table,
                headers="keys",
                tablefmt="simple_outline",
                disable_numparse=True,
                colalign=(
                    "center",
                    "center",
                    "decimal",
                    "center",
                    "decimal",
                    "left",
                ),
            )
            length = len(tmp.splitlines()[0])
            text += "\n"
            text += "Thermal Conductivity".center(length)
            text += "\n"
            text += tmp
            text += "\n"

            printer.normal(__(text, indent=8 * " ", wrap=False, dedent=False))

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
        """Create a plot of the orrelation functions and any fit to them.

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

    def process_run(self, run, run_dir):
        """Get the heat fluxes from the run and do initial processing.

        Parameters
        ----------
        run : int
            The run number
        run_dir : pathlib.Path
            The toplevel directory of the run.
        """
        paths = sorted(run_dir.glob("**/HeatFlux.npz"))

        if len(paths) == 0:
            raise RuntimeError(f"There is no heat flux data for run {run}.")
        elif len(paths) > 1:
            raise NotImplementedError(
                f"Cannot handle multiple HeatFlux files from run {run}."
            )
        with np.load(paths[0]) as data:
            V = float(data["V"])
            T = float(data["T"])
            timestep = float(data["timestep"])

            J = data["J"]

        self.V.append(V)
        self.T.append(T)
        self.timestep.append(timestep)

        # We need properties like the temperature and volume.
        T = Q_(T, "K")
        V = Q_(V, "Å^3")
        dt = Q_(timestep, "fs")
        k_B = Q_("k_B")
        Jsq = Q_("W^2/m^4")

        # Limit the lengths of the data
        n = J.shape[1]
        m = min(n // 20, 10000)

        # Create the Helfand moments
        constants = Jsq * V * dt**2 / (2 * k_B * T**2)
        M = create_helfand_moments(J, m=m) * constants.m_as("W/m/K*ps")
        self.Ms.append(M)

        # And correlation functions for Green-Kubo method
        constants = V / (k_B * T**2) * Jsq * dt
        Jcf, integral = create_correlation_functions(J, m=m)
        self.Jcfs.append(Jcf)
        self.GK_integrals.append(integral * constants.m_as("W/m/K"))
        # Merge the results and get the averages and errors
        tmp = np.stack(self.Ms)
        self.M = np.average(tmp, axis=0)
        self.M_err = np.std(tmp, axis=0)
        tmp = np.stack(self.Jcfs)
        self.Jcf = np.average(tmp, axis=0)
        self.Jcf_err = np.std(tmp, axis=0)
        tmp = np.stack(self.GK_integrals)
        self.GK_integral = np.average(tmp, axis=0)
        self.GK_integral_err = np.std(tmp, axis=0)

    def set_id(self, node_id=()):
        """Sequentially number the subnodes"""
        self.logger.debug("Setting ids for subflowchart {}".format(self))
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            self.set_subids(self._id)
            return self.next()

    def set_subids(self, node_id=()):
        """Set the ids of the nodes in the subflowchart"""
        node = self.subflowchart.get_node("1").next()
        n = 1
        while node is not None:
            node = node.set_id((*node_id, str(n)))
            n += 1
