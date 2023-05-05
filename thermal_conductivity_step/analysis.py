# -*- coding: utf-8 -*-

"""Routines to help do Green-Kubo and Helfand moments analysis."""
import warnings

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import curve_fit, OptimizeWarning
import statsmodels.tsa.stattools as stattools

tensor_labels = [
    ("xx", "red", "rgba(255,0,0,0.1)"),
    ("yy", "green", "rgba(0,255,0,0.1)"),
    ("zz", "blue", "rgba(0,0,255,0.1)"),
    ("xy", "rgb(127,127,0)", "rgba(127,127,0,0.1)"),
    ("xz", "rgb(255,0,255)", "rgba(255,0,255,0.1)"),
    ("yz", "rgb(0,255,255)", "rgba(0,255,255,0.1)"),
]
a1 = 0
tau1 = 0


def moving_average(a, n):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1 :] / n


def frequencies(data):
    ps = np.abs(np.fft.rfft(data)) ** 2

    # time_step = 1 / 100
    freqs = np.fft.rfftfreq(data.size)
    idx = np.argsort(freqs)

    for i in idx[0:500]:
        print(f"{freqs[i]:.4f} {ps[i]:12.0f}")


def exp_1(t, a, tau):
    return a * (1 - np.exp(-t / tau))


def exp_1a(t, a2, tau2):
    global a1, tau1
    return a1 * (1 - np.exp(-t / tau1)) + a2 * (1 - np.exp(-t / tau2))


def exp_2(t, a1, tau1, a2, tau2):
    return a1 * (1 - np.exp(-t / tau1)) + a2 * (1 - np.exp(-t / tau2))


def axb(x, a, b):
    return a * x + b


def fit_h(x, a, b, tau):
    return a * (tau * np.exp(-x / tau) + x) + b


def acf_err(acf):
    """The standard error of the autocorrelation function.

    Copied from statsmodels.tsa.stattools.acf
    """
    nobs = acf.shape[0]
    varacf = np.ones_like(acf) / nobs
    varacf[0] = 0
    varacf[1] = 1.0 / nobs
    varacf[2:] *= 1 + 2 * np.cumsum(acf[1:-1] ** 2)
    std = np.sqrt(varacf)
    return std


def ccf_err(acf1, acf2):
    """The standard error of crosscorrelation functions.

    Copied from statsmodels.tsa.stattools.acf
    """
    nobs = acf1.shape[0]
    varccf = np.ones_like(acf1) / nobs
    varccf[0] = 0
    varccf[1] = 1.0 / nobs
    varccf[2:] *= 1 + 2 * np.cumsum(acf1[1:-1] * acf2[1:-1])
    std = np.sqrt(varccf)
    return std


def create_correlation_functions(J, m=None):
    """Create the correlation functions of the heat fluxes.

    Parameters
    ----------
    J : numpy.ndarray(3, n)
        The heat fluxes in x, y, and z

    Returns
    -------
    numpy.ndarray(6, n)
        The correlation functions
    """

    if m is None:
        m = J.shape[1]

    ccf = np.zeros((6, m))
    integral = np.zeros((6, m))

    ccf[0] = stattools.ccovf(J[0], J[0])[:m]
    ccf[1] = stattools.ccovf(J[1], J[1])[:m]
    ccf[2] = stattools.ccovf(J[2], J[2])[:m]
    ccf[3] = stattools.ccovf(J[0], J[1])[:m]
    ccf[4] = stattools.ccovf(J[0], J[2])[:m]
    ccf[5] = stattools.ccovf(J[1], J[2])[:m]

    integral[0] = cumulative_trapezoid(ccf[0], initial=0.0)
    integral[1] = cumulative_trapezoid(ccf[1], initial=0.0)
    integral[2] = cumulative_trapezoid(ccf[2], initial=0.0)
    integral[3] = cumulative_trapezoid(ccf[3], initial=0.0)
    integral[4] = cumulative_trapezoid(ccf[4], initial=0.0)
    integral[5] = cumulative_trapezoid(ccf[5], initial=0.0)

    return ccf, integral


def create_helfand_moments(J, m=None):
    """Create the Helfand moments from heat fluxes.

    Parameters
    ----------
    J : numpy.ndarray(3, n)
        The heat fluxes in x, y, and z
    m : int
        The length of the Helfand moments wanted

    Returns
    -------
    numpy.ndarray(6, m)
        The Helfand moments
    """

    n = J.shape[1]
    if m is None:
        m = min(n // 20, 10000)

    M = np.zeros((6, m))
    for i in range(n - m):
        Ix = cumulative_trapezoid(J[0, i : m + i], initial=0.0)
        Iy = cumulative_trapezoid(J[1, i : m + i], initial=0.0)
        Iz = cumulative_trapezoid(J[2, i : m + i], initial=0.0)

        M[0, :] += Ix * Ix
        M[1, :] += Iy * Iy
        M[2, :] += Iz * Iz
        M[3, :] += Ix * Iy
        M[4, :] += Ix * Iz
        M[5, :] += Iy * Iz

    M /= n - m

    return M


def fit_green_kubo_integral(y, xs, sigma=None):
    """Find the best a * (1 - exp(-tau/t)) form.

    Parameters
    ----------
    y : [float] or numpy.ndarray()
        The integral of the correlation functions

    xs : [float]
        The time (x) coordinate

    sigma : [float] or numpy.ndarray()
        Optional standard error of y

    Returns
    -------
    a : [float]
        The list of coefficients
    tau : [float]
        The time constants
    a_err : [float]
        The standard error of the coefficients.
    tau_err : [float]
        The standard error of the time constants
    n : int
        The point in the x (time) vector that is 3*tau
    """
    # frequencies(y)
    dx = xs[1] - xs[0]
    width = int(0.4 / dx)
    ma = moving_average(y, width)

    # Find the initial, steep rise. Ignore first couple points.
    for i in range(2, len(ma) // 2):
        if ma[i] > 0.9 * ma[i + width]:
            break
    npts = i + width
    tau_guess = xs[npts]
    a_guess = ma[i]
    npts *= 20
    if 2 * npts > y.size:
        npts = y.size // 2

    sigma0 = sigma + np.average(sigma[0 : y.size // 10])

    # First fit with single exponential
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        warnings.simplefilter("ignore", OptimizeWarning)
        popt1, pcov1, infodict1, msg1, ierr1 = curve_fit(
            exp_1,
            xs[1:npts],
            y[1:npts],
            full_output=True,
            sigma=sigma0[1:npts],
            absolute_sigma=True,
            p0=[a_guess, tau_guess],
        )

    err = np.sqrt(np.diag(pcov1)).tolist()
    a1 = popt1[0]
    tau1 = popt1[1]

    # find a range for the next fit to keep it a "good" part of the data
    i_min = int(tau1 / dx / 2)  # tau / 2
    i_max = int(tau1 / dx * 20)  # 20 * tau1

    # And add a second exponential
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            warnings.simplefilter("ignore", OptimizeWarning)
            popt, pcov, infodict, msg, ierr = curve_fit(
                exp_2,
                xs[i_min:i_max],
                y[i_min:i_max],
                full_output=True,
                sigma=sigma[i_min:i_max],
                absolute_sigma=True,
                p0=[a1, tau1, 0.1 * a1, 5 * tau1],
            )
        err = np.sqrt(np.diag(pcov)).tolist()
        a = [popt[0], popt[2]]
        a_err = [err[0], err[2]]
        tau = [popt[1], popt[3]]
        tau_err = [err[1], err[3]]
        kappa = a[0] + a[1]
        kappa_err = a_err[0] + a_err[1]
        if kappa > 2 * a1:  # Shouldn't change this much!
            err = np.sqrt(np.diag(pcov1)).tolist()
            a = [popt1[0]]
            a_err = [err[0]]
            tau = [popt1[1]]
            tau_err = [err[1]]
            kappa = a[0]
            kappa_err = a_err[0]
    except Exception:
        err = np.sqrt(np.diag(pcov1)).tolist()
        a = [popt1[0]]
        a_err = [err[0]]
        tau = [popt1[1]]
        tau_err = [err[1]]
        kappa = a[0]
        kappa_err = a_err[0]

    # Find point for 3 * tau
    tmp = max(tau)
    nf = int(tmp / dx) * 20
    if nf > len(xs):
        nf = len(xs)

    fy = []
    if len(a) == 1:
        for t in xs[1:nf]:
            fy.append(exp_1(t, a[0], tau[0]))
    else:
        for t in xs[1:nf]:
            fy.append(exp_2(t, a[0], tau[0], a[1], tau[1]))

    return kappa, kappa_err, a, a_err, tau, tau_err, xs[1:nf], fy


def fit_helfand_moment(y, xs, sigma=None, start=1):
    """Find the best linear fit to longest possible segment.

    Parameters
    ----------
    y : [float] or numpy.ndarray()
        The Helfand moment

    xs : [float]
        The time (x) coordinate

    sigma : [float] or numpy.ndarray()
        Optional standard error of y

    Returns
    -------
    slope : float
        The fit slope.
    stderr : float
        The 95% standard error of the slope
    xs : [float]
        The x values (time) for the fit curve
    ys : [float]
        The y values for the fit curve.
    """
    dx = xs[1] - xs[0]

    # We know the curves curve near the origin, so ignore the first part
    i = int(start / dx)
    if i > len(y):
        i = int(1 / dx)

    popt, pcov, infodict, msg, ierr = curve_fit(
        axb,
        xs[i:],
        y[i:],
        full_output=True,
        sigma=sigma[i:],
        absolute_sigma=True,
    )
    slope = float(popt[0])
    b = float(popt[1])
    err = float(np.sqrt(np.diag(pcov)[0]))

    ys = []
    for x in xs[i:]:
        ys.append(axb(x, slope, b))

    return slope, err, xs[i:], ys


def get_helfand_slope(y, xs, sigma=None):
    """Get the instantaneous slope of the Helfand moments

    Parameters
    ----------
    y : [float] or numpy.ndarray()
        The Helfand moment

    xs : [float]
        The time (x) coordinate

    sigma : [float] or numpy.ndarray()
        Optional standard error of y

    Returns
    -------
    slope : float
        The fit slope.
    stderr : float
        The 95% standard error of the slope
    xs : [float]
        The x values (time) for the fit curve
    ys : [float]
        The y values for the fit curve.
    """
    n = len(y)

    dx = xs[1] - xs[0]
    slope = np.zeros_like(y)
    slope[2:] = (y[2:] - y[0 : n - 2]) / (2 * dx)
    new_x = xs
    new_sigma = np.zeros_like(y)
    new_sigma[2:] = sigma[2:] + sigma[0 : n - 2]
    new_sigma[0] = new_sigma[3]
    new_sigma[1] = new_sigma[3]
    return slope, new_x, new_sigma


def plot_correlation_functions(
    figure, CCF, ts, err=None, fit=None, labels=tensor_labels
):
    """Create a plot for the heat flux cross-correlation functions.

    Parameters
    ----------
    figure : seamm_util.Figure
        The figure that contains the plots.

    CCF : numpy.mdarray(6, m)
        The cross-correlation functions in W^2/m^4

    ts : [float]
        The times associated with the CCF, in ps

    err : numpy.ndarray(6, m)
        The standard errors of the CCF (optional)

    fit :
        The fit parameters for any fit of the CFF (optional)
    """
    plot = figure.add_plot("CCF")

    x_axis = plot.add_axis("x", label="t (ps)")
    y_axis = plot.add_axis("y", label="J0Jt (W^2/m^4)", anchor=x_axis)
    x_axis.anchor = y_axis

    for i, tmp in enumerate(labels):
        label, color, colora = tmp
        if fit is not None:
            hover = f"{label} = {fit[i]['kappa']:.3f} ± {fit[i]['stderr']:.3f} W/m/K"
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"fit{label}",
                hovertemplate=hover,
                x=fit[i]["xs"],
                xlabel="t",
                xunits="ps",
                y=fit[i]["ys"],
                ylabel=f"fit{label}",
                yunits="W/m/K*ps",
                color=color,
                dash="dash",
                width=3,
            )
        if err is not None:
            errs = np.concatenate((CCF[i] + err[i], CCF[i, ::-1] - err[i, ::-1]))
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"±{label}",
                x=ts + ts[::-1],
                xlabel="t",
                xunits="ps",
                y=errs.tolist(),
                ylabel=f"±{label}",
                yunits="W^2/m^4",
                color=colora,
                fill="toself",
                visible="legendonly",
            )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name=f"J0Jt{label}",
            x=ts,
            xlabel="t",
            xunits="ps",
            y=CCF[i].tolist(),
            ylabel=f"J0Jt{label}",
            yunits="W^2/m^4",
            color=color,
        )
    return plot


def plot_GK_integrals(
    figure, x, ts, err=None, fit=None, labels=tensor_labels, _range=None
):
    """Create a plot of the Green-Kubo integrals

    Parameters
    ----------
    figure : seamm_util.Figure
        The figure that contains the plots.

    x : numpy.mdarray(6, m)
        The cumulative integrals in W/m/K

    ts : [float]
        The times associated with the integrals, in ps

    err : numpy.ndarray(6, m)
        The standard errors of the integrals (optional)

    fit :
        The fit parameters for any fit of the CFF (optional)

    _range : [float]
        The range of the x-axis in terms of units in x.
    """
    plot = figure.add_plot("GKI")

    x_axis = plot.add_axis("x", label="t (ps)")
    y_axis = plot.add_axis("y", label="Kappa (W/m/K)", anchor=x_axis)
    x_axis.anchor = y_axis

    if _range is not None:
        x_axis["range"] = _range

    for i, tmp in enumerate(labels):
        label, color, colora = tmp
        if fit is not None and i < len(fit):
            hover = f"{label} = {fit[i]['kappa']:.3f} ± {fit[i]['stderr']:.3f} W/m/K"
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"fit{label}",
                hovertemplate=hover,
                x=fit[i]["xs"],
                xlabel="t",
                xunits="ps",
                y=fit[i]["ys"],
                ylabel=f"fit{label}",
                yunits="W/m/K",
                color=color,
                dash="dash",
                width=3,
            )

        if err is not None:
            errs = np.concatenate((x[i] + err[i], x[i, ::-1] - err[i, ::-1]))
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"±{label}",
                x=ts + ts[::-1],
                xlabel="t",
                xunits="ps",
                y=errs.tolist(),
                ylabel=f"±{label}",
                yunits="W/m/K",
                color=colora,
                fill="toself",
                visible="legendonly",
            )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name=f"K{label}",
            x=ts,
            xlabel="t",
            xunits="ps",
            y=x[i].tolist(),
            ylabel=f"K{label}",
            yunits="W/m/K",
            color=color,
        )
    return plot


def plot_helfand_moments(figure, M, ts, err=None, fit=None, labels=tensor_labels):
    """Create a plot for the Helfand moments.

    Parameters
    ----------
    figure : seamm_util.Figure
        The figure that contains the plots.

    M : numpy.mdarray(6, m)
        The Helfand moments, in W/m/K*ps

    ts : [float]
        The times associated with the moments, in ps
    """
    plot = figure.add_plot("M")

    x_axis = plot.add_axis("x", label="Time (ps)")
    y_axis = plot.add_axis("y", label="M (W/m/K*ps)", anchor=x_axis)
    x_axis.anchor = y_axis

    for i, tmp in enumerate(labels):
        label, color, colora = tmp
        if fit is not None:
            hover = f"κ{label} = {fit[i]['kappa']:.3f} ± {fit[i]['stderr']:.3f} W/m/K"
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"fit{label}",
                hovertemplate=hover,
                x=fit[i]["xs"],
                xlabel="t",
                xunits="ps",
                y=fit[i]["ys"],
                ylabel=f"fit{label}",
                yunits="W/m/K*ps",
                color=color,
                dash="dash",
                width=3,
            )
        if err is not None:
            errs = np.concatenate((M[i] + err[i], M[i, ::-1] - err[i, ::-1]))
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"±{label}",
                x=ts + ts[::-1],
                xlabel="t",
                xunits="ps",
                y=errs.tolist(),
                ylabel=f"±{label}",
                yunits="W/m/K*ps",
                color=colora,
                fill="toself",
                visible="legendonly",
            )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name=f"M{label}",
            x=ts,
            xlabel="t",
            xunits="ps",
            y=M[i, :].tolist(),
            ylabel=f"M{label}",
            yunits="W/m/K*ps",
            color=color,
        )
    return plot


def plot_helfand_slopes(
    figure, x, ts, err=None, fit=None, labels=tensor_labels, _range=None
):
    """Create a plot of the slope of the Helfand moments

    Parameters
    ----------
    figure : seamm_util.Figure
        The figure that contains the plots.

    x : numpy.mdarray(6, m)
        The slope in W/m/K

    ts : [float]
        The times associated with the points, in ps

    err : numpy.ndarray(6, m)
        The standard errors of the points (optional)

    fit :
        The fit parameters for any fit of the slope (optional)

    _range : [float]
        The range of the x-axis in terms of units in x.
    """
    plot = figure.add_plot("Slope")

    x_axis = plot.add_axis("x", label="t (ps)")
    y_axis = plot.add_axis("y", label="Kappa (W/m/K)", anchor=x_axis)
    x_axis.anchor = y_axis

    if _range is not None:
        x_axis["range"] = _range

    for i, tmp in enumerate(labels):
        label, color, colora = tmp
        if fit is not None and i < len(fit):
            hover = f"{label} = {fit[i]['kappa']:.3f} ± {fit[i]['stderr']:.3f} W/m/K"
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"fit{label}",
                hovertemplate=hover,
                x=fit[i]["xs"],
                xlabel="t",
                xunits="ps",
                y=fit[i]["ys"],
                ylabel=f"fit{label}",
                yunits="W/m/K",
                color=color,
                dash="dash",
                width=3,
            )
        if err is not None:
            errs = np.concatenate((x[i] + err[i], x[i][::-1] - err[i][::-1]))
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"±{label}",
                x=ts + ts[::-1],
                xlabel="t",
                xunits="ps",
                y=errs.tolist(),
                ylabel=f"±{label}",
                yunits="W/m/K",
                color=colora,
                fill="toself",
                visible="legendonly",
            )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name=f"K{label}",
            x=ts,
            xlabel="t",
            xunits="ps",
            y=x[i].tolist(),
            ylabel=f"K{label}",
            yunits="W/m/K",
            color=color,
        )
    return plot
