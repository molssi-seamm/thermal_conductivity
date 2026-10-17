#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Regression tests of the analysis on a heat flux with a known thermal conductivity.

Each component of the flux is an Ornstein-Uhlenbeck process, whose autocorrelation
function is σ² exp(-t/τ), so its Green-Kubo integral is σ² τ and
κ = V σ² τ / (k_B T²).
"""

import numpy as np
import pytest

from seamm_md_analysis import isotropic_average
from seamm_util import Q_
from thermal_conductivity_step.thermal_conductivity import ThermalConductivity

T = 300.0  # K
V = 8000.0  # Å^3
dt = Q_(2.0, "fs")
tau = Q_(0.1, "ps")
sigma = 1.0e10  # W/m^2


def expected_kappa():
    k_B = Q_("k_B")
    value = Q_(V, "Å^3") * Q_(sigma, "W/m^2") ** 2 * tau / (k_B * Q_(T, "K") ** 2)
    return value.m_as("W/m/K")


def ou_flux(n, rng):
    """Three independent Ornstein-Uhlenbeck series of n samples."""
    a = np.exp(-(dt / tau).m_as(""))
    noise = rng.normal(size=(3, n)) * sigma * np.sqrt(1 - a**2)
    J = np.empty((3, n))
    J[:, 0] = rng.normal(size=3) * sigma
    for i in range(1, n):
        J[:, i] = a * J[:, i - 1] + noise[:, i]
    return J


@pytest.fixture(scope="module")
def analyzed():
    """Five runs of 20 ps, processed as the step does."""
    rng = np.random.default_rng(12345)
    tc = ThermalConductivity.__new__(ThermalConductivity)
    tc._reset()
    for run in range(1, 6):
        tc.process_run(run, (ou_flux(10000, rng), dt, T, V))
    tc._combine()
    ts = (np.arange(tc.M.shape[1]) * dt.m_as("ps")).tolist()
    gk_ts = (np.arange(tc.GK_integral.shape[1]) * dt.m_as("ps")).tolist()
    per_sample, off_diagonal = tc._per_sample(ts, gk_ts, [1.0] * 6)
    return tc, per_sample, off_diagonal


def test_kappa_by_each_method(analyzed):
    _, per_sample, _ = analyzed
    kappa = expected_kappa()
    for method, values in per_sample.items():
        assert values.shape == (5, 3)
        assert np.all(np.isfinite(values)), method
        # The mean is within 15%, and within its 95% interval
        estimate = isotropic_average(values)
        assert estimate.n == 5
        assert estimate.value == pytest.approx(kappa, rel=0.15), method
        assert abs(estimate.value - kappa) < estimate.ci95, method


def test_off_diagonal_terms_are_zero(analyzed):
    _, _, off_diagonal = analyzed
    assert off_diagonal.shape == (5, 3)
    text = ThermalConductivity._off_diagonal_check(off_diagonal)
    assert text.startswith("The off-diagonal terms are zero within their errors")


def test_off_diagonal_check_flags_nonzero_terms():
    rng = np.random.default_rng(1)
    values = 1.0 + 0.01 * rng.normal(size=(5, 3))
    text = ThermalConductivity._off_diagonal_check(values)
    assert "not zero" in text
    assert ThermalConductivity._off_diagonal_check(values[:1]) == ""


def test_curves_are_averaged(analyzed):
    tc, _, _ = analyzed
    assert tc.M.shape == tc.Ms[0].shape
    assert np.allclose(tc.M, np.mean(tc.Ms, axis=0))
    assert np.allclose(tc.GK_integral_err, np.std(tc.GK_integrals, axis=0))


def test_short_block_is_refused():
    tc = ThermalConductivity.__new__(ThermalConductivity)
    tc._reset()
    with pytest.raises(RuntimeError, match="Block 1 is only"):
        tc.process_run(1, (np.ones((3, 1000)), dt, T, V), what="Block")
