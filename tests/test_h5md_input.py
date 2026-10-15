#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests that Thermal Conductivity reads md.h5md."""

import numpy as np
import pytest

from seamm_md_analysis import H5MDWriter

from thermal_conductivity_step.thermal_conductivity import ThermalConductivity

cell = np.diag([10.0, 20.0, 30.0])


def write(path, periodic=True, volume=False):
    steps = np.arange(0, 100, 2)
    times = steps * 1.0
    with H5MDWriter(path) as w:
        w.particles(ids=[1, 2], cell=cell if periodic else None, periodic=periodic)
        flux = np.column_stack([steps, -steps, 2 * steps]).astype(float)
        w.observable("heat_flux", steps, times, flux, "W m-2")
        w.observable("temperature", steps, times, np.full(len(steps), 300.0), "K")
        if volume:
            w.observable(
                "volume", steps, times, np.full(len(steps), 5000.0), "Angstrom3"
            )


def test_reads_flux_temperature_and_cell_volume(tmp_path):
    write(tmp_path / "md.h5md")
    tc = ThermalConductivity.__new__(ThermalConductivity)
    J, dt, T, V = tc._read_h5md(1, tmp_path)
    assert J.shape == (3, 50)
    assert np.allclose(J[1], -J[0])
    assert dt.m_as("fs") == pytest.approx(2.0)
    assert T == pytest.approx(300.0)
    assert V == pytest.approx(6000.0)


def test_npt_volume_is_used(tmp_path):
    write(tmp_path / "md.h5md", volume=True)
    tc = ThermalConductivity.__new__(ThermalConductivity)
    assert tc._read_h5md(1, tmp_path)[3] == pytest.approx(5000.0)


def test_non_periodic_is_refused(tmp_path):
    write(tmp_path / "md.h5md", periodic=False)
    tc = ThermalConductivity.__new__(ThermalConductivity)
    with pytest.raises(RuntimeError, match="not periodic"):
        tc._read_h5md(1, tmp_path)


def test_no_h5md_returns_none(tmp_path):
    tc = ThermalConductivity.__new__(ThermalConductivity)
    assert tc._read_h5md(1, tmp_path) is None
