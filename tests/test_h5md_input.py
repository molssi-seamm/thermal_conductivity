#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests that Thermal Conductivity reads md.h5md."""

import numpy as np
import pytest

from seamm_md_analysis import H5MDWriter

from thermal_conductivity_step.thermal_conductivity import ThermalConductivity
from seamm_util import Q_

cell = np.diag([10.0, 20.0, 30.0])


def write(path, periodic=True, volume=False, species=None):
    steps = np.arange(0, 100, 2)
    times = steps * 1.0
    with H5MDWriter(path) as w:
        w.particles(ids=[1, 2], cell=cell if periodic else None, periodic=periodic)
        flux = np.column_stack([steps, -steps, 2 * steps]).astype(float)
        w.observable("heat_flux", steps, times, flux, "W m-2")
        w.observable("temperature", steps, times, np.full(len(steps), 300.0), "K")
        if species is not None:
            k = len(species)
            values = np.arange(len(steps) * k * 3, dtype=float).reshape(-1, k, 3)
            w.observable(
                "species_mass_flux",
                steps,
                times,
                values,
                "g mol-1 Angstrom fs-1",
                attrs={"species": species},
            )
        if volume:
            w.observable(
                "volume", steps, times, np.full(len(steps), 5000.0), "Angstrom3"
            )


def test_reads_flux_temperature_and_cell_volume(tmp_path):
    write(tmp_path / "md.h5md")
    J, dt, T, V, S, species = ThermalConductivity._read_h5md(tmp_path / "md.h5md")
    assert S is None and species is None
    assert J.shape == (3, 50)
    assert np.allclose(J[1], -J[0])
    assert dt.m_as("fs") == pytest.approx(2.0)
    assert T == pytest.approx(300.0)
    assert V == pytest.approx(6000.0)


def test_reads_species_fluxes(tmp_path):
    write(tmp_path / "md.h5md", species=["O", "[Na+]", "[Cl-]"])
    J, _, _, _, S, species = ThermalConductivity._read_h5md(tmp_path / "md.h5md")
    assert species == ["O", "[Na+]", "[Cl-]"]
    # All but the first species, as k x 3 x n
    assert S.shape == (2, 3, 50)
    assert S[0, 0, 0] == 3.0 and S[1, 2, 1] == 9 + 8.0
    assert J.shape == (3, 50)


def test_one_species_is_not_a_mixture(tmp_path):
    write(tmp_path / "md.h5md", species=["O"])
    result = ThermalConductivity._read_h5md(tmp_path / "md.h5md")
    assert result[4] is None and result[5] is None


def test_npt_volume_is_used(tmp_path):
    write(tmp_path / "md.h5md", volume=True)
    assert ThermalConductivity._read_h5md(tmp_path / "md.h5md")[3] == pytest.approx(
        5000.0
    )


def test_non_periodic_is_refused(tmp_path):
    write(tmp_path / "md.h5md", periodic=False)
    with pytest.raises(RuntimeError, match="not periodic"):
        ThermalConductivity._read_h5md(tmp_path / "md.h5md")


def test_no_heat_flux_returns_none(tmp_path):
    with H5MDWriter(tmp_path / "md.h5md") as w:
        w.particles(ids=[1, 2], cell=cell, periodic=True)
        steps = np.arange(0, 100, 2)
        w.observable("temperature", steps, steps * 1.0, np.full(50, 300.0), "K")
    assert ThermalConductivity._read_h5md(tmp_path / "md.h5md") is None


def test_blocks_split_one_run():
    J = np.ones((3, 40000))
    dt = Q_(2.0, "fs")
    blocks = ThermalConductivity._split_into_blocks((J, dt, 300.0, 5000.0, None), 4)
    assert len(blocks) == 4
    assert all(b[0].shape == (3, 10000) for b in blocks)
    assert blocks[2][2:] == (300.0, 5000.0, None)
    with pytest.raises(RuntimeError, match="too few"):
        ThermalConductivity._split_into_blocks(
            (J, Q_(10.0, "fs"), 300.0, 5000.0, None), 25
        )
    S = np.arange(2 * 3 * 40000.0).reshape(2, 3, 40000)
    blocks = ThermalConductivity._split_into_blocks((J, dt, 300.0, 5000.0, S), 4)
    assert blocks[1][4].shape == (2, 3, 10000)
    assert np.array_equal(blocks[1][4], S[..., 10000:20000])


def test_blocks_too_short():
    """Each block must be long enough for the Green-Kubo fit, 10 ps."""
    J = np.ones((3, 10000))
    with pytest.raises(RuntimeError, match="too short"):
        ThermalConductivity._split_into_blocks(
            (J, Q_(2.0, "fs"), 300.0, 5000.0, None), 4
        )


def test_requirements():
    tc = ThermalConductivity.__new__(ThermalConductivity)
    requirements = tc.md_requirements()
    (need,) = requirements.needs
    assert need.names == ["heat_flux"]
    assert need.max_steps == 2
    assert requirements.periodic
    assert requirements.dynamics


def test_external_engine_is_recognized(tmp_path):
    """A run whose forces came from an MDI engine (an MLFF) has no per-atom energies
    or stresses, so its heat flux is incomplete (ChemAI job 5357, 2026-10-21)."""
    from thermal_conductivity_step.thermal_conductivity import external_engine

    write(tmp_path / "md.h5md")
    assert external_engine(tmp_path / "md.h5md") is None  # no fixes recorded
    with H5MDWriter(tmp_path / "mdi.h5md") as w:
        w.particles(ids=[1, 2], cell=cell, periodic=True)
        w.parameters({"fixes": ["mdi/qm", "nvt"]})
    assert external_engine(tmp_path / "mdi.h5md") == "mdi/qm"
    with H5MDWriter(tmp_path / "ff.h5md") as w:
        w.particles(ids=[1, 2], cell=cell, periodic=True)
        w.parameters({"fixes": ["nvt", "shake"]})
    assert external_engine(tmp_path / "ff.h5md") is None
