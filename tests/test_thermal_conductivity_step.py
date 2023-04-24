#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `thermal_conductivity_step` package."""

import pytest  # noqa: F401
import thermal_conductivity_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = thermal_conductivity_step.ThermalConductivity()
    assert str(type(result)) == (
        "<class 'thermal_conductivity_step.thermal_conductivity"
        ".ThermalConductivity'>"
    )
