# -*- coding: utf-8 -*-

"""
thermal_conductivity_step
A SEAMM plug-in for Thermal Conductivity
"""

# Bring up the classes so that they appear to be directly in
# the thermal_conductivity_step package.

from .thermal_conductivity import ThermalConductivity
from .thermal_conductivity_parameters import ThermalConductivityParameters
from .thermal_conductivity_step import ThermalConductivityStep
from .tk_thermal_conductivity import TkThermalConductivity

from .metadata import metadata

from .analysis import (
    create_correlation_functions,
    create_helfand_moments,
    fit_green_kubo_integral,
    fit_helfand_moment,
    plot_correlation_functions,
    plot_GK_integrals,
    plot_helfand_moments,
)

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
