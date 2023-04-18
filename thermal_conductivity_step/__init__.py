# -*- coding: utf-8 -*-

"""
thermal_conductivity_step
A SEAMM plug-in for Thermal Conductivity
"""

# Bring up the classes so that they appear to be directly in
# the thermal_conductivity_step package.

from .thermal_conductivity import ThermalConductivity  # noqa: F401, E501
from .thermal_conductivity_step import ThermalConductivityStep  # noqa: F401, E501
from .tk_thermal_conductivity import TkThermalConductivity  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
