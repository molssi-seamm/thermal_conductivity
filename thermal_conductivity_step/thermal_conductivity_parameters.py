# -*- coding: utf-8 -*-
"""
Control parameters for the Thermal Conductivity step in a SEAMM flowchart
"""

from seamm_md_analysis.step_support import source_parameters
import logging

import seamm

logger = logging.getLogger(__name__)


class ThermalConductivityParameters(seamm.Parameters):
    """
    The control parameters for Thermal Conductivity.

    The MD to analyze is run earlier in the flowchart -- in a Loop for independent
    runs -- or in other jobs (see seamm_md_analysis.step_support.source_parameters).
    The thermal conductivity is calculated from the heat flux by the Green-Kubo and
    Helfand methods.
    """

    # Dropped with the own-MD route (phase 8): ignored in older flowcharts
    retired = (
        "nruns",
        "errors",
        "approach",
        "T",
        "ehex",
        "nlayers",
        "deltaT",
        "swap frequency",
        "n to swap",
    )

    parameters = {
        **source_parameters("the thermal conductivity"),
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("ThermalConductivityParameters.__init__")

        super().__init__(
            defaults={**ThermalConductivityParameters.parameters, **defaults}, data=data
        )
