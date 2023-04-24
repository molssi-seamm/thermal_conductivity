# -*- coding: utf-8 -*-
"""
Control parameters for the Thermal Conductivity step in a SEAMM flowchart
"""
import logging

import seamm

logger = logging.getLogger(__name__)


class ThermalConductivityParameters(seamm.Parameters):
    """
    The control parameters for Thermal Conductivity.

    You need to replace the "time" entry in dictionary below these comments with the
    definitions of parameters to control this step. The keys are parameters for the
    current plugin,the values are dictionaries as outlined below.

    Examples
    --------
    ::

        parameters = {
            "time": {
                "default": 100.0,
                "kind": "float",
                "default_units": "ps",
                "enumeration": tuple(),
                "format_string": ".1f",
                "description": "Simulation time:",
                "help_text": ("The time to simulate in the dynamics run.")
            },
        }

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.
        Each key of the dictionary is a dictionary that contains the
        the following keys:

    parameters["default"] :
        The default value of the parameter, used to reset it.

    parameters["kind"] : enum()
        Specifies the kind of a variable. One of  "integer", "float", "string",
        "boolean", or "enum"

        While the "kind" of a variable might be a numeric value, it may still have
        enumerated custom values meaningful to the user. For instance, if the parameter
        is a convergence criterion for an optimizer, custom values like "normal",
        "precise", etc, might be adequate. In addition, any parameter can be set to a
        variable of expression, indicated by having "$" as the first character in the
        field. For example, $OPTIMIZER_CONV.

    parameters["default_units"] : str
        The default units, used for resetting the value.

    parameters["enumeration"]: tuple
        A tuple of enumerated values.

    parameters["format_string"]: str
        A format string for "pretty" output.

    parameters["description"]: str
        A short string used as a prompt in the GUI.

    parameters["help_text"]: str
        A longer string to display as help for the user.

    See Also
    --------
    ThermalConductivity, TkThermalConductivity, ThermalConductivityParameters,
    ThermalConductivityStep
    """

    parameters = {
        "approach": {
            "default": "Green-Kubo",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "Green-Kubo",
                "Reverse Non-Equilibrium MD (RNEMD)",
                "Thermostatted Temperature Difference",
                "Heat Removal / Addition",
            ),
            "format_string": "",
            "description": "Approach:",
            "help_text": (
                "The approach or method for determining the thermal conductivity."
            ),
        },
        "T": {
            "default": 298.15,
            "kind": "float",
            "default_units": "K",
            "format_string": ".2f",
            "description": "Temperature:",
            "help_text": "The temperature for the thermal conductivity calculation.",
        },
        "nruns": {
            "default": "20",
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of runs to average:",
            "help_text": "The number for separate runs to average.",
        },
        "errors": {
            "default": "continue to next run",
            "kind": "string",
            "default_units": "",
            "enumeration": (
                "continue to next run",
                "exit the thermal conductivity step",
                "stop the job",
            ),
            "format_string": "s",
            "description": "On errors",
            "help_text": "How to handle errors in the runs",
        },
        "ehex": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Use enhanced heat exchange algorithm:",
            "help_text": (
                "Whether to use the enhanced algorithm, which conserves energy "
                " noticeably better than the original algorithm at little extra cost."
            ),
        },
        "nlayers": {
            "default": "50",
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of layers:",
            "help_text": (
                "How many layers to slice the system into for determining the "
                "temperature profile."
            ),
        },
        "deltaT": {
            "default": "20",
            "kind": "float",
            "default_units": "K",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Target temperature difference:",
            "help_text": (
                "The target difference in temperatures between the hot and cold "
                "zones."
            ),
        },
        "swap frequency": {
            "default": "100",
            "kind": "float",
            "default_units": "fs",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Momenta swap frequency:",
            "help_text": "The frequency to swap momenta in the RNEMD approach.",
        },
        "n to swap": {
            "default": "1",
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of momenta to swap:",
            "help_text": "The number of momenta to swap in the RNEMD approach.",
        },
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
