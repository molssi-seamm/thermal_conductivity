# -*- coding: utf-8 -*-

"""The graphical part of a Thermal Conductivity step"""

import pprint  # noqa: F401
import tkinter as tk
import tkinter.ttk as ttk

from .thermal_conductivity_parameters import ThermalConductivityParameters
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw


class TkThermalConductivity(seamm.TkNode):
    """
    The graphical part of a Thermal Conductivity step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in Thermal Conductivity_parameters.py

    See Also
    --------
    ThermalConductivity, TkThermalConductivity,
    ThermalConductivityParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        namespace="org.molssi.seamm.tk",
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.namespace = namespace
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )
        self.create_dialog()

    def create_dialog(self):
        """
        Create the dialog. A set of widgets will be chosen by default
        based on what is specified in the Thermal Conductivity_parameters
        module.

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkThermalConductivity.reset_dialog
        """

        frame = super().create_dialog(title="Thermal Conductivity", widget="notebook")
        # make it large!
        screen_w = self.dialog.winfo_screenwidth()
        screen_h = self.dialog.winfo_screenheight()
        w = int(0.9 * screen_w)
        h = int(0.8 * screen_h)
        x = int(0.05 * screen_w / 2)
        y = int(0.1 * screen_h / 2)

        self.dialog.geometry(f"{w}x{h}+{x}+{y}")

        # Add a frame for the flowchart
        notebook = self["notebook"]
        flowchart_frame = ttk.Frame(notebook)
        self["flowchart frame"] = frame
        notebook.add(flowchart_frame, text="Flowchart", sticky=tk.NSEW)

        self.tk_subflowchart = seamm.TkFlowchart(
            master=flowchart_frame,
            flowchart=self.node.subflowchart,
            namespace=self.namespace,
        )
        self.tk_subflowchart.draw()

        # Fill in the control parameters
        # Shortcut for parameters
        P = self.node.parameters

        # conductivity frame to isolate widgets
        frame = self["conductivity_frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Thermal Conductivity",
            labelanchor="n",
            padding=10,
        )

        for key in ThermalConductivityParameters.parameters:
            self[key] = P[key].widget(frame)

        # and binding to change as needed
        self["approach"].combobox.bind(
            "<<ComboboxSelected>>", self.reset_conductivity_frame
        )

        # and lay them out
        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Layout the widgets in the dialog.

        The widgets are chosen by default from the information in
        Thermal Conductivity_parameter.

        This function simply lays them out row by row with
        aligned labels. You may wish a more complicated layout that
        is controlled by values of some of the control parameters.
        If so, edit or override this method

        Parameters
        ----------
        widget : Tk Widget = None

        Returns
        -------
        None

        See Also
        --------
        TkThermalConductivity.create_dialog
        """

        row = 0

        self["conductivity_frame"].grid(row=row, column=0, sticky=tk.EW, pady=10)
        row += 1
        self.reset_conductivity_frame()

        return row

    def reset_conductivity_frame(self, widget=None):
        """Layout the widgets in the conductivity frame
        as needed for the current state"""

        approach = self["approach"].get()

        frame = self["conductivity_frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0

        # Main controls
        widgets = []
        for key in ("approach", "T", "nruns", "errors"):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.W)
            widgets.append(self[key])
            row += 1

        widgets0 = []
        if "RNEMD" in approach:
            for key in ("nlayers", "swap frequency", "n to swap"):
                self[key].grid(row=row, column=1, sticky=tk.W)
                widgets0.append(self[key])
                row += 1
        elif "Thermostatted" in approach:
            for key in ("nlayers", "deltaT"):
                self[key].grid(row=row, column=1, sticky=tk.W)
                widgets0.append(self[key])
                row += 1
        elif "Heat Removal" in approach:
            for key in ("nlayers", "ehex"):
                self[key].grid(row=row, column=1, sticky=tk.W)
                widgets0.append(self[key])
                row += 1

        sw.align_labels(widgets, sticky=tk.E)
        if len(widgets0) > 0:
            sw.align_labels(widgets0, sticky=tk.E)

        frame.columnconfigure(0, minsize=50)

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkThermalConductivity.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
