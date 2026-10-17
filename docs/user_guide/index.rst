.. _user-guide:

**********
User Guide
**********

The Thermal Conductivity step calculates the thermal conductivity κ from the
fluctuations of the heat flux J in equilibrium molecular dynamics. By Green and Kubo,

.. math::

   \kappa_{ab} = \frac{V}{k_B T^2} \int_0^\infty \langle J_a(t) J_b(0) \rangle \, dt

with V the volume and T the temperature. The step reports κ in three ways, which
should agree within their errors:

Green-Kubo
    The plateau of the running integral of the autocorrelation function of the heat
    flux, fitted by one or, if justified by the data, two rising exponentials.
Helfand moments
    The slope at long times of the Helfand moments, the mean square of the integral
    of the heat flux.
Helfand derivative
    The plateau of the derivative of the Helfand moments, fitted like the Green-Kubo
    integral.

The heat flux must be sampled every two timesteps to follow its fastest oscillations;
sampling it more coarsely gives wrong results for molecules with fast vibrations, such
as water. The step needs a periodic system, since the heat flux is divided by the
volume of the cell.

Where the MD comes from
=======================

The step analyzes MD that ran before it in the flowchart, or the ``md.h5md`` files of
earlier runs. For independent runs, put the MD in a Loop -- in parallel, with each
run starting from the same structure ("Discard the structure" in the MD sub-steps)
and new velocities -- and analyze "the last MD run of each loop iteration", the
default. A typical flowchart:

.. code-block:: text

    Model Chemistry (FF@...), Packmol, minimization
    NPT equilibration; Set Cell to its mean volume; NVT equilibration
    Loop over the replicas (parallel)
        LAMMPS: velocities, NVT decorrelation, NVT production sampling the heat flux
    Thermal Conductivity

**MD to analyze** chooses the runs:

the last MD run of each loop iteration
    The run of each iteration of the loop around the MD, analyzed as independent
    runs; without a loop, the one run. This is the default.
the last MD run
    The last MD run earlier in the flowchart that sampled the heat flux. One
    production run can feed several analysis steps.
files
    The ``md.h5md`` files given, each an independent run; ``job://<n>/<path>`` reads a
    file of job ``n``.

The MD must sample the heat flux every two timesteps: "Heat flux" in the LAMMPS NVT
sub-step (the Heat Flux sub-step does too, but runs NVE). Checking the flowchart before it runs
reports what is missing, with a fix that turns the sampling on. It also warns when the
analysis is inside the loop holding its MD, when the loop changes the temperature or
pressure (the runs are then not replicas of one state), and when the thermostat
perturbs the dynamics (Langevin, Andersen) or there is none (NVE): use NVT with a
Nosé-Hoover or CSVR thermostat for the production runs. Before the MD runs, the step
warns if the system is small for the analysis (fewer than about 200 molecules).

Each run, and each block of a single run, must be at least 10 ps long, and much longer
than the correlation time of the heat flux; a few hundred picoseconds per run is
typical for a liquid.

Convergence and errors
======================

Each method is fitted to each run, and the three diagonal components are averaged
within each run. κ is the mean of these averages over the runs, and its error is from
their spread: the table gives the 95% confidence interval (Student's t times the
standard error) and says so. A single run is divided into **Blocks for errors from one
run** (5 by default), analyzed in the same way. With one run in one block there is no
error estimate.

For an isotropic run -- a cubic cell with nothing directional in it, such as walls or
an electric field -- the table gives one value per method, and a line checks that the
diagonal components agree within their errors (χ²) and that the off-diagonal terms
κ_xy, κ_xz and κ_yz, from the Helfand moments, are zero within theirs. Otherwise the
diagonal components are given too. The values of each run or block, with the
components, are in the step's detailed output (``step.out``).

The step graphs the correlation functions of the heat flux (``JJ.graph``), the
Green-Kubo integrals (``GK_integral.graph``), the Helfand moments
(``HelfandMoments.graph``) and their derivatives (``HelfandDerivatives.graph``),
averaged over the runs, with their fits.

Results
=======

* κ with its standard error by each method (``kappa``, ``kappa (Helfand)`` and
  ``kappa (Helfand derivative)``), and the diagonal components by Green-Kubo
  (``kappa_xx``, ``kappa_yy``, ``kappa_zz``);
* the state of the runs: temperature, pressure, volume and density with their errors,
  and for runs at constant volume the equilibrium density at the target pressure when
  an NPT run before them gives the bulk modulus;
* the number of runs.

κ by Green-Kubo and from the Helfand moments, with their errors, can be stored as
properties; all the results can be put in variables and tables.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
