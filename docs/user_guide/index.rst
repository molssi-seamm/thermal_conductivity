.. _user-guide:

**********
User Guide
**********

The Thermal Conductivity step calculates the thermal conductivity from the
fluctuations of the heat flux in equilibrium molecular dynamics. Its subflowchart
contains the molecular dynamics, normally a LAMMPS step ending with the Heat Flux
sub-step, and is run several times, each run starting from new velocities. The heat
flux is sampled every two timesteps, which is needed to follow its fastest
oscillations; sampling it more coarsely gives wrong results for molecules with fast
vibrations, such as water.

The step reports the thermal conductivity in three ways, which should agree within
their errors:

Green-Kubo
    The plateau of the running integral of the autocorrelation function of the heat
    flux, fitted by one or, if justified by the data, two rising exponentials.
Helfand Moments
    The slope of the Helfand moments, the integral of the heat flux squared, at
    long times. This also gives the off-diagonal elements of the tensor, which
    should be zero within their errors for a fluid.
Helfand Derivative
    The plateau of the derivative of the Helfand moments, fitted like the
    Green-Kubo integral.

Each value is fitted to the average over the runs. Its error, given as a 95%
interval (twice the standard error), comes from the spread of the values fitted to
each run separately, so at least two runs are needed for errors, and more runs give
better ones. With a single run no errors are reported.

The step needs a periodic system, since the heat flux is divided by the volume of the
cell; it stops before running any MD if the structure is not periodic. It reads the
heat flux, the temperature and the volume from the ``md.h5md`` file that the LAMMPS
step writes for each MD run. Runs made with older versions of the LAMMPS step, which
have no ``md.h5md``, are read from their text files as before.

..
    <remove the dots above and this line and unindent the toctree to expose it>
    Contents:

    .. toctree::
       :glob:
       :maxdepth: 2
       :titlesonly:

       *

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
