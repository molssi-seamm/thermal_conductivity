=======
History
=======
2026.10.18 -- Bugfix: the thermal conductivity of mixtures was too large
    * In a mixture the heat flux includes the energy carried by the interdiffusion of
      the species, which is not conduction of heat, so κ was too large: for 1 m NaCl in
      SPC/E water about 1.15 W/m/K instead of 0.81. The step now removes it, using the
      mass fluxes of the species that the LAMMPS step samples with the heat flux
      (Bertossa, Grasselli, Ercole and Baroni, Phys. Rev. Lett. 122, 255901 (2019)), by
      both Green-Kubo and the Helfand moments. The correction for one run or block is
      too noisy, so it uses the averages over the runs or blocks, with errors by the
      jackknife. The notes give the uncorrected value too.
    * For a mixture whose MD did not sample the mass fluxes, made with an older LAMMPS
      step, the step warns that κ is too large.

2026.10.17 -- Analyze MD run before the step; one isotropic κ with labelled errors
    * The step analyzes MD run before it in the flowchart -- the runs of a Loop around
      the MD, as independent replicas ("the last MD run of each loop iteration", the
      default) -- or md.h5md files; it no longer runs its own MD in a sub-flowchart.
      Older flowcharts that stored its settings for that can still be read; the check
      before the job runs says what replaces them.
    * The check before the job runs also warns when the analysis is inside the loop
      holding its MD, when the loop changes the temperature or pressure, and when the
      thermostat does not suit the analysis, with fixes; the step warns before the MD
      when the system is small.
    * Errors are the 95% confidence interval, Student's t for the number of runs or
      blocks, and the output says so; stored errors are standard errors.
    * The state of the runs -- temperature, pressure, volume and density with their
      errors, and the equilibrium density where an NPT run before them gives the bulk
      modulus -- is reported and stored. The values of each run or block are in
      step.out.
    * κ is the average of the three diagonal components within each run, with a check
      that they agree and that the off-diagonal terms are zero; κ and its components
      are stored as results.
    * Bugfix: in a loop, the step kept the data of earlier iterations.

2026.10.16 -- Bugfix: a short run stopped the analysis with an error
    * A run, or block of a run, shorter than about 10 ps stopped the analysis with an
      IndexError. Each must now be at least 10 ps long, and the step says so.

2026.10.15.2 -- Analyze MD run earlier in the flowchart
    * New "MD to analyze" choice: run the MD in the sub-flowchart (as before), the last
      MD run earlier in the flowchart, the last MD run of each loop iteration, or
      md.h5md files, including other jobs'. A single run is divided into blocks to
      estimate the errors.
    * Before the job runs, the step checks that the MD will sample the heat flux every
      two timesteps, and the job does not start if it does not.
    * The dialog shows only the settings that apply, and hides the sub-flowchart when
      it is not used.
    * When the run (or block) is too short to fit the Helfand moments, the output now
      says so; before, those rows were silently missing.

2026.10.15.1 -- Reads the H5MD files; needs a periodic system
    * Reads the heat flux, temperature and volume from the ``md.h5md`` files that the
      LAMMPS step now writes, and from the older text files for runs made before.
    * Stops before running any MD if the system is not periodic, since the heat flux is
      divided by the volume of the cell.
    * The temperature is now averaged over the same samples as the heat flux, the whole
      run. It was averaged only after the equilibration detected in the run, while the
      heat flux was used whole. The thermal conductivity changes very little, e.g. by
      0.05% for liquid argon.
    * Bugfix: a run whose temperature and volume were missing stopped with an obscure
      error; the message now says what is missing. The sampling interval was taken to be
      in fs even when the file gave it in ps.

2026.10.15 -- Bugfix: crash with coarse sampling; unreliable fits
    * The Green–Kubo fit failed when the heat flux was sampled less often than every
      0.4 ps. It no longer fails.
    * The Green–Kubo fit could accept a second exponential with a negative amplitude,
      or one too slow to level off within the data, and give a value far from the
      plateau of the integral. Such fits are now rejected in favour of a single
      exponential. The error of a two-exponential fit, which could be enormous when
      the two time constants were nearly equal, now accounts for their correlation.
    * The errors were from the fits to the curves averaged over the runs, which
      treat their strongly correlated points as independent: the Green-Kubo and
      Helfand-moment errors were far too small (e.g. +/- 0.0004 W/m/K where the runs
      differed by 0.03), and the Helfand-derivative errors too large. The error of
      each value is now the standard error of the values fitted to each run, and the
      95% interval is twice that.
    * With a single run the Helfand-moment fits failed and printed 1.00; they now give
      the slope. With a single run there is no spread between runs to estimate
      errors from, so none are printed.
    * The slope of the Helfand moments ("Helfand Derivative") was shifted by one
      sample, and its errors were in the wrong units. Both are corrected.
    * The analysis now uses the shared ``seamm_md_analysis`` library, also used by the
      Diffusivity step, instead of its own copy. The Helfand moments are computed much
      faster (by FFT): for a million samples, 0.4 s instead of about 2 minutes.
    * Thermal conductivities from the LAMMPS Heat Flux sub-step need lammps-step
      2026.10.15 or later: earlier versions computed the heat flux wrongly.
    * The step no longer needs pandas or statsmodels, and its tests run with uv in CI
      using the package's declared dependencies.

2026.3.1 -- Internal: switching from deprecated library pkg_resources to importlib

2024.6.3 -- Bugfix: handling of options for subflowchart
    * Fixed a bug where the options for the subflowchart were not being parsed
      correctly.
      
2024.3.22 -- Updated for new scheme for running background tasks.

2024.1.5 -- Bugfix: Thermal conductivity
    * If the Helfand moments fit in the thermal conductivity step failed it stopped the
      entire job. This is fixed, as well as some of the underlying causes for
      convergence issues.
      
2023.5.29 -- Converged with general approach for trajectory analysis

2023.5.6 -- Bugfix
    * Fixed an error handling Nan's and Inf's that caused a crash
    * Added the predictions from the derivatives of the Helfand moments to the output to
      give a better feel for the quality of the results.
      
2023.5.5 -- Improved analysis
    * Considerable improvements to the analysis, results now seem solid
    * Fixed issues with fitting the linear portion of the Helfand moments
    * Added plot of the slope from the Helfand moments, which is similar to the
      Green-Kubo integral.
    * Cleaned up both the output and graphs.
      
2023.4.24 -- Initial working version
    * Initial tests seem to work but needs more thorough testing.
    * Needs documentation!
      
2023.4.18 (2023-04-18)
    * Plug-in created using the SEAMM plug-in cookiecutter.
