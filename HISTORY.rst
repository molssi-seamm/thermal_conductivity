=======
History
=======
2026.10.14 -- Bugfix: crash with coarse sampling of the heat flux
    * The Green–Kubo fit failed when the heat flux was sampled less often than every
      0.4 ps. It no longer fails.
    * The analysis now uses the shared ``seamm_md_analysis`` library, also used by the
      Diffusivity step, instead of its own copy. The results are the same, and the
      Helfand moments are computed much faster (by FFT): for a million samples, 0.4 s
      instead of about 2 minutes.
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
