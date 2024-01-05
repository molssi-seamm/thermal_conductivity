=======
History
=======
2024.1.5 -- Bugfix: Thermal conductivity
    * If the Helfand moments fit in the thermal conductivity step failed it stopped the
      entire job. This is correct, as well as some of the undelying causes for
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
