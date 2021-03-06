# Scripts undertaking power analysis to detect difference in correlations in MZ and DZ twin groups #
https://bitbucket.org/carlhiggs/bca_rp2/src/master/scripts/

## Conduct simulations
### Main analysis script
https://bitbucket.org/carlhiggs/bca_rp2/src/master/scripts/corr_power_simulation_analysis.r
### Config file (requires set up with specification of SQL connection detail, if using sql related functions)

https://bitbucket.org/carlhiggs/bca_rp2/src/master/scripts/config.yml.README.txt
## Prepare plots and associated analyses (e.g. spline interpolation to estimate power given sample size)
https://bitbucket.org/carlhiggs/bca_rp2/src/master/scripts/corr_power_plots.R

## Interactive power calculator app
A prototype of an interactive calculator of power for detecting difference in correlations.  

See: https://anaestheteick.shinyapps.io/corr_power_app/

Current 'proof of concept' approach uses Fisher's Z analytic formula only; however, this can be expanded to draw upon underlying database of pre-processed results and display interpolated plot types as included in report.  The intent is that a researcher planning a study for which they will use comparison of Pearson correlations can specify their key parameters (hypothesised population correlations in each group, distributions, sample sizes etc) and interactively to determine power given parameters to inform their approach.


https://bitbucket.org/carlhiggs/bca_rp2/src/master/scripts/corr_power_app/server.R

https://bitbucket.org/carlhiggs/bca_rp2/src/master/scripts/corr_power_app/ui.R

## Additional
### Early script testing time to process various functions
https://bitbucket.org/carlhiggs/bca_rp2/src/master/scripts/corr_timetest.r

### C++ file with alternate routines for random number draws (for RCCP gtv function; thanks to Koen Simons)
https://bitbucket.org/carlhiggs/bca_rp2/src/master/scripts/test.cpp

### A Stata do file with early notes experimenting with correlation formulae related to twin data using mata
https://bitbucket.org/carlhiggs/bca_rp2/src/master/scripts/stata_notes.do.txt

Carl Higgs 11 June 2018
carlhiggs@gmail.com