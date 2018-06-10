# Scripts undertaking power analysis to detect difference in correlations in MZ and DZ twin groups #

## Conduct simulations
### Main analysis script
corr_power_simulation_analysis.r
### Config file (requires set up with specification of SQL connection detail, if using sql related functions)
config.yml.README.txt	

## Prepare plots and associated analyses (e.g. spline interpolation to estimate power given sample size)
corr_power_plots.R	

## Interactive power calculator app
A prototype of an interactive calculator of power for detecting difference in correlations.  Current 'proof of concept' approach uses Fisher's Z analytic formula only; however, this can be expanded to draw upon underlying database of pre-processed results and display interpolated plot types as included in report.  The intent is that a researcher planning a study for which they will use comparison of Pearson correlations can specify their key parameters (hypothesised population correlations in each group, distributions, sample sizes etc) and interactively to determine power given parameters to inform their approach.
https://anaestheteick.shinyapps.io/corr_power_app/

server.R
ui.R

## Additional
### Early script testing time of various functions
corr_timetest.r	

### C++ file with alternate routines for random number draws (for RCCP gtv function; thanks to Koen Simons)
test.cpp	       

### A stata do file with early notes experimenting with correlation formulae using mata
stata_notes.do.txt	

Carl Higgs 11 June 2018
carlhiggs@gmail.com