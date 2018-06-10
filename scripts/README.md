# Scripts undertaking power analysis to detect difference in correlations in MZ and DZ twin groups #

## Conduct simulations
### Main analysis script
corr_power_simulation_analysis.r
### Config file (requires set up with specification of SQL connection detail, if using sql related functions)
config.yml.README.txt	

## Prepare plots and associated analyses (e.g. spline interpolation to estimate power given sample size)
corr_power_plots.R	

## Additional
### Early script testing time of various functions
corr_timetest.r	

### C++ file with alternate routines for random number draws (for RCCP gtv function; thanks to Koen Simons)
test.cpp	       

### A stata do file with early notes experimenting with correlation formulae using mata
stata_notes.do.txt	

Carl Higgs 11 June 2018
carlhiggs@gmail.com