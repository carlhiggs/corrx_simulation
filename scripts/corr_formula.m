% Project: PowerMzDzDiffR - bca_rp2
% Date: 28 Feb 2018
% Author:  Carl Higgs
% Purpose: Sketch file for matlab code to get a feel for the language


% Inferences about \rho\hat based on more than one sample
% From: Burkett - Correlation_Regression_and_Analysis_of_Variance.
% http://www.academia.edu/34260699/Correlation_Regression_and_Analysis_of_Variance

r = [.7  .9]
n = [19 25]
z = atanh(r)
meanzeta = z - 5*r/(2*n)
varzeta = (n - 3/2 + (5/2)*(1-r.^2)).^-1
postvz = (1/varzeta(1) + 1/varzeta(2))^-1
postmz = postvz*sum(meanzeta./varzeta)
rhohat = tanh(postmz)
rhohat

ci95zeta = [postmz - 1.96*sqrt(postvz)  postmz + 1.96*sqrt(postvz)]
ci95rho  = [tanh(ci95zeta(1)) tanh(ci95zeta(2))]
ci95rho