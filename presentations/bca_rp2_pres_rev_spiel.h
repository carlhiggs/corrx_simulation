1. Research question
Through my research project, I planned and carried out a power analysis of the detection of differences in correlations between identical and non-identical twin pairs, as a pre-cursor step to undertaking a classic twin study.

The classic twin study exploits the differing degrees of genetic relatedness in identical twins and non-identical twins in order to draw inferences on the heritability of traits

In broad terms, heritability is the degree to which variation in a trait or phenotype, such as propensity to gain body weight, or become a centenarian, can be attributed to genetics.

Identical twins arise from the same zygote, or fertilised egg, and are genetically very similar.  Non-identical twins arise from fertilisation of two seperate eggs, and are as genetically alike as ordinary siblings.  

2. 
We might assume that variation in a trait occuring in identical twins must be due to environmental, not genetic factors. 
If we make a some more assumptions, understanding that these likely do not strictly hold in practice we can undertake a classic twin study.  

3. 
the comparison of phenotypic traits within identical and non-identical twin pair samples allows for partitioning the variance in traits into that attributable to  shared environment, individual environment or to genetics.  This allows us to for example have a better understanding of the mechanics of health and disease processes so that we can develop intervention measures which appropriately target the hypothesised causal mechanisms.
[see Falconer 183, 184]
[1 min 15 Secs]
You will have noticed the clustered data structure here. 
Contemporary twin studies use mixed effects and structural equation modelling to evaluate these clustered differences in variance components for a particular trait between mono and dizygotic twins.  A recent article reported on the development of R commands which can be used to estimate the power to detect a difference in these.  

However, a routine preliminary step undertaken by researchers in this field is the calculation and comparison of the more basic Pearson correlations as a precursor to considering intraclass correlations.  This research project focuses on the reported needs of researchers undertaking this early analysis step to consider the power to detect differences in Pearson correlations for bivariate normal and non-normally distributed data between mono and dizygotic twins.

4. A skip through the literature  
The statistical treatment of correlation was popularised by Francis Galton.  In context of broad social interest in eugenics, Galton described methods which could be used to describe the co-relatededness of variables sourced from closely related family members.  

Karl Pearson was a keen follower of Galton's research, and writing in the context of inheritance and natural selection in 1895 makes reference to what he called Galton's function, describing it as a coefficient of correlation.  

William Gossett and other big names in statistics made contributions to the understanding of what we now know as Pearson's correlation. 

Ronald Fisher observed that a geometric transformation of the correlation cofficient using its inverse hyperbolic tangent could be used to approximate a normal distribution.  This has the effect of mapping the distribution of correlation of coefficients from a domain of negative 1 through positive one to negative infinity through positive infinity. Being a simple and accurate approximation of the normal distribution, this transformation known as Fisher's Z is ubiquitous in statistical treatment of correlation coefficients, for example when seeking to compare their differences.

Florence Nightingale David was a protege of Pearson's who suggested that she prepare a volume of numerically accurate tables and interpolated plots of the distribution of the sample correlation coefficient r given n and the population correlation rho which could act as a standard against which to judge approximations such as that of Fisher.  This visualisation of David's of the chance of rejecting the null hypothesis when true given alpha and rho and n was an influence on this project's presentation of results.

In the context of genetics and kinship, Douglas Falconer outlined methods of using comparison of identical and non-identical twins for estimating heritability, and noted some of the assumptions that this involves.

Michael Neale and colleagues developed methods and software to facilitate the analysis of variance components in twin studies using structural equation modelling.  Brad Verhulst, building on the work of Peter Visscher, developed functions for power analyses in this variance component modelling context, for example detecting a difference in genetic correlations.

This project is concerned with power to detect a difference in Pearson correlations as part of preliminary analysis before variance component modelling.

5. Power for difference in correlations 
The consideration of the how large an effect size must be in order to be reliably detected given sample sizes is a key step in planning a study and developing a statistical analysis plan.  This involves a compromise between type 1 and type 2 error thresholds, respectively the expected proportion of null hypotheses to be rejected when true and not rejected when false.  These could be chosen to suit the requirements of a particular study, but for historical reasons the usuals consensus is for 5% and 20%, and this is the parameterisation adopted for the results presented here.Power is the complement of type 2 error.

A test statistic for the difference in correlations can be calculated as the difference in Fisher's Z transformed values divded by the approximate standard error of the difference.  This is a classic formulation, present in David's 1938 book and still used in programs in Stata and R; assuming approximate bivariate normality it is regarded as a very good approximation.

If you take the difference between a normal reference score such as 1.96 --which approximately corresponds to our alpha level of 0.05-- and the absolute value of the test statistic, and then find the probability of observing a value of at least this magnitude on the normal distribution; 1 minus this value is the power.

Using a simulation approach we take our hypothesis tests and apply them to draws from simulated data designed to mimic our samples through parameterisation using the hypothesised underlying bivariate population distributions.  

So where in our formula we might plug in hypothesised sample coefficients of 0.2 and 0.5, in the simulation we use these values to represent the true correlations in the underlying population from which we draw our samples.  Over a large number of simulations of bivariate mz and dz twin data the proportion of hypothesis tests returning p-values lower than our type 1 error threshold is our power estimate.

6. Simulations
The tests we identified and developed for inclusion were 

The simulation equivalent of the Fisher's Z test already discussed
The GVT test involves transforms the simulated sample correlations into gneralised pivotal quantities the difference of which is then use to calculate a test statistic and p-value, which should be robust to departures from normality.

The signed log likelihood ratio is asymptotically normally distributed and a test is formulated as the signed difference in sample correlation coefficients multiplied by the square root of the sum of respective coefficients log-likelihoods.

The permutation test is a non-parametric approach which compares the absolute difference of the Z-transformed sample correlations with correlations from a series of group members ship permuatations using the sample rank orders as values.

Zou's confidence interval approach evaluates whether zero lies within the lower and upper bounds of the interval estimate of the difference in correlations, returning 1 if so or otherwise zero.  Over a run of simulations this would be expected to return identical results to the Fisher Z test, which it is theoretically identical to.  The main reason for inclusion is how much faster it is in doing this.

7. Preliminary results: simulation times
Reporting on the simulation efficiency investigation which led to revised plan with coarser resolution subset of parameters, and dropping permutation test --- the final line is what I actually did, the line above it is what I originally planned (and revised, for obvious reasons!)

8. Main findings (and slide 9; its really like 8a and 8b...)
I will point out that a summary table of average performance of tests is not so informative (although I could provide this, and will in report).  Instead I provide two examples of how the simulation architecture we have developed can be applied for specific questions --- how will such and such a test perform given parameterisation (distribution, N, MzDz ratio, etc); or, given a fixed parameterisation set allowing N to vary how do the various test approaches differ in recommended sample size?

The images here are developed using a monotonic increasing spline function; I trialled a number of approaches (e.g. prediction using polynomial regression, loess curves, etc) and this was most faithful, plausible and reliable across the range of parameters.  Surprisingly nice fit, given 7 points per test.  (this is a limitation really).

Anyway, here I explain the first plot as inspired by the FN David images displayed earlier; however, the right image (inspired by Stata) is more applied and useful for comparison across tests.  The Stata reference is also deliberate as I will contrast the Fisher Z formula performance (used in Stata two correlation power calculator) with next image using skewed data to illustrate the impact of extreme data irregularities (yes, in normal practice we would transform our data to better meet assumptions; but what would happen if we didn't? some suggestion here).

The SLR test appears a strong performer - but its systematic elevation at tails of distribution is concerning -- is this a systematic bias, and perhaps innacurate rather than truly improved performance?  I need to think and investigate this more as I think its a surprising results.

10. Limitations 
Well quite a few, I'm happy to be honest here (its not end of semester after all).  The project did take me longer than I anticipated, but these things happen.  There is more we can do, and I have some ideas about how these things can be done.

11. Sum up and next steps
I'll talk about how we have created both the architecture for a process, as well as a database of simulation scenarios that can be interrogated.  Both can be expanded as required.  I'll mention the sketch of a web app we have developed, and intent to incorporate lookup of the pre-processed database into this to allow comparison of tests using archived simulation results on the fly.