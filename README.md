# weighted_propscores
Code for the simulation study and applied example from "Propensity Scores in Convenience Samples."
 
## Simulation Study from Section 3

Code to generate the plots for Figures 1, 2, 3, 4, and B1 is in plots_for_simstudy.R
The results needed for these plots are generated with the following files:

+ simstudy_vary_corr.R
+ simstudy_vary_T.R
+ simstudy_vary_TK.R
+ simstudy_vary_X2.R
+ simstudy_vary_X3.R
+ simstudy_vary_X4.R
+ simstudy_vary_X5.R
+ simstudy_vary_propmodel.R

You can run these files to generate the results and plots_for_study.R 
has comments telling you which results are needed to create each plot.

The results and plot for Figure 5 are created with simstudy_withSE.R

## Applied Example from Section 4

Code to clean the NACC data is in applied_example_clean_NACC.R and code to create Table 1 and Figure 6 is in applied_example_run_analysis.R

Note that you will need to use NHANES as a representative dataset and code to clean that data is in the estweight_simulationstudy repository. NACC data needs to be requested from https://naccdata.org/
