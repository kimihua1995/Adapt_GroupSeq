# Adapt_GroupSeq
R code for the article "Multiplicity Control in Clinical Trials with Adaptive Selection Followed by Group-Sequential Testing", including:

1) A mock example, `Mock Example.Rmd`, to illustrate the step-by-step calculation of the group-sequential p-value;
2) All codes of the simulation studies.

The main code is in `sim.R`, with 5 source function codes:

`source_combine.R`: includes the main function for the simulation with data generation, test statistics calculation, implement of the proposed method and the comparater method (Sugitani et al., 2016).

`source_datsim.R`: functions related to data generation.

`source_proposed_design.R`: functions related to the proposed method.

`source_graphical.R`: functions related to the comparater method (Sugitani et al., 2016).

`source_other.R`: functions including calculating test statistics, correlation matrix, get the IA times, etc.

In addition, the code `information_derivation.R` shows how to derive the information fraction based on observed number of events using simulation method.

