# AH-Metabolic-Modeling

SIMPLY DOWNLOAD FILE, OPEN SCRIPT, AND HIT RUN!

# ABSTRACT:
Acute-on-chronic diseases are marked by progressive, often clinically silent, deterioration bifurcated by sudden shifts in compensatory state. Predicting this decompensation is a critical clinical need, as it enables intervention when treatment tends to be most effective. Alcohol-associated liver disease (ALD) is one such disease: its acute decompensation, alcohol-associated hepatitis (AH), carries a short-term mortality of 30-50%. To understand how the liver approaches this threshold, we build on metabolic control analysis (MCA) and extend it from single pathways to interacting metabolic modules, yielding a reduced order, flux-centric framework in which module velocities, rather than metabolite concentrations, are the dynamical variables. A fixed total velocity is shared across modules through a stoichiometry derived coupling network: steady state capacity fractions play the role of control coefficients, and metabolite variables act as elasticity load signals, so perturbation cannot create capacity, only reallocate it. This sidesteps the parameter nonidentifiability and numerical stiffness that limit kinetic models at organ scale. We specify the framework for hepatic metabolism in AH, coarse-graining it into glutathione synthesis, the methionine cycle, the urea cycle, and fatty acid biosynthesis, with baseline fluxes inferred from 89 patient liver transcriptomes. Canonical ethanol phenotypes (reduced methylation, suppressed urea flux, and enhanced lipogenesis) emerge from resource competition in 82-$94% of simulated patients and degrade when inter-module coupling is removed. Repeated exposure drives a qualitative transition toward joint oxidative and metabolic failure rather than gradual decline. Finally, simulated corticosteroid treatment restores stability in 72-82% of severe AH trajectories, bracketing the clinical responder rate of ~75\% without calibration, with non-response determined by failure class and treatment timing.

# FILES:
(1) Validation: Contains all files necessary to run initial concordance search and mouse model experiment. Outputs included.
(2) Data: contains script to construct adjacency matrix (with necessary associated file) along with script needed to compute rate estimations of flux dynamics
(3) UQ: Contains all files necessary to run UQ of one-hit and two-hit ethanol perturbations. Outputs included.
(4) Treatment: Contains all files necessary to simulate and plot corticosteroid treatment on the system. Outputs included.
(5) MultiHit: Contains all files necessary to run the one to four-hit ethanol perturbation simulations. Outputs included. 

