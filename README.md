# Myosin-cluster-analysis
This Matlab script is used to quantitatively analyze myosin clusters in cells. Myosin clusters were endogenously tagged and imaged with spinning-disk confocal microscopy. The manuscript is in preparation.

myosinPeakInt_main_v2.m is the main script that runs the analysis. The workflow is as follows: subtract background, feature identification, feature filtering, intensity fitting. After identifying myosin clusters and calculating their intensity, the results are saved and can be loaded into myosinPeakInt_analysis.m for further analysis.

This analysis routine uses the particle identification code based on John Crocker and David Grier's collidal particle tracking [1]. The Matlab implementation of their work is done by Maria Kilfoil.

Tested on Matlab 2021a and 2021b. Requires the Image Processing Toolbox.

Refernce: 
[1] J. C. Crocker and D. G. Grier, "Methods of digital video microscopy for colloidal studies," Journal of Colloid and Interface Science 179, 298-310 (1996).
