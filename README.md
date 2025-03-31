This repository contains code used in article "Isotropy testing in spatial point patterns: nonparametric versus parametric replication under misspecification" by J.J Pypkowski, A.M. Sykulski, and J.S. Martin (https://doi.org/10.48550/arXiv.2411.19633).

Folder "functions" contains the most up-to-date versions of functions used in our study. Files in this folder also contian comments that help understand how the functions work. 

Folders "GIBBS", "LGCP", and "PLCP" contain scripts used to run the simulation study in our article; subfolders "short" contain the scripts' counterparts for smaller r_max=0.05 (and both smaller and larger r_max when replication using Thomas process was used). These scripts take a very long time to run! If you want to test them quickly, remember to reduce the number of patterns, number of each pattern's replicates, and, if applicable, the number of iterations for algorihms such as stochastic reconstruction or inference for PLCP. Folder "GIBBS" also contains a script "Gibbs_big_2.R" simulating patterns from Gibbs processes used in the article. The output of this file is saved in .csv format. The output of the remianing scripts in the three folders is not saved in this repository due to memory constraints. 

Files "output_process.R" and "output_process_short.R" in folder "plotting" process the output files of the simulation scripts, and produce rejection rates saved in the .csv format folder "plotting". File "plotting_paper.R" produces plots presented in the article. 

Folder "ambrosia dumosa - application" contains scripts used in the analysis of the Ambrosia Dumosa data in the final section of our article. The data can be found online at [https://doi.org/10.5063/AA/connolly.205.1](url). The outputs of all scripts are also saved in the folder.

Lastly, file "test_your_pattern.R" can be used by users who wish to apply our test (with either tiling or stochastic reconstruction) to their data. Please follow the instructions in the script to run the code. 

Please contact the corresponding author of the article if you have any questions regarding the code.
