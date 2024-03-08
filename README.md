# dinosaur_diversification
This repository contains the data and code associated with Allen BJ, Volkova Oliveira MV, Stadler T, Vaughan TG, Warnock RCM (2024) Mechanistic phylodynamic models do not provide conclusive evidence that dinosaurs were in decline before their final extinction, Cambridge Prisms: Extinction, in press.

/data/change_times.txt
Describes the break points (in time) used in the piecewise constant analyses.

/data/dinosaur_ages.csv
The age range data for all dinosaurs, downloaded from the Paleobiology Database.

/data/tip_constraints.csv
The cleaned age range data corresponding to tips in the phylogenies.

/R code/BDSKY_post_processing.R
Code to process and plot the log data from the BDSKY analyses.

/R code/Coalescent_post_processing.R
Code to process and plot the log data from the coalescent analyses.

/R code/PBDB_tip_constraints.R
Code to convert the age data from the Paleobiology Database into a suitable format for setting tip constraints in BEAST2.

/R code/Wrange_trees.R
Code to clean, split, wrangle and plot the phylogenies used, ready for input into BEAST2.

/trees/*.tree
Files containing the cleaned, fixed tree topologies inputted into each BEAST2 analysis.

/trees/Tree modifications log.xslx
A list of all modifications made to the phylogenies used, as taken from the supplementary files of their respective papers.

/XML/BDSKY.xml
The XML file used to run the fossilised birth-death skyline analyses in BEAST2.

/XML/PiecewiseCoalescent.xml
The XML file used to run the piecewise coalescent analyses in BEAST2.
