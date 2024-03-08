# Bethany Allen (9th Jan 2023)
# Code to clean PBDB age ranges into tip constraints

library(tidyverse)

#Read in occurrence data
occurrences <- read_csv("trees/dinosaur_ages.csv")

#Add underscore to taxon names
occurrences$taxon_name <- gsub(" ", "_", occurrences$taxon_name)

#Remove taxa without age data
no_NA <- filter(occurrences, !is.na(firstapp_max_ma))

#Transform full age range relative to K-Pg boundary
no_NA$firstapp_max_ma <- no_NA$firstapp_max_ma - 66.0
no_NA$lastapp_min_ma <- no_NA$lastapp_min_ma - 66.0

#Write as tip constraints
write_csv(no_NA, "trees/tip_constraints.csv")
