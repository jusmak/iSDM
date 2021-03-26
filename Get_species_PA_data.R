wd <- 'C:/Users/OMISTAJA/OneDrive - Yale University/Hummingbird_project'
prevalence_temp <- .02
library(dplyr)
survey_sites <- read.csv(paste(wd, '/Data/PA_observations/Juan_parra_checklists/Sites_8Feb2011.csv', sep = ''))
survey_species <- read.csv(paste(wd, '/Data/PA_observations/Juan_parra_checklists/SpeciesxSite8Feb2011.csv', sep = ''))
survey_species= survey_species[,c(1:2,5)]
survey_jp = left_join(survey_sites, survey_species)

#check prevalence for each species
sp_unique <- matrix(unique(survey_species$Spname))
sp_prevalence <- matrix(apply(sp_unique, 1, function(x) sum(survey_jp$Spname==x)/nrow(survey_sites)))

#print species into a csv
sp_prevalent <- sp_unique[sp_prevalence>prevalence_temp]
write.table(sp_prevalent, file = paste(wd, '/Data/PA_observations/Sp_prevalent_2.txt', sep = ''),
          row.names = FALSE, col.names = FALSE)
