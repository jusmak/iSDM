survey_sites <- read.csv('/project/fas/jetz/data/diego/Hummingbirds/Survey_Data/code/Juan_parra_checklists/Sites_8Feb2011.csv')
survey_species <- read.csv('/project/fas/jetz/data/diego/Hummingbirds/Survey_Data/code/Juan_parra_checklists/SpeciesxSite8Feb2011.csv')
survey_species= survey_species[,c(1:2,5)]
survey_jp = left_join(survey_sites, survey_species)