## MAKE A POST ABOUT THIS!!!

#   Script to clean up PDF scrape of Sarawak Gazetteer
#   Most of cleaning was performed using specifications in tabula:
#       system("java -jar tabula-1.0.2-jar-with-dependencies.jar -t --area %10,0,85,100 --columns 245,320,435 --pages 16-100 ~/Documents/GradSchool/Research/CircuitTheory_Borneo/distanceToPopulation/Revision2_PopulationSources/Sarawak_Gazetteer.pdf -o ~/Documents/GradSchool/Research/CircuitTheory_Borneo/distanceToPopulation/Revision2_PopulationSources/Sarawak_Gazetteer_tabula_allpages.csv")

# Created by MDeith on June 9, 2019; last edited June 9, 2019

setwd("~/Documents/GradSchool/Research/CircuitTheory_Borneo/distanceToPopulation/Revision2_PopulationSources/")
o_csv <- read.csv("Sarawak_Gazetteer_tabula_allpages.csv", skip = 15, colClasses = "character", na.strings = "")

head(o_csv)
for(i in 1:nrow(o_csv)){
  if(sum())
}