## Program to prepare data on species traits
## Gabriel Munoz
## Functional Asymmetry in Ecological Networks
## Nov 2024


## Load and subset trait data 
library(dplyr)
library(ggplot2)

# Load data on palm traits 1.0 

palm_trait_data <- read.csv("00_Data/01_species_traits/01_palms/PalmTraits_1.0.txt", 
                            sep = '\t')

# Load data from Elton database 

mamm_traits <- read.csv("00_Data/01_species_traits/02_mammals/Elton_Traits/Elton Traits/MamFuncDat.txt", 
                        sep = '\t')
# Load data from Body size 

mamma_body_mass <- read.csv("00_Data/01_species_traits/02_mammals/Mammal_BodyMass-2023/Mammal_BodyMass/MammalMassSandom2013.csv")

# subset to interaction relevant traits 

## Palms 

palm_trait_data <- palm_trait_data %>% 
  select(SpecName, PalmTribe, accGenus, accSpecies, PalmTribe, Acaulescent, Erect,
         MaxStemHeight_m, AverageFruitLength_cm) %>% 
  filter(Acaulescent != 2) %>% 
  mutate(across(all_of(matches(c('MaxStemHeight_m',
                                 'AverageFruitLength_cm'))), 
                ~log1p(.)
  ))


## Mammals 


mamm_traits <- mamm_traits %>% 
  select(Scientific,MSWFamilyLatin, Diet.Inv,Diet.Vend, Diet.Vect, Diet.Vfish, 
         Diet.Vunk, Diet.Scav, Diet.Fruit,Diet.Nect, Diet.Seed, Diet.PlantO,
         Activity.Nocturnal,Activity.Crepuscular, Activity.Diurnal, BodyMass.Value) %>% 
  filter(Diet.Fruit != 0, 
         MSWFamilyLatin != 'Phyllostomidae') %>% 
  mutate(BodyMass.Value = log1p(BodyMass.Value))


## write final trait data 

## Palms 

write.csv(mamm_traits,
          "00_Data/01_species_traits/final_mammal_trait.csv", 
          row.names = F)
## Mammals 

write.csv(palm_trait_data,
          "00_Data/01_species_traits/final_palm_trait.csv", 
          row.names = F)
