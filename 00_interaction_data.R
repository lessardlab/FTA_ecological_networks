## Program to read and prepare interaction data
## Functional Assymmetry in Ecological Networks 
## Gabriel Munoz
## Nov. 2024


int_data <- read.csv("00_Data/02_species_interactions/PalmDryadRepo/Munozetal2019/PalmFrugDatasetOCT2018.csv")

# filter to include mammal species from original dataset


int_data <- int_data %>% 
  filter(frugClass == 'MAMMAL', 
         biogeographicRegion == 'Neotropics')
  
unique(int_data$PALM)
dim(int_data)

unique(int_data$FRUGIVORE) %>% length()

saveRDS(int_data, '00_Data/02_species_interactions/final_int_data.RDS')
