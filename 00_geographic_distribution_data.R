## Program to load and generate species distribution data
## Gabriel Munoz
## Functional Asymmetry in Ecological Networks
## Nov. 2024 


## Dependencies 


library(tidyverse) # for data manipulation 
library(sf) # for simple feature manipulation 

# assum planar geodesics (to fix double vertices issue with range polygons)
sf::sf_use_s2(FALSE)


# load data directory
palm_all_files <- list.files("00_Data/00_species_distribution/Palm-distribution-ranges/Shapefiles/", 
full.names = T)
# get  shp files
palm_shp_files <- palm_all_files[str_detect(palm_all_files, ".shp")]
# filter xml out
palm_shp_files <- palm_shp_files[!str_detect(palm_shp_files, ".xml")]


# get map of biogeographic dominions on the Neotropics 

neotropics <- st_read('00_Data/03_Landscape/Morrone_Neotropics/Lowenberg_Neto_2014.shp')

grid <- st_make_grid(neotropics, cellsize = c(1, 1), what = "polygons", 
                     crs = sf::st_crs(st_read(palm_shp_files[1])))
# Convert the grid to a simple feature collection
grid <- st_sf(grid)


# load data directory
all_files <- list.files("00_Data/00_species_distribution/TERRESTRIAL_MAMMALS/", full.names = T)
# get  shp files
shp_files <- all_files[str_detect(all_files, ".shp")]
# filter xml out
shp_files <- shp_files[!str_detect(shp_files, ".xml")]



# read polygons for palms 


palm_ranges <- palm_shp_files %>% 
  purrr::map(function(shapefile){
    
    st_read(shapefile)
  })


# define function to grid palm ranges 
get_grid_palms <- function(shapefile, grid){
 
    # Intersect the shapefile with the grid
  result <- st_intersection(shapefile, grid)

  
  return(result)
  
}

# get a safe version of the function 

get_grid_palms_safe <- safely(get_grid_palms)

# apply the function for all palm species 
palm_grids <- palm_ranges  %>% 
  purrr::map(function(shapefiles){
    
    get_grid_palms_safe(shapefiles, grid)
    
    
    
  })

saveRDS(palm_grids, "00_Data/00_species_distribution/gridded_palm_data.RDS")


# read mammal shapes 
mamm <- st_read(shp_files)

grid <- st_make_grid(neotropics, cellsize = c(1, 1), what = "polygons", 
                     crs = sf::st_crs(mamm))
# Convert the grid to a simple feature collection
grid <- st_sf(grid)


lobstr::tree(mamm$id_no)
# read in trait data 
mammal_traits <- read.csv("00_Data/01_species_traits/final_mammal_trait.csv")

# filter range data to those species with trait data (frugivore species)
mamm1 <- mamm %>%
  filter(mamm$binomial %in% mammal_traits$Scientific)
# Crop the resulting data for the neotropics
mamm2 <- st_crop(mamm1, grid)


# mammal grids 
mammal_grids <- 1:nrow(mamm2) %>% 
  map(function(species){
    
    # apply function to grid mammals 
 st_intersection(mamm2$geometry[species], grid)
    
  })

mammal_grids <- mammal_grids %>% set_names(mamm2$binomial)

saveRDS(mammal_grids, "00_Data/00_species_distribution/gridded_mammal_data.RDS")

