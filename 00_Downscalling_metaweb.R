## Implementation program: Functional Asymmetry in Ecological Networs
## Gabriel Munoz 

## As per the article draft:
## Dependencies 

library(tidyverse)
library(sf)
library(raster)
source("00_functions_nov24.R")
## Methods:


### Geographic distribution data

# get grid level assemblages for palms and mammals

palm_grids <-
    readRDS("/00_Data/00_species_distribution/gridded_palm_data.RDS")
mammal_grids <-
    readRDS("00_Data/00_species_distribution/gridded_mammal_data.RDS")

### Trait data

palm_trait_data <- read.csv('00_Data/01_species_traits/final_palm_trait.csv')
mammal_trait_data <- read.csv('00_Data/01_species_traits/final_mammal_trait.csv')

### Interaction data

# get available interaction data for the Neotropics
int_data <- readRDS('00_Data/02_species_interactions/final_int_data.RDS')
palm_trait_data <- palm_trait_data %>% 
  filter(palm_trait_data$SpecName %in% int_data$PALM )
mammal_trait_data <- mammal_trait_data %>% 
  filter(mammal_trait_data$Scientific %in% int_data$FRUGIVORE )
# filter to match species between databasees
int_data <- int_data %>% 
  filter(int_data$PALM %in% palm_trait_data$SpecName,
         int_data$FRUGIVORE %in% mammal_trait_data$Scientific)


### Species pool 

neotropics <- st_read('00_Data/03_Landscape/Morrone_Neotropics/Lowenberg_Neto_2014.shp')
grid <- st_make_grid(neotropics, cellsize = c(1, 1), what = "polygons")
# Convert the grid to a simple feature collection
grid <- st_sf(grid)

### Environmental data

# download climatic data
WCLim <- raster::getData("worldclim", var="bio",res=10)

cropMask <- function(raster,prov){
  ## crop and mask
  r2 <- crop(raster, extent(prov))
  r3 <- mask(r2, prov)
  return(r3)
}

# crop data to the neotropics
WCLim <- cropMask(WCLim, neotropics)

# Separate variables of interest

Temp <- WCLim[[1]]
Prec <- WCLim[[12]]
PrecSe <- WCLim[[15]]
IsoTer<- WCLim[[3]]
TempSeaso<- WCLim[[14]]

Temp <- aggregate(Temp, 1/0.17)
Prec <- aggregate(Prec, 1/0.17)
PrecSe <- aggregate(PrecSe, 1/0.17)
TempSeaso <- aggregate(TempSeaso, 1/0.17)


## Statistical analyses

##=========================================================================
# Building a probabilistic metaweb from aggregated binary observations 
##=========================================================================


## Step 1: make binary matrix
N <- int_data %>% 
  xtabs(~PALM + FRUGIVORE, .)
N[N>1] <- 1

# Step 2:  get L matrix

pTRLQ <- palm_trait_data %>% 
  column_to_rownames("SpecName") %>% 
  dplyr::select(Acaulescent, Erect, MaxStemHeight_m, AverageFruitLength_cm)
  
# Step 3: get Q matrix
mTRLQ <- mammal_trait_data %>% 
  column_to_rownames("Scientific") %>%
  dplyr::select(!MSWFamilyLatin)

## Step 4: Reorder the matrix 
nRLQ <- N[rownames(pTRLQ),rownames(mTRLQ)]

## Step 5: prepare object to fit latent trait models 
Ng <- cassandRa::CreateListObject(nRLQ)

## Step 6: fit all models (time consuming)
# latent_network_models <- cassandRa::FitAllModels(Ng)

## Read in fitted result 
latent_network_models <- readRDS('C:/Users/gabri/Documents/PhD/00_Chapter_palms_mammal_interactions/R-analisis/00_Data/04_models/latent_net_mod.RDS')

## Step 7. Compare models with Youden's J 

### list all models 
probNet <- list(latent_network_models$SBM_ProbsMat, 
                latent_network_models$C_ProbsMatrix, 
                latent_network_models$M_ProbsMatrix, 
                latent_network_models$B_ProbsMat)


TestYJ(probNet[[1]],latent_network_models$obs, 2)
## Apply the to TestYJ function to compute the youdens J statistic
YJtestin <- lapply(1:4, function(i) TestYJ(probNet[[i]],latent_network_models$obs, 100))

## Rearrange the resulting dataset and name variables appropiatetly

YJtestin <- YJtestin |>
  set_names(c('SBM', 'Cent', 'Matc', 'Match_Cent')) |>
  imap(~{
    .x |>
      mutate(id = .y )
  }) |>
  bind_rows()


## Compare aggregate measures across models 
YJtestin_agg <- 
YJtestin %>%
  group_by(id) %>% 
  summarize(sens = median(sens), 
            speci = median(speci), 
            yj = median(YJ))


# Refit SBM
SBMs <- cassandRa::FitSBM(Ng)

## Fitted group associations for palms 
PalmNet <- data.frame(Ng$HostNames,
                      "SBMs.SB_H" = SBMs$SBM1$SB_H)

## Fitted group associations for mammals 
MammNet <- data.frame(Ng$WaspNames,
                      "SBMs.SB_W" = SBMs$SBM1$SB_W)

PalmNet <- data.frame(PalmNet,
                      palm_trait_data[match(PalmNet$Ng.HostNames,
                                        palm_trait_data$SpecName),])
MammNet <- data.frame(MammNet,
                      mammal_trait_data[match(MammNet$Ng.WaspNames,
                                          mammal_trait_data$Scientific),])


## Downscaling the continental metaweb to generate grid-cell level networks


## Step 1: multinomial logistic regression models aimed to predict the species level SBM model results


library(tidymodels)

MammNet <- MammNet %>% 
  mutate(SBMs.SB_W = as.factor(SBMs.SB_W))
{
  # Split the data into training and testing sets 
  
  data_split <- initial_split(MammNet, prop = 0.80)  # 75% of the data goes to the training set
  
  # Extract the training set
  data_train <- training(data_split)
  
  # Extract the testing set
  data_test <- testing(data_split)
  
  ## Define  a list of recipes 
  
  rec_list <- list(
    "recipe1" = recipe(SBMs.SB_W ~  BodyMass.Value + Diet.Fruit, data = MammNet) %>%
      step_log(BodyMass.Value, base = 10),
    "recipe2" = recipe(SBMs.SB_W ~  BodyMass.Value + Diet.Fruit, data = MammNet),
    "recipe3" = recipe(SBMs.SB_W ~ BodyMass.Value + Diet.Fruit, data = MammNet) %>%
      step_log(BodyMass.Value, base = 10),
    "recipe4" = recipe(SBMs.SB_W ~ BodyMass.Value + Diet.Fruit, data = MammNet)
  )
  
  
  # Specify model using parsnip
  model_spec <- multinom_reg() %>%
    set_engine("nnet") %>%
    set_mode("classification")
  
  
  # Create a list of workflows using purrr::map
  workflows <- map(rec_list, ~workflow() %>%
                     add_recipe(.x) %>%
                     add_model(model_spec) %>%
                     fit(data_train))
  }



# Calculate the AUC for each model
aucs_mam <-  map_df(workflows, ~.x %>%
                      augment(data_test) %>%
                      roc_auc(truth = SBMs.SB_W, .pred_1:.pred_7))



# Augment the workflows and calculate the ROC data
roc_data_mammal <- map(workflows, ~.x %>%
                  augment(data_test) %>%
                  roc_curve(truth = SBMs.SB_W, .pred_1:.pred_7))

# Combine the ROC data into a single data frame
roc_data_mammal_combined <- bind_rows(roc_data_mammal, .id = "Model")



### Repeat for palms 


PalmNet <- PalmNet %>% 
  mutate(SBMs.SB_H = as.factor(SBMs.SB_H))

# Split the data into training and testing sets 
{
  data_split <- initial_split(PalmNet, prop = 0.80)  
  # 75% of the data goes to the training set
  
  # Extract the training set
  data_train <- training(data_split)
  
  # Extract the testing set
  data_test <- testing(data_split)
  
  ## Define  a list of recipes 
  
  rec_list <- list(
    "recipe1" = recipe(SBMs.SB_H ~  MaxStemHeight_m + AverageFruitLength_cm, data = data_train) %>%
      step_log(AverageFruitLength_cm, base = 10, offset = 1),
    "recipe2" = recipe(SBMs.SB_H ~  MaxStemHeight_m + AverageFruitLength_cm, data = data_train),
    "recipe3" = recipe(SBMs.SB_H ~  MaxStemHeight_m, data = data_train) %>%
      step_log(MaxStemHeight_m, base = 10, offset = 1),
    "recipe4" = recipe(SBMs.SB_H ~  AverageFruitLength_cm, data = data_train) %>% 
      step_log(AverageFruitLength_cm, base = 10, offset = 1),
    "recipe5" = recipe(SBMs.SB_H ~ AverageFruitLength_cm, data = data_train)
  )
  
  
  # Specify model using parsnip
  model_spec <- multinom_reg() %>%
    set_engine("nnet") %>%
    set_mode("classification")
  
  
  # Create a list of workflows using purrr::map
  workflows <- map(rec_list, ~workflow() %>%
                     add_recipe(.x) %>%
                     add_model(model_spec) %>%
                     fit(data_train))
  }


# Augment the workflows and calculate the ROC data
roc_data_palm <- map(workflows, ~.x %>%
                  augment(data_test) %>%
                  roc_curve(truth = SBMs.SB_H,.pred_1:.pred_7 ))



# Combine the ROC data into a single data frame
roc_data_palm_combined <- bind_rows(roc_data_palm, .id = "Model")

## Mammals 
## refit using nnet 
library(nnet)
refit_mammal <- nnet::multinom(SBMs.SB_W ~ log(BodyMass.Value) + Diet.Fruit, data = MammNet)
colnames(MammNet)

# Assuming your data frame is named df
# Specify the outcome variable
outcome_var <- "SBMs.SB_W"

# Specify the columns to exclude
exclude_vars <- c("Ng.WaspNames","Scientific", "MSWFamilyLatin", outcome_var)

# Create the formula dynamically
predictor_vars <- setdiff(names(MammNet), exclude_vars)
formula <- as.formula(paste(outcome_var, "~", paste(predictor_vars, collapse = " + ")))

## Mammals 
# Fit the multinomial logistic regression model
refit_mammal <- multinom(formula, data = MammNet)
Pred_trait_data <- data.frame(MammNet, "pred" = predict(refit_mammal))

## palms
## refit using nnet 
refit_palm <- nnet::multinom(SBMs.SB_H ~  MaxStemHeight_m + AverageFruitLength_cm + Acaulescent + Erect,
 data = PalmNet[complete.cases(PalmNet),])
Pred_trait_data_palm <- data.frame(PalmNet[complete.cases(PalmNet),], "pred" = predict(refit_palm))

var_im_palm <- t((caret::varImp(refit_palm)))
colnames(var_im_palm) <- str_remove(colnames(var_im_palm), 'PalmTribe')


PalmPreds <- data.frame("spNamePalm" = palm_traits$SpecName,
                        "group" = predict(refit_palm, palm_traits, allow.new.levels = T))


# predict mammal

mammPreds <- data.frame("spNameMam" = MammNet$Scientific,
                        "group" = predict(refit_mammal,MammNet))


# get the assemblages 

palm_grids <- readRDS("00_Data/00_species_distribution/gridded_palm_data.RDS")
mammal_grids <- readRDS("00_Data/00_species_distribution/gridded_mammal_data.RDS")


palm_grids <- palm_grids %>% set_names(str_replace(str_remove(basename(palm_shp_files), '.shp'),'_', " "))
palm_grids <- keep(palm_grids,~ !is.null(.x$result))

# Get the centroids of gridded ranges
sf::sf_use_s2(FALSE)

centroids_mammals <- mammal_grids %>% imap(~st_centroid(.x) %>% 
                                             st_coordinates() %>%
                                             data.frame() %>% 
                                             mutate(id = .y, 
                                                    area = st_area(.x))) %>% 
  bind_rows()


get_palm_centroids <-  function(palm_grids){
  palm_grids %>% imap(~st_centroid(.x$result) %>% 
                        st_coordinates() %>%
                        data.frame() %>% 
                        mutate(id = .y, 
                               area = st_area(.x$result), 
                               X1 = NULL,
                               x2 = NULL)) %>% 
    bind_rows() 
}

## Make a safe version of the ftion to avoid errors
safe_get_palm_centroids <- safely(get_palm_centroids)

centroids_palms <- safe_get_palm_centroids(palm_grids)

#make assemblages for all species 

all_assemblages <- centroids_mammals %>% 
  rbind(centroids_palms$result %>% dplyr::select(!X2))


# round to 2 decimals 
all_assemblages <- all_assemblages %>% 
  mutate(taxa = case_when(id %in% palm_traits$SpecName ~ 'palm', 
                          id %in% mammal_traits$Scientific ~ 'mammals',
                          TRUE~NA_character_), 
         grid_id = paste0(X,'_', Y)) 

# transform centroids to features
all_assemblages <- st_as_sf(all_assemblages, coords = c('X', 'Y'), crs = st_crs(grid))
# set right crs
all_assemblages <- st_set_crs(all_assemblages,value = st_crs(grid) )
# intersect back with grid
int <- st_intersects(all_assemblages$geometry, grid)
# add grid id
all_assemblages$grid <- unlist(int)

all_preds_sbm <- rbind(PalmPreds %>% setNames(c('id', 'SBM_G')), mammPreds %>% setNames(c('id', 'SBM_G')))

mammPreds$group  |> unique()

# join trait data
all_assemblages <- all_assemblages %>% 
  left_join(all_preds_sbm, c('id'))


# count species numbers
table_taxa_grid <- all_assemblages %>% 
  split(.$grid) %>% 
  imap(~{
    (table(.$taxa)) %>% data.frame() %>% mutate(id = .y)
  })


table_taxa_grid <- table_taxa_grid %>%
  bind_rows()

richtab <- xtabs(Freq~id+Var1, table_taxa_grid) 
richtab <- richtab[(richtab[,1]>5 & richtab[,2] > 5),] %>% rownames()

### filter those grids with at least 5 species

all_assemblages_prunned <- all_assemblages %>% 
  filter(grid %in% richtab)


## Read in final product (gridded assemblages)
all_assemblages_prunned <- readRDS( '00_Data/02_species_interactions/Metaweb.RDS')

