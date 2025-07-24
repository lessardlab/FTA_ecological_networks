library(furrr)
library(tidyverse)
source("C:/Users/gabri/Documents/PhD/00_Chapter_palms_mammal_interactions/R-analisis/Pruned_scripts/00_functions_nov24.R")

#==========================================================================================================
# Calculate functional diversity and network structure metrics
#============================================================================================================
all_assemblages_prunned2 <- all_assemblages_prunned
all_assemblages_prunned2$geometry <- NULL
all_assemblages_prunned2 |> group_by(SBM_G)



plan(multisession, workers = 4)

my_net__output <- 
  all_assemblages_prunned2[1:1000,] %>% split(.$grid) %>% 
  furrr::future_map(function(grid) {
    cal_net_metric_safe(grid, SBMs)
  }) 




my_networks <- readRDS('00_Data/02_species_interactions/final-networks-grid_prevised.RDS')
my_networks <- keep(my_networks, ~is.null(.x$error))

## Extract results 

fta_obs <-  my_networks |>  map(~.x$result) |>
  map(~.x$fta |> unlist()) |>
  bind_rows() |>  mutate('grid' = names(my_networks))

fr_palms <-  my_networks |>   map(~.x$result) |>
  map(~.x$fr_palm |> unlist()) |>
  bind_rows() |>  mutate('grid' = names(my_networks))

fr_mammals <-   my_networks |>  map(~.x$result) |>
  map(~.x$fr_mammals |> unlist()) |>  bind_rows() |>
  mutate('grid' = names(my_networks))

fr_norm_palms <-  my_networks |>  map(~.x$result) |>
  map(~.x$fr_norm_palm |> unlist()) |>
  bind_rows() |> mutate('grid' = names(my_networks))

fr_norm_mammals <- my_networks |> 
  map(~.x$result) |>
  map(~.x$fr_norm_mammals |> unlist()) |>
  bind_rows() |>
  mutate('grid' = names(my_networks))

h2_grid <-  my_networks |> 
  map(~.x$result) |>
  map(~.x$h2 |> unlist()) |>
  unlist() |>
  data.frame() |> 
  setNames('h2') |>
  rownames_to_column('grid')

h2_grid$h2 |> range(na.rm = T)
#=================
## NOT PART OF EC MONOGRAPH PAPER 
#=======================================================================================================
# Run null model replicates 
#=======================================================================================================

## add biogeographic region based on xy 

xy_sf <- st_as_sf(all_assemblages_prunned, coords = c("cord_x", "cord_y"), crs = st_crs(neotropics))

## make sure crs match 

xy_sf <- st_set_crs(xy_sf, st_crs(neotropics))

all_assemblages_prunned_biog <- st_join(xy_sf, neotropics)
all_assemblages_prunned_biog <- all_assemblages_prunned_biog |> 
  group_by(id, taxa, grid, Dominions) |>  slice(1)
all_assemblages_prunned_biog2 <-all_assemblages_prunned_biog
all_assemblages_prunned_biog2$geometry <- NULL

grids_to_sample <- all_assemblages_prunned$grid |> unique() 
cl <- NULL
## open a parallel cluster to run safe_expected_values in parallel 

library(parallel)
library(foreach)
library(doParallel)

# Register parallel backend
cl <- makeCluster(10)
registerDoParallel(cl)
# Export variables and libraries to each cluster
clusterExport(cl, c('all_assemblages_prunned_biog', 'grids_to_sample', 'SBMs', 'calc_net_metric2', 'get_expected_val','safe_expected_values'))
clusterEvalQ(cl, {
  library(tidyverse)
  library(sf)
  library(cassandRa)
  library(vegan)
})
# parallel::stopCluster(cl)
gsample <- (grids_to_sample)
nrep <- 999
## uncomment to activate 

### my_null_result_full <- foreach(grid = gsample, 
###        .packages = c('tidyverse', 'sf', 'cassandRa', 'vegan')) %dopar% {
###         safe_expected_values(nrep, all_assemblages_prunned_biog, grid, SBMs)
  
###                 }

my_null_result_full <- readRDS('00_Data/02_species_interactions/null-networks-grid_final_all2.RDS')
names(my_null_result_full) <- gsample

#=======================================================================================================
# Preparing data to model 
#=======================================================================================================
mres <- keep(my_null_result_full, ~ !is.null(.x$result)) 
mres <- mres |> map(~.x$result)
mres_zscore <- mres 

fta_expected_mean <-   mres_zscore |> 
  map(~.x['fta',] |> bind_rows() |> colMeans()) |>
  bind_rows() |>  mutate(grid = names(mres_zscore))

fta_expected_sd <-   mres_zscore |> 
  map(~.x['fta',] |> bind_rows() |> apply(2,sd)) |>
  bind_rows() |>  mutate(grid = names(mres_zscore))

fr_palm_mean <-   mres_zscore |> 
  map(~.x['fr_palm',] |> bind_rows() |> colMeans()) |>
  bind_rows() |>  mutate(grid = names(mres_zscore))

fr_mammals_mean <-   mres_zscore |> 
  map(~.x['fr_mammals',] |> bind_rows() |> colMeans()) |>
  bind_rows() |>  mutate(grid = names(mres_zscore))

fr_norm_palm_mean <-   mres_zscore |> 
  map(~.x['fr_norm_palm',] |> bind_rows() |> colMeans()) |>
  bind_rows() |>  mutate(grid = names(mres_zscore))

fr_norm_mammals_mean <-  mres_zscore |> 
  map(~.x['fr_norm_mammals',] |> bind_rows() |> colMeans()) |>
  bind_rows() |>  mutate(grid = names(mres_zscore))

h2_mean <-  mres_zscore |>  map(~.x['h2',]|> unlist() |> mean(na.rm = T) ) |>
  unlist() |> data.frame() |>  setNames('h2_x')|>
  mutate(grid = names(mres_zscore))

h2_sd <- mres_zscore |> map(~.x['h2',]|> unlist() |> sd(na.rm = T) ) |>
  unlist() |>  data.frame() |>  setNames('h2sd')|>
  mutate(grid = names(mres_zscore))

h2_obs <- mres_zscore |> imap(~.x['h2',] |> unlist() |> data.frame() |> setNames('h2_obs') |>
         mutate(grid = .y) |>  mutate(rep = 1:50))|> bind_rows()

# Define the original vector
original_vector <- 1:length(mres_zscore) 

# Define the size of each smaller vector
size <- 100

# Split the original vector
split_vectors <- split(original_vector, gl(ceiling(length(original_vector) / size), size, length(original_vector)))

full_fta_expected <-
  
  split_vectors |>
  map(function(split_vec){
    split_vec |> 
      map(~{
        1:50 |>
          map(~{
            expand.grid(
              mres_zscore[[1]]['fr_norm_palm',.x] |> bind_rows() |> as.matrix(),
              mres_zscore[[1]]['fr_norm_mammals',.x] |> bind_rows() |> as.matrix()
            ) 
            
          }) |>
          bind_rows() |>
          mutate(lab = bind_rows(
            replicate(50,
                      expand.grid(matrix(rep(1:7),ncol = 7, byrow = T),
                                  matrix(rep(1:7),ncol = 7, byrow = T)),
                      simplify = FALSE)) |> 
              mutate(label = paste0('p', Var1, 'm', Var2)) |> 
              dplyr::pull(label)  ) |>
          mutate(fta = abs(Var1 - Var2)) |>
          mutate(grid = grids_to_sample[.x]) |>
          group_by(lab) |>
          mutate(h2_obs = mres_zscore[[.x]]['h2',] |>unlist()) |>
          ungroup()
      })  |>
      bind_rows()
  }) |>
  bind_rows()
#============================================================================
#================= 



fr_palms_with_climate <- add_clim_data(fr_norm_palms, coordinates_grid)

fr_norm_palm_mean_with_climate <- add_clim_data(fr_norm_palm_mean, coordinates_grid)

#=======================================================================================================
# Fit generalized additive models
#=======================================================================================================

mamal_fr_to_model <- fr_palms_with_climate |> pivot_wider(names_from = clim_var, values_from = clim_val)

palms_fr_to_model <- fr_mammals_with_climate |> pivot_wider(names_from = clim_var, values_from = clim_val)

mamal_fr_to_model <- readRDS('00_Data/mammal_fr.RDS')

palms_fr_to_model <- readRDS('00_Data/palm_fr.RDS')

#=======================================================
# Fit models relating climate to FR 
#=======================================================

library(mgcv)

model_gam_fr_mammal <- gam(obs_ab ~ (Temp) + (Prec) +(TS) + (PS) + 
                             s(Temp, by = as.factor(SBM_G)) +
                                 s(Prec, by = as.factor(SBM_G)) +
                                     s(TS, by = as.factor(SBM_G)) + 
                                        s(PS, by = as.factor(SBM_G)) ,
                 data = mamal_fr_to_model)

summary(model_gam_fr_mammal)

model_gam_fr_palm <- gam(obs_ab ~ (Temp) + (Prec) + (TS) + (PS) + 
                   s(Temp, by = as.factor(SBM_G)) + 
                   s(Prec, by = as.factor(SBM_G)) +
                   s(TS, by = as.factor(SBM_G)) + 
                   s(PS, by = as.factor(SBM_G)),
                 data = palms_fr_to_model)

#=======================================================
## Fit models relating climate to FTA 
#=======================================================

full_fta_val_wt_clim <- add_clim_data2(full_fta_val, coordinates_grid)

unique(full_fta_val_wt_clim$lab) |> unique()




which(full_fta_val_wt_clim$fta != 0) |> length() / 
(which(full_fta_val_wt_clim$fta == 0) |> length() + 
   which(full_fta_val_wt_clim$fta != 0) |> length())

full_fta_val_wt_clim$grid |> unique() |> length()



## remove infinite and na values from zscore

full_fta_val_wt_clim <-  full_fta_val_wt_clim |>  filter(!is.infinite(zscore), !is.na(zscore))

# full_fta_val_wt_clim <- readRDS('00_Data/04_models/fta_data_to_model.RDS')

# Testing the influence of interaction strenght instead of categorical guild identity 

int_strenght <- SBMs$SBM1$Omega_rs %>% reshape2::melt()
int_strenght$lab <- paste0('p',int_strenght$Var1,'m',int_strenght$Var2 )
int_strenght <- int_strenght |> rename('int_str' = 'value')
full_fta_val_wt <- full_fta_val_wt_clim |>  left_join(int_strenght, by = 'lab')


# full_fta_val_wt <- readRDS('00_Data/04_models/fta_full_with_str.RDS')


model_gam_str <- gam(fta ~
                  s(Temp) + s(Temp, int_str) +  # Interaction with continuous variable
                   s(PS) + s(PS, int_str) +  # Interaction with continuous variable
                   s(TS) + s(TS, int_str) +
                   s(Prec) + s(Prec,int_str) +
                  s(int_str),  # To model the overall influence of strength,
                 data = full_fta_val_wt)

full_fta_val_wt |> 
  summarize( fta_mean = sum(fta), 
             int_str = mean(int_str), 
             fta_w = fta_mean * int_str,
             .by = lab)|>
  arrange(fta_w)



full_fta_val_wt <- full_fta_val_wt |> left_join(h2_grid, by = 'grid', suffix = c('zscore','obs'))

full_fta_val_wt$h2obs |> range(na.rm = T)

full_fta_val_wt <- full_fta_val_wt |>  group_by(grid) |>
  summarize(sum_fta = sum(fta*int_str), 
            h2 = mean(h2obs, na.rm = T), 
            Temp = mean(Temp, na.rm = T), 
            TS = mean(TS, na.rm = T), 
            PS = mean(PS, na.rm = T), 
            Prec = mean(Prec, na.rm = T)) 

## Get the influence of asymmetry in network specialization

model_gam_h2 <- readRDS('02_Outputs/model_gam_h2.RDS')

# Fit a GAM with lab as a fixed effect
model_gam_h2 <- gam(h2 ~ (sum_fta) +
                      s(Temp) + 
                      s(PS) +
                      s(TS) + 
                      s(Prec),
                    data = full_fta_val_wt)
summary(model_gam_h2)


#========================================================
# Ecological significance test with null replicates
#=============================================================

### Observed relationship
## compute the differences for all combinations 


full_fta <-
  1:nrow(fr_norm_palms) |>
  map(function(row){
    expand.grid(
      fr_norm_palms[row,] |> dplyr::select(!grid) |> as.numeric(),
      fr_norm_mammals[row,] |> dplyr::select(!grid) |> as.numeric()) |>
      mutate(lab = expand.grid(1:7, 1:7) |> 
               mutate(label = paste0('p', Var1, 'm', Var2)) |> 
               dplyr::pull(label) ) |>
      mutate(fta = abs(Var1 - Var2)) |>
      mutate(grid = fr_norm_palms$grid[row])
    
    
    
  }) |>
  bind_rows()


full_fta_val_wt <- 
full_fta_val_wt |> 
  mutate(sum_fta_scaled = scale(sum_fta))






## Get the distribution of expected coefficients between simulated fta and h2 

full_fta_expected <- 
full_fta_expected |>
  group_by(grid,lab) |>
  mutate(rep = rep(1:50))


h2_mod_fta <- 
  full_fta_expected |>
  group_by(grid, rep) |>
  summarise(mean_fta = mean(fta, na.rm = T), 
            sd_fta = sd(fta, na.rm = T),
            h2_obs = mean(h2_obs))|>
  group_map(~{
    lm(h2_obs~mean_fta+sd_fta, data = .x) })


h2_mod_fta_coef <- 
h2_mod_fta |>
  map(~coef(.x)) |> 
  bind_rows() 

h2_mod_fta |> map(~vegan::RsquareAdj(.x)[[1]]) |> unlist()