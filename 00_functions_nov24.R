
# Defining a function to compute Youden's J statistic

#' @title TestYJ 
#' @description Function to compute Youden's J statistic
#' @param probNet A matrix with the probabilities of interaction between species
#' @param obs A binary vector with the observed interactions
#' @param n The number of iterations to compute the statistic
#' @return A data frame with the sensitivity, specificity, and Youden's J statistic for each iteration

TestYJ <- function(probNet,obs, n){
  
  ## parameters to debug 
  
   probNet <- as.matrix(probNet)
   obs <- as.matrix(obs)
   obs[obs>0] <- 1
  # n <- 2
  
  
  
  
  sq <- seq(range(probNet)[1],
            range(probNet)[2], 
            diff(range(probNet))/n)
  sens <- c()
  speci <- c()
  YJ <- c()
  
  for (i in 1:n) {
    prob10 <- ifelse(probNet > sq[i], 1,0)

    Ttab <- prop.table(table(obs,prob10))

    sens[i] <- Ttab[4]/c(Ttab[4] + Ttab[2])
    speci[i] <- Ttab[1]/c(Ttab[1] + Ttab[3])
    
    YJ[i] <- sens[i] + speci[i] - 1
    
  }
  
  ret <- data.frame(sens, speci, YJ)
  ret
  return(ret)
  
}


# A function to calculate metrics for a grid
calc_net_metric <- function(grid_test, SBMs){
  
  
  #grid_test <- all_assemblages_prunned %>% split(.$grid) %>% pluck(1)
  
  a <- grid_test %>% 
    distinct(id, taxa, SBM_G) %>% 
    split(.$taxa)
  
  ### NOte here the switch between mammals and palms started :P 
  
  fr_palm <- table(a[[1]]$SBM_G)
  
  fr_mammals <- table(a[[2]]$SBM_G)
  
  ## compute normalized asymmetry 
  
  fr_norm_palm <- fr_palm/sum(fr_palm)
  
  fr_norm_mammals <- fr_mammals/sum(fr_mammals)
  
  fta <- abs(fr_norm_palm - fr_norm_mammals)
  
  ## compute specialization
  n <- expand.grid(pluck(a,'mammals', 'id'), pluck(a,'palm', 'id')) 
  
  
  areas <- grid_test %>% 
    split(.$taxa) %>% 
    map(~{
      .x %>% 
        mutate(area = as.numeric(area)) %>% 
        group_by(id) %>% 
        summarize(area_sum = sum(area))
      
    }) %>% 
    bind_rows() 
  
  n <- n %>% 
    left_join(pluck(a,'mammals'), by = c('Var1' = 'id')) %>% 
    left_join(pluck(a,'palm'), by = c('Var2' = 'id')) %>% 
    left_join(areas, by = c('Var1' = 'id') ) %>% 
    left_join(areas, by = c('Var2' = 'id') )
  
  n$intPro <- sapply(1:length(n$Var1), function(i) 
    (SBMs$SBM1$Omega_rs[n$SBM_G.x[i], n$SBM_G.y[i]]))
  
  n <- n %>% 
    mutate(int_area = ((area_sum.x ) / sum(area_sum.x,area_sum.y)) * ((area_sum.y ) / sum(area_sum.x,area_sum.y)))
  
  n$n_geog_dist <- st_distance(x = n$geometry.x,y = n$geometry.y) %>% 
    diag()
  
  rescaled_dist <- scales::rescale(as.numeric(n$n_geog_dist),c(0,1))
  
  n$int_final <- scales::rescale(n$intPro,c(0,1)) * scales::rescale(n$int_area,c(0,1)) * (1 -rescaled_dist)
  
  netT <- xtabs(int_final~Var1 + Var2, n)
  
  h2 <-  cassandRa::RarefyNetwork(netT,
                                  abs_sample_levels = 100,
                                  metrics = c("H2"))
  
  
  h2 <- h2$H2 |> median()
  
  
  return(list(fr_palm = fr_palm, 
              fr_mammals = fr_mammals, 
              fr_norm_palm = fr_norm_palm, 
              fr_norm_mammals = fr_norm_mammals, 
              fta = fta, 
              netT = netT, 
              h2 = h2))
}    


cal_net_metric_safe <- safely(calc_net_metric)

# A function to calculate metrics for a grid (for zcores)
calc_net_metric2 <- function(grid_test, SBMs){
  
  
  #grid_test <- expected_comm
  
  a <- grid_test %>% 
    distinct(id, taxa, SBM_G) %>% 
    split(.$taxa)
  
  
  fr_palm <- table(a[[1]]$SBM_G)
  
  fr_mammals <- table(a[[2]]$SBM_G)
  
  ## compute normalized asymmetry 
  
  fr_norm_palm <- fr_palm/sum(fr_palm)
  
  fr_norm_mammals <- fr_mammals/sum(fr_mammals)
  
  fta <- abs(fr_norm_palm - fr_norm_mammals)
  
  ## compute specialization
  n <- expand.grid(pluck(a,'mammals', 'id'), pluck(a,'palm', 'id')) 
  
  
  
  n <- n %>% 
    left_join(pluck(a,'mammals'), by = c('Var1' = 'id')) %>% 
    left_join(pluck(a,'palm'), by = c('Var2' = 'id')) 
  
  n$intPro <- sapply(1:length(n$Var1), function(i) 
    (SBMs$SBM1$Omega_rs[n$SBM_G.x[i], n$SBM_G.y[i]]))
  
  netT <- xtabs(intPro~Var1 + Var2, n)
  
  h2 <-  cassandRa::RarefyNetwork(netT,
                                  abs_sample_levels = 100,
                                  metrics = c("H2"))
  
  
  h2 <- h2$H2 |> median()
  
  
  return(list(fr_palm = fr_palm, 
              fr_mammals = fr_mammals, 
              fr_norm_palm = fr_norm_palm, 
              fr_norm_mammals = fr_norm_mammals, 
              fta = fta, 
              netT = netT, 
              h2 = h2))
}    




# Define function to compute expected values 

get_expected_val <- function(all_assemblages_prunned_biog, grid_to_sample, SBMs){
  
  biog_to_sample <- (all_assemblages_prunned_biog$Dominions[all_assemblages_prunned_biog$grid == grid_to_sample ] |> 
                       table() |> sort(decreasing = T))[1] |> names()
  
  expected_comm <- 
    all_assemblages_prunned_biog |>
    filter(Dominions %in% biog_to_sample) %>% 
    split(.$taxa) |>
    map(~ .x %>% 
          group_by(id) |> 
          slice(1) |> 
          ungroup() |>
          slice_sample(n = 10)) |>
    bind_rows()
  
  
  
  return(calc_net_metric2(expected_comm, SBMs))
  
  
}



## Make safe version of the function

safe_expected_values <- function(n_rep, all_assemblages_prunned_biog, grids_to_sample, SBMs){
  
  replicate(n_rep,get_expected_val(all_assemblages_prunned_biog, grids_to_sample, SBMs))
  
}

safe_expected_values <- safely(safe_expected_values)




add_clim_data <- function(z_score_table, coordinates_grid){
  # z_score_table <- fr_palms
  
  
  # z_score_table <- z_scores_across_biogeog
  grid_coords <- as.data.frame(coordinates_grid[as.numeric(z_score_table$grid),])
  z_score_table <-  cbind(z_score_table,grid_coords)
  
  ## make xy dataframe as simple feature
  grid_coords <- st_as_sf(grid_coords, coords = c("X", "Y"), crs = st_crs(neotropics))
  
  
  ## add id to the models 
  
  clim_var <- 
    data.frame(
      'Temp' =  raster::extract(Temp, sf::st_as_sf(grid_coords)),
      'Prec' =  raster::extract(Prec,sf::st_as_sf(grid_coords)),
      'TS' =  raster::extract(TempSeaso,sf::st_as_sf(grid_coords)),
      'PS' =  raster::extract(PrecSe,sf::st_as_sf(grid_coords))
    )
  
  ggsf <- data.frame(z_score_table,clim_var)
  ggsf <- na.omit(ggsf)
  
  ggsf <- ggsf |>
    reshape2::melt(id.vars = c('grid', 'Temp', 'Prec', 'TS', 'PS'), value.name = c('obs_ab'), variable.name = 'SBM_G') |>
    reshape2::melt(id.vars = c('grid','SBM_G', 'obs_ab'), value.name = c('clim_val'), variable.name = 'clim_var') |> 
    filter(!SBM_G %in% c('X', 'Y'))
  
  return(ggsf)
  
}


add_clim_data2 <- function(z_score_table, coordinates_grid){
  # z_score_table <- fr_palms
  
  
  # z_score_table <- z_scores_across_biogeog
  grid_coords <- as.data.frame(coordinates_grid[as.numeric(z_score_table$grid),])
  z_score_table <-  cbind(z_score_table,grid_coords)
  
  ## make xy dataframe as simple feature
  grid_coords <- st_as_sf(grid_coords, coords = c("X", "Y"), crs = st_crs(neotropics))
  
  
  ## add id to the models 
  
  clim_var <- 
    data.frame(
      'Temp' =  raster::extract(Temp, sf::st_as_sf(grid_coords)),
      'Prec' =  raster::extract(Prec,sf::st_as_sf(grid_coords)),
      'TS' =  raster::extract(TempSeaso,sf::st_as_sf(grid_coords)),
      'PS' =  raster::extract(PrecSe,sf::st_as_sf(grid_coords))
    )
  
  ggsf <- data.frame(z_score_table,clim_var)
  ggsf <- na.omit(ggsf)
  
  
  return(ggsf)
  
}