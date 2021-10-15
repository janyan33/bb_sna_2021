## setwd("C:/Users/janya/Desktop/R/bedbugs/bb_sna_2021")

library(tidyverse)
library(asnipe)
library(igraph)
library(ggplot2); theme_set(theme_classic())
library(lme4)
library(glmmTMB)
library(assortnet)
library(janitor)

##################### INPUTTING AND ORGANIZING DATA #######################

## Data for aggregation-based networks
groups_agg <- read.csv("data/bbsna_aggregations.csv") %>%  
              filter(Network != "N/A") %>% 
              remove_empty("cols")
              
groups_agg_reps <- split(groups_agg, groups_agg$Network)


## Data for shelter-based networks
groups_shelter <- read.csv("data/bbsna_shelter_groups.csv") %>% 
                  remove_empty("cols")

groups_shelter_reps <- split(groups_shelter, groups_shelter$Network)


#################### VISUALIZING SOCIAL NETWORKS ############################

## Function for turning group data into igraph objects
func_igraph <- function(rep_groups){
               group_list <- strsplit(rep_groups$Members, " ")
               gbi_matrix <- get_group_by_individual(group_list, data_format = "groups")
               ibi_matrix <- get_network(gbi_matrix, data_format = "GBI")
               ibi_matrix <- ibi_matrix[order(rownames(ibi_matrix)) , order(colnames(ibi_matrix))] # alphabetical order
               igraph <- graph_from_adjacency_matrix(ibi_matrix, diag = FALSE, weighted = TRUE, mode = "undirected")
               igraph <- set_vertex_attr(igraph, "sex", 
                         value = ifelse(V(igraph)$name %in% LETTERS[1:12], "Male", "Female"))
               strength <- strength(igraph)
               igraph <- set_vertex_attr(igraph, "strength", value = strength)
return(igraph)
}

## Function for visualizing networks
func_plot_network <- function(igraph_object, node_size){
                     V(igraph_object)$color <- ifelse(V(igraph_object)$sex == "Female", "sandybrown", "skyblue3")
                     V(igraph_object)$size <- V(igraph_object)$strength*node_size
                     V(igraph_object)$label.color <- "white"
                     E(igraph_object)$width <- E(igraph_object)$weight*6
                     plot(igraph_object, edge.color = "dimgrey") #, vertex.label = NA) allows us to turn on/off vertex labels
}

## Visualizing aggregation-based networks
igraph_objects_agg <- lapply(groups_agg_reps, func_igraph)
lapply(igraph_objects_agg, node_size = 15, func_plot_network)

## Visualizing shelter-based networks
igraph_objects_shelter <- lapply(groups_shelter_reps, func_igraph)
lapply(igraph_objects_shelter, node_size = 12, func_plot_network)


######################## TESTING ASSORTATIVITY ############################

## Function for turning group data into ibi_matrices (data needed in this format to use the assortment.discrete function)
func_ibi <- function(rep_groups){
            group_list <- strsplit(rep_groups$Members, " ")
            gbi_matrix <- get_group_by_individual(group_list, data_format = "groups")
            ibi_matrix <- get_network(gbi_matrix, data_format = "GBI")
            ibi_matrix <- ibi_matrix[order(rownames(ibi_matrix)) , order(colnames(ibi_matrix))] # alphabetical order
return(ibi_matrix)
}

## Function that runs the assortativity permutation test
func_permute_assort <- function(ibi_matrix){
                       sex_table <- as.data.frame(colnames(ibi_matrix)) %>% 
                       rename("ID" = "colnames(ibi_matrix)") %>% 
                       mutate(sex = ifelse(ID %in% LETTERS[1:12], "Male", "Female"))
                       # Getting observed assortativity score
                       obs_assort_index <- assortment.discrete(ibi_matrix, types = sex_table$sex, weighted = TRUE)$r    
  
                       # Setting up the permutation
                       n_sim <- 999; set.seed(33)
                       sim_assort_index <- numeric(n_sim)
                       
                       # Loop that creates a shuffled ibi matrix and calculates new assort index for each iteration         
                       for (i in 1:n_sim){
                           new_names <- sample(colnames(ibi_matrix))
                           colnames(ibi_matrix) <- new_names
                           rownames(ibi_matrix) <- new_names
    
                           sex_table_new <- as.data.frame(colnames(ibi_matrix)) %>% 
                           rename("ID" = "colnames(ibi_matrix)") %>% 
                           mutate(sex = ifelse(ID %in% LETTERS[1:12], "Male", "Female"))
    
                           sim_assort_index[i] <- assortment.discrete(ibi_matrix, types = sex_table_new$sex, 
                                                  weighted = TRUE)$r  
                       }
                      sim_assort_index <- c(sim_assort_index, obs_assort_index)
                      
                      # Results in histogram
                      hist(sim_assort_index, breaks = 29, xlim = c(min(sim_assort_index), max(sim_assort_index)), 
                                             ylim = c(0, 180), col = "aliceblue")
                      lines(x = c(obs_assort_index, obs_assort_index), y = c(0, 150), col = "red", lty = "dashed", lwd = 2)
  
                      # Computing the p-value
                      if (obs_assort_index >= mean(sim_assort_index)) { 
                          p <- 2*mean(sim_assort_index >= obs_assort_index) } else {
                          p <- 2*mean(sim_assort_index <= obs_assort_index)
                      }
                      
                      list <- list("p-value" = p, "observed assortativity score" = obs_assort_index)
return(list)
}

## Assortativity for aggregation-based networks
ibi_objects_agg <- lapply(groups_agg_reps, func_ibi)
lapply(ibi_objects_agg, func_permute_assort)

## Assortativity for shelter-based networks
ibi_objects_shelter <- lapply(groups_shelter_reps, func_ibi)
lapply(ibi_objects_shelter, func_permute_assort)











