##setwd("C:/Users/janya/Desktop/R/bedbugs/bb_sna_2021")

library(tidyverse)
library(asnipe)
library(igraph)
library(ggplot2); theme_set(theme_classic())
library(lme4)
library(glmmTMB)
library(assortnet)
library(intergraph)
library(DHARMa)
library(emmeans)
library(janitor)
source("scripts/igraphplot2.R") ## ADD LATER

## Importing and organizing data
mate_mount_data <- read.csv("data/mate_mount_data.csv") %>% 
                   select(-patch_partner)

mate_data <- mate_mount_data %>% 
             filter(behaviour == "mating") 

mate_data_reps <- split(mate_data, mate_data$network)

mount_data <- mate_mount_data %>% 
              filter(behaviour == "mount")

mount_data_reps <- split(mount_data, mount_data$network)

mount_matrices <- readRDS("mount_matrices.rds")
mating_matrices <- readRDS("mating_matrices.rds")

## Number of matings per individual per rep
table(mate_data_reps[[1]]$focal_individual)
table(mate_data_reps[[2]]$focal_individual)
table(mate_data_reps[[3]]$focal_individual)
table(mate_data_reps[[4]]$focal_individual)
table(mate_data_reps[[5]]$focal_individual)
table(mate_data_reps[[6]]$focal_individual)

table(mate_data_reps[[1]]$social_partner)
table(mate_data_reps[[2]]$social_partner)
table(mate_data_reps[[3]]$social_partner)
table(mate_data_reps[[4]]$social_partner)
table(mate_data_reps[[5]]$social_partner)
table(mate_data_reps[[6]]$social_partner)

## Number of mounts performed
table(mount_data_reps[[1]]$focal_individual)
table(mount_data_reps[[2]]$focal_individual)
table(mount_data_reps[[3]]$focal_individual)
table(mount_data_reps[[4]]$focal_individual)
table(mount_data_reps[[5]]$focal_individual)
table(mount_data_reps[[6]]$focal_individual)

## Number of mounts received 
table(mount_data_reps[[1]]$social_partner)
table(mount_data_reps[[2]]$social_partner)
table(mount_data_reps[[3]]$social_partner)
table(mount_data_reps[[4]]$social_partner)
table(mount_data_reps[[5]]$social_partner)
table(mount_data_reps[[6]]$social_partner)

## Function that creates mating/mounting networks
func_matrix_to_igraph <- function(matrix, mode, behaviour){
                         igraph <- graph_from_adjacency_matrix(matrix, diag = FALSE, weighted = TRUE, mode = mode)
                         igraph <- set_vertex_attr(igraph, "sex", 
                                   value = ifelse(V(igraph)$name %in% LETTERS[1:12], "Male", "Female"))
                         strength <- strength(igraph, mode = "in")
                         igraph <- set_vertex_attr(igraph, behaviour, value = strength)
                         V(igraph)$color <- ifelse(V(igraph)$sex == "Female", "sandybrown", "skyblue3")
                         V(igraph)$label.color <- "white"
                         V(igraph)$size <- V(igraph)$behaviour*3.5
                         E(igraph)$width <- E(igraph)$weight*1.5
                         plot(igraph, edge.color = "dimgrey", layout = layout_nicely(igraph))
return(igraph)
}

func_matrix_to_igraph(mating_matrices[[1]], mode = "undirected", behaviour = "mating")









