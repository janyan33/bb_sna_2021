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
source("scripts/functions.R") ## ADD LATER
source("scripts/igraphplot2.R") ## ADD LATER


###################### INPUTTING AND ORGANIZING DATA ########################
## Data for aggregation-based networks
groups <- read.csv("data/bbsna_aggregations.csv") %>% 
          filter(Network != "N/A")

rep_list <- split(groups, groups$Network) # creates a list of replicates

## Import and clean raw attribute data
attr <- read.csv("data/bbsna_attributes_raw.csv") %>% 
        filter(network != "prelim") %>% 
        filter(exclude != "yes")

################# CREATING AGGREGATION NETWORKS ######################
## Function that creates a list of igraph objects (1 per replicate)
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

## Using func_igraph
igraph_objects <- lapply(rep_list, func_igraph)

################ CALCULATING STRENGTH VALUES FOR EACH INDIVIDUAL ####
func_attr <- function(igraph_objects){
             new_attr <- data.frame()
             for (i in 1:length(igraph_objects)){
                        attr_i <- subset(attr, network == i & exclude != "yes")
                        igraph_objects[[i]] <- set_vertex_attr(igraph_objects[[i]], "size", value = attr_i$size)
                        igraph_objects[[i]] <- set_vertex_attr(igraph_objects[[i]], "treatment", value = attr_i$treatment)
                        igraph_objects[[i]] <- set_vertex_attr(igraph_objects[[i]], "network", value = attr_i$network)
                        igraph_objects[[i]] <- set_vertex_attr(igraph_objects[[i]], "block", value = attr_i$block)
                        (new_attr <- rbind(new_attr, vertex_attr(igraph_objects[[i]])))
             }
return(new_attr)
}

attr_strength <- func_attr(igraph_objects)

attr_strength$network <- as.factor(attr_strength$network)
attr_strength$block <- as.factor(attr_strength$block)
attr_strength$treatment <- as.factor(attr_strength$treatment)
attr_strength$sex <- as.factor(attr_strength$sex)
attr_strength$size <- as.numeric(attr_strength$size)

############### SIGNIFICANCE TESTING MALE VS. FEMALE STRENGTH ######
strength_glmm <- glmer(data = attr_strength, strength ~ sex + treatment + (1|size) + (1|block), family = Gamma(link = "log"))
summary(strength_glmm)

# Permutation function (creates new shuffled igraphs, one per rep)
func_permute_igraph <- function(rep_list) { 
                       group_list <- strsplit(rep_list$Members, " ")
                       gbi_matrix <- get_group_by_individual(group_list, data_format = "groups")
                       ibi_matrix <- get_network(gbi_matrix, data_format = "GBI")
                       
                       #shuffle names 
                       new_names <- sample(colnames(ibi_matrix))
                       colnames(ibi_matrix) <- new_names
                       rownames(ibi_matrix) <- new_names
  
                       igraph <- graph_from_adjacency_matrix(ibi_matrix, diag = FALSE, weighted = TRUE, mode = "undirected")
                       igraph <- set_vertex_attr(igraph, "sex", 
                       value = ifelse(V(igraph)$name %in% LETTERS[1:12], "Male", "Female"))
                       strength <- strength(igraph)
                       igraph <- set_vertex_attr(igraph, "strength", value = strength)
                       return(igraph)
}  

#
n_sim <- 999
set.seed(33)
sim_coefs <- numeric(n_sim)

for (i in 1:n_sim_1){
  random_igraphs <- lapply(rep_list_groups, func_permute_igraph)
  # Runs the glm on the new shuffled igraph objects; save coefs
  sim_coefs_1[i] <- func_random_model_p1(random_igraphs) 
}
# Plot histogram 
sim_coefs_1 <- c(sim_coefs_1, obs_coefs_mat[2])
hist(sim_coefs_1, main = "Prediction 1", xlab = "Coefficient value for sexMale", col = "slategray2", breaks = 40)
lines(x = c(obs_coefs_mat[2], obs_coefs_mat[2]), y = c(0, 270), col = "red", lty = "dashed", lwd = 2) 

# Obtain p-value
if (obs_coefs_mat[2] >= mean(sim_coefs_1)) {
  pred1_p <- 2*mean(sim_coefs_1 >= obs_coefs_mat[2]) } else {
    pred1_p <- 2*mean(sim_coefs_1 <= obs_coefs_mat[2])
  }



