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
library(car)
source("scripts/functions.R") 

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
               degree <- degree(igraph)
               igraph <- set_vertex_attr(igraph, "strength", value = strength)
               igraph <- set_vertex_attr(igraph, "degree", value = degree)
return(igraph)
}

## Using func_igraph
igraph_objects <- lapply(rep_list, func_igraph)
func_plot_network(igraph_objects[[1]], node_size = 10)

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
attr_strength$treatment <- relevel(attr_strength$treatment, "two")

############### SIGNIFICANCE TESTING MALE VS. FEMALE STRENGTH ######
## Linear model with strength as response variable
strength_glmm <- lmer(data = attr_strength, log(strength) ~ sex*treatment + (1|block))
summary(strength_glmm)
Anova(strength_glmm)
plot(simulateResiduals(strength_glmm))

## Linear model with degree as response variable
degree_glmm <- lmer(data = attr_strength, degree ~ sex*treatment + (1|block))
summary(degree_glmm)
Anova(degree_glmm)
plot(simulateResiduals(degree_glmm)) # NOT GOOD

e_strength_2 <- emmeans(degree_glmm, c("sex", "treatment"))
pairs(e_strength_2)
plot(e_strength_2)

## Extracting model coefficients

# Model coefficient for main effect of sex
main_effect_observed <- summary(strength_glmm)$coefficients[2,3]

# Contrasts model
e_strength <- emmeans(strength_glmm, c("sex", "treatment"))
pairs(e_strength)
plot(e_strength)

# Males vs. females for two shelter treatment
two_t_ratio <- as.data.frame(pairs(e_strength))[1,5]  

# Males vs. females for twelve shelter treatment
twelve_t_ratio <- as.data.frame(pairs(e_strength))[6,5]

#################### STRENGTH VISUALIZATION #############################
## Boxplot
levels(attr_strength$treatment) <- c("two shelter", "twelve shelter")

ggplot(data = attr_strength, aes(y = strength, x = treatment, fill = sex, color = sex)) + 
      geom_boxplot(alpha = 0.6, outlier.shape = NA) + 
      scale_fill_manual(values = c("sandybrown", "skyblue3")) + 
      scale_color_manual(values = c("sandybrown", "skyblue3")) +
      geom_point(position = position_jitterdodge(), alpha = 0.8, size = 0.5)
#      geom_jitter(alpha = 0.5, width = 0.2)

#################### DEGREE VISUALIZATION #############################
## Boxplot
levels(attr_strength$treatment) <- c("two shelter", "twelve shelter")

ggplot(data = attr_strength, aes(y = degree, x = treatment, fill = sex, color = sex)) + 
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + 
  scale_fill_manual(values = c("sandybrown", "skyblue3")) + 
  scale_color_manual(values = c("sandybrown", "skyblue3")) +
  geom_point(position = position_jitterdodge(), alpha = 0.8, size = 0.5)
#      geom_jitter(alpha = 0.5, width = 0.2)






######### PERMUTATION TEST FOR MAIN EFFECT OF SEX #############
n_sim <- 999
set.seed(33)
sim_coefs_1 <- numeric(n_sim)

for (i in 1:n_sim){
  # Creates new igraph objects where the nodes are shuffled
  random_igraphs <- lapply(rep_list, func_permute_igraph)
  # Runs the glm on the new shuffled igraph objects; save coefs
  sim_coefs_1[i] <- func_random_model_p1(random_igraphs, statistic = "main") 
}

# Plot histogram 
sim_coefs_1 <- c(sim_coefs_1, main_effect_observed)
hist(sim_coefs_1, main = "Main effect of sex", xlab = "t-value value for sexMale", col = "azure2", breaks = 40)
lines(x = c(main_effect_observed, main_effect_observed), y = c(0, 270), col = "red", lty = "dashed", lwd = 2) 

# Obtain p-value
if (main_effect_observed >= mean(sim_coefs_1)) {
  pred1_p <- 2*mean(sim_coefs_1 >= main_effect_observed) } else {
    pred1_p <- 2*mean(sim_coefs_1 <= main_effect_observed)
  }

# Add p-value to histogram
#text(x = 0.17, y = 40, "p = 0.006")


######### PERMUTATION TEST FOR TWO SHELTER CONTRAST ############
set.seed(33)
contrast_two_sims <- numeric(n_sim)

for (i in 1:n_sim){
  # Creates new igraph objects where the nodes are shuffled
  random_igraphs <- lapply(rep_list, func_permute_igraph)
  # Runs the glm on the new shuffled igraph objects; save coefs
  contrast_two_sims[i] <- func_random_model_p1(random_igraphs, statistic = "contrast_two") 
}

# Plot histogram 
contrast_two_scores <- c(contrast_two_sims, two_t_ratio)
hist(contrast_two_scores, main = "Two shelter contrast", xlab = "t-value value for sexMale", col = "azure2", breaks = 30)
lines(x = c(two_t_ratio, two_t_ratio), y = c(0, 270), col = "red", lty = "dashed", lwd = 2) 

# Obtain p-value
if (two_t_ratio >= mean(contrast_two_scores)) {
  con_2_p <- 2*mean(contrast_two_scores >= two_t_ratio) } else {
    con_2_p <- 2*mean(contrast_two_scores <= two_t_ratio)
  }

# Add p-value to histogram
#text(x = 0.17, y = 40, "p = 0.006")


######### PERMUTATION TEST FOR TWELVE SHELTER CONTRAST ############
set.seed(33)
contrast_twelve_sims <- numeric(n_sim)

for (i in 1:n_sim){
  # Creates new igraph objects where the nodes are shuffled
  random_igraphs <- lapply(rep_list, func_permute_igraph)
  # Runs the glm on the new shuffled igraph objects; save coefs
  contrast_twelve_sims[i] <- func_random_model_p1(random_igraphs, statistic = "contrast_twelve") 
}

# Plot histogram 
contrast_twelve_scores <- c(contrast_twelve_sims, twelve_t_ratio)
hist(contrast_twelve_scores, main = "Twelve shelter contrast", xlab = "t-value value for sexMale", col = "azure2", breaks = 20)
lines(x = c(twelve_t_ratio, twelve_t_ratio), y = c(0, 270), col = "red", lty = "dashed", lwd = 2) 

# Obtain p-value
if (twelve_t_ratio >= mean(contrast_twelve_scores)) {
  con_12_p <- 2*mean(contrast_twelve_scores >= twelve_t_ratio) } else {
    con_12_p <- 2*mean(contrast_twelve_scores <= twelve_t_ratio)
  }

# Add p-value to histogram
#text(x = 0.17, y = 40, "p = 0.126")








