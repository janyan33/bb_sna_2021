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
library(car)
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


################# ANALYZING AND VISUALING MATING AND MOUNTING ##################
attr <- read.csv("data/bbsna_attributes_full.csv", stringsAsFactors = TRUE)
attr$network <- as.factor(attr$network)
attr$treatment <- relevel(attr$treatment, "two")

attr_fem <- attr %>%  
            filter(sex == "Female")

attr_fem$block <- as.factor(attr_fem$block)

############################## MOUNTING ANALYSES ############################### 
## Mounting rate figure
ggplot(data = attr_fem, aes(y = mounts_in_rate, x = treatment)) + geom_boxplot() + theme(text = element_text(size = 20)) + 
       geom_jitter(position = position_jitter(width = 0.15, height = 0), color = "sandybrown", alpha = 0.5, size = 1.5) + 
       labs(y = "", x = "") + scale_y_continuous(breaks = c(5, 10, 15))

## Mounting rate model
mount_model <- lmer(data = attr_fem, mounts_in_rate ~ treatment + (1|block))
summary(mount_model)
Anova(mount_model, type = 3)

# Diagnostic plots
plot(mount_model)
residuals_mount_model <- simulateResiduals(mount_model)
plot(residuals_mount_model)


############################## MATING ANALYSES ################################
## Mating rate figure
ggplot(data = attr_fem, aes(y = mating_rate, x = treatment)) + geom_boxplot() + theme(text = element_text(size = 20)) + 
       geom_jitter(position = position_jitter(width = 0.15, height = 0), color = "sandybrown", alpha = 0.5, size = 1.5) + 
       labs(y = "", x = "") + scale_y_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3))

## Mounting rate model
mating_model <- lmer(data = attr_fem, mating_rate ~ treatment + (1|block))
summary(mating_model)
Anova(mating_model, type = 3)

# Diagnostic plots
plot(mating_model)
residuals_mating_model <- simulateResiduals(mating_model)
plot(residuals_mating_model)

############################ MATING and STRENGTH ##############################
## Mating
ggplot(data = attr_fem, aes(y = mating_rate, x = strength, color = block)) + 
       geom_point() + facet_grid(~treatment) + geom_smooth(method = "lm") + 
       scale_color_manual(values = c("darkred", "darkorange2", "goldenrod2")) + 
       labs(y = "Insemination rate", x = "Strength") + 
       theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))

strength_mate_model <- lmer(data = attr_fem, mating_rate ~ strength*treatment + (1|block))
summary(strength_mate_model)
Anova(strength_mate_model)

plot(strength_mate_model)
plot(simulateResiduals(strength_mate_model))

## Mounting
ggplot(data = attr_fem, aes(y = mounts_in_rate, x = strength, color = block)) + 
       geom_point() + facet_grid(~treatment) + geom_smooth(method = "lm") + 
       scale_color_manual(values = c("darkred", "darkorange2", "goldenrod2")) + 
       labs(y = "Rate of mounts received", x = "Strength") + 
       theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))

strength_mount_model <- lmer(data = attr_fem, mounts_in_rate ~ strength*treatment + (1|block))
summary(strength_mount_model)
Anova(strength_mount_model)

plot(strength_mount_model)
plot(simulateResiduals(strength_mount_model))










