## setwd("C:/Users/janya/Desktop/R/bedbugs/bb_sna_2021")

library(tidyverse)
library(data.table)
library(netdiffuseR)
library(janitor)

## Loading in dataframe with all observations
all_data <- read.csv("data/bbsna_raw_combined.csv")
all_data <- remove_empty(all_data, which = c("rows", "cols"), quiet = TRUE)


## Getting all misspelled behaviours and IDs
sort(unique(all_data$behaviour)) 
sort(unique(all_data$focal_individual))
sort(unique(all_data$social_partner))
table_focals <- as.data.frame(table(all_data$focal_individual)) # Number of each mispelling 
table_partners <- as.data.frame(table(all_data$social_partner)) 


## Fixing mispelled behaviours
patch_table_behaviour <- read.csv("extra/patch_table.csv")

patch <- all_data %>% 
          left_join(patch_table_behaviour, by = "behaviour")

all_data <- patch %>% 
            mutate(behaviour = ifelse(is.na(patch_behaviour), behaviour, as.character(patch_behaviour))) %>% 
            select(-patch_behaviour)

sort(unique(all_data$behaviour))

## Fixing mispelled focal individuals
patch_table_focals <- read.csv("extra/ID_key_focal.csv")

patch_focals <- all_data %>%
                left_join(patch_table_focals, by = "focal_individual")

all_data <- patch_focals %>% 
            mutate(focal_individual = ifelse(is.na(patch_focal), "unknown", as.character(patch_focal))) %>% 
            select(-patch_focal)


## Fixing mispelled partner individuals
patch_table_partners <- read.csv("extra/ID_key_partner.csv")

patch_partners <- all_data %>%
                  left_join(patch_table_partners, by = "social_partner")

all_data <- patch_partners %>% 
            mutate(social_partner = ifelse(is.na(patch_partner), "unknown", as.character(patch_partner)))

mate_mount_data <- all_data %>%  
                   filter(behaviour == "mating" | behaviour == "mount")

sort(unique(mate_mount_data$focal_individual))
sort(unique(mate_mount_data$social_partner))


## CREATING MATING AND MOUNTING MATRICES
## Creating a function that turns data into edgelists and then into interaction matrices
func_mount_mat <- function(all_data, behav) {
                  all_data <- all_data %>% 
                              filter(behaviour == behav) %>% 
                              filter(focal_individual != "unknown") %>% 
                              select(c(focal_individual, patch_partner)) %>% 
                              mutate(edge_weight = 1)
                  mount_edgelist <- aggregate(data = all_data, edge_weight ~ focal_individual + patch_partner, FUN = sum)
                  mount_matrix <- edgelist_to_adjmat(mount_edgelist[1:2], w = mount_edgelist$edge_weight, 
                                                     undirected = FALSE)
return(as.matrix(mount_matrix))
}

## Applying function to all replicates and storing matrices as a list object
rep_list <- split(all_data, all_data$network)
mount_matrices <- lapply(rep_list, func_mount_mat, "mount") # Creates list of mounting matrices
mating_matrices <- lapply(rep_list, func_mount_mat, "mating") # Creates list of mating matrices

# Remove Q from rep 1 because she died after day 1
mount_matrices[[1]] <- mount_matrices[[1]][-17, -17]
mating_matrices[[1]] <- mating_matrices[[1]][-17, -17]

# Remove I from rep 5 because he had a deformity and could not mate
mount_matrices[[5]] <- mount_matrices[[5]][-9, -9]
mating_matrices[[5]] <- mating_matrices[[5]][-9, -9]

saveRDS(mount_matrices, "mount_matrices.rds")
saveRDS(mating_matrices, "mating_matrices.rds")


