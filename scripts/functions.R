## R-script for all functions created

## FUNCTION 1: Turns raw data from bbsna_aggregations into igraph objects (one per replicate)
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

## FUNCTION 2: Turns raw data from bbsna_aggregations into ibi matrices (one per replicate)
# Useful for feeding into assortment.discrete 
func_ibi <- function(rep_groups){
  group_list <- strsplit(rep_groups$Members, " ")
  gbi_matrix <- get_group_by_individual(group_list, data_format = "groups")
  ibi_matrix <- get_network(gbi_matrix, data_format = "GBI")
  ibi_matrix <- ibi_matrix[order(rownames(ibi_matrix)) , order(colnames(ibi_matrix))] # alphabetical order
  return(ibi_matrix)
}

## FUNCTION 3: Assigns additional attributes (size, replicate, treatment) to each node on each igraph object
# Input: Takes a lists of igraph objects
# Ourput: Creates a dataframe of node attributes

func_attr <- function(igraph_objects){
  new_attr <- data.frame()
  for (i in 1:length(igraph_objects)){
    attr_i <- subset(attr, network == i & notes != "died" & notes != "weird")
    igraph_objects[[i]] <- set_vertex_attr(igraph_objects[[i]], "size", value = attr_i$thorax.mm)
    igraph_objects[[i]] <- set_vertex_attr(igraph_objects[[i]], "treatment", value = attr_i$treatment)
    igraph_objects[[i]] <- set_vertex_attr(igraph_objects[[i]], "network", value = attr_i$network)
    igraph_objects[[i]] <- set_vertex_attr(igraph_objects[[i]], "block", value = attr_i$block)
    (new_attr <- rbind(new_attr, vertex_attr(igraph_objects[[i]])))
  }
  return(new_attr)
}

## FUNCTION 4: Visualizing the social networks (this is useful for detecting errors in other parts of the code)
# Input: An igraph object
# Output: The SNA graph where size = strength*10
func_plot_network <- function(igraph_object, node_size){
  V(igraph_object)$color <- ifelse(V(igraph_object)$sex == "Female", "sandybrown", "skyblue3")
  V(igraph_object)$size <- V(igraph_object)$strength*node_size
  V(igraph_object)$label.color <- "white"
  E(igraph_object)$width <- E(igraph_object)$weight*6
  plot(igraph_object, edge.color = "dimgrey") #, vertex.label = NA) allows us to turn on/off vertex labels
}

## FUNCTION 5: Shuffle node labels
# Input: aggregations raw data
# Output: igraph objects with nodes randomized and sexes assigned
func_permute_igraph <- function(rep_list_group) { 
  group_list <- strsplit(rep_list_group$Members, " ")
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

## FUNCTION 6: Assigns the rest of the node attributes to the randomized nodes using left_join; runs the glm
# Input: a list of randomized igraph objects
# Steps: creates a dataframe combining the attributes from the different igraph objects
       # then runs the glm using this new dataframe 
# Output: the coefficient from the glm for effect of sex on strength
func_random_model_p1 <- function(random_igraphs){
  sim_attr <- data.frame()
  for (i in 1:length(random_igraphs)){
      attr_i <- attr %>% 
      filter(network == i) %>% 
      rename("name" = "ID") %>% 
      select(c("name", "thorax.mm", "network", "treatment", "block"))
      
      new_attr <- as.data.frame(vertex_attr(random_igraphs[[i]])) %>% 
                  left_join(attr_i, by = "name")
      sim_attr <- rbind(sim_attr, new_attr)
  }
  sim_model <- glmer(strength ~ sex*treatment + (1|block) + (1|thorax.mm), data = sim_attr, family = Gamma(link="log"))
  sim_coef_mat <- as.matrix(fixef(sim_model))
  sim_coef <- sim_coef_mat[2]
  return(sim_coef)
}

## FUNCTION 7: Assortativity permutation
# Input: An ibi matrix
# Output: Histogram of permutation + the observed assortativity score + the p-value
func_permute_assort <- function(ibi_matrix){
     sex_table <- as.data.frame(colnames(ibi_matrix)) %>% 
                  rename("ID" = "colnames(ibi_matrix)") %>% 
                  mutate(sex = ifelse(ID %in% LETTERS[1:12], "Male", "Female"))
     # Getting observed assortativity score
     obs_assort_index <- assortment.discrete(ibi_matrix, 
                                             types = sex_table$sex, weighted = TRUE)$r    
     n_sim_2 <- 999 # Setting up the permutation
     set.seed(33)
     sim_assort_index <- numeric(n_sim_2)
     
     for (i in 1:n_sim_2){ #Runs the permutation
        new_names <- sample(colnames(ibi_matrix))
        colnames(ibi_matrix) <- new_names
        rownames(ibi_matrix) <- new_names
    
        sex_table_new <- as.data.frame(colnames(ibi_matrix)) %>% 
                         rename("ID" = "colnames(ibi_matrix)") %>% 
                         mutate(sex = ifelse(ID %in% LETTERS[1:12], "Male", "Female"))

        sim_assort_index[i] <- assortment.discrete(ibi_matrix, types = sex_table_new$sex, 
                                                   weighted = TRUE)$r  
     }
     if (obs_assort_index >= mean(sim_assort_index)) { #Computes the p-value
         p <- 2*mean(sim_assort_index >= obs_assort_index) } else {
         p <- 2*mean(sim_assort_index <= obs_assort_index)
     }
     sim_assort_index <- c(sim_assort_index, obs_assort_index)
     list <- list("p-value" = p, "observed assortativity score" = obs_assort_index)
     hist(sim_assort_index, breaks = 29, xlim = c(min(sim_assort_index), max(sim_assort_index)), 
          ylim = c(0, 180), col = "aliceblue")
     lines(x = c(obs_assort_index, obs_assort_index), y = c(0, 150), col = "red", lty = "dashed", lwd = 2)
     return(list)
}

## FUNCTION 8: Calculates # of matings (equivalent to in-strength of mating network)
               # Also converts matrix to igraph object and plots the network bc why not 
# Input: A mating matrix
# Output: An igraph matrix with node attributes "matings" and "sex" assigned and the SNA graph
func_matrix_to_igraph <- function(matrix){
  igraph <- graph_from_adjacency_matrix(matrix, diag = FALSE, weighted = TRUE, mode = "directed")
  igraph <- set_vertex_attr(igraph, "sex", 
                            value = ifelse(V(igraph)$name %in% LETTERS[1:12], "Male", "Female"))
  strength <- strength(igraph, mode = "in")
  igraph <- set_vertex_attr(igraph, "matings", value = strength)
  V(igraph)$color <- ifelse(V(igraph)$sex == "Female", "red", "blue")
  V(igraph)$label.color <- "white"
  V(igraph)$size <- V(igraph)$matings*3.5
  E(igraph)$width <- E(igraph)$weight*1.5
  plot(igraph, edge.color = "dimgrey", layout = layout_nicely(igraph))
  return(igraph)
}

## FUNCTION 9: Creates igraph objects where only the female nodes are shuffled among themselves
# Input: A matrix
# Output: A randomized igraph object
func_permute_igraph_females <- function(matrix) { 
  #shuffle names 
  names <- colnames(matrix)
  m_names <- subset(names, names %in% LETTERS[1:12])
  f_names <- sample(subset(names, names %in% LETTERS[13:24]))
  new_names <- c(m_names, f_names)
  colnames(matrix) <- new_names; rownames(matrix) <- new_names
  
  igraph <- graph_from_adjacency_matrix(matrix, diag = FALSE, weighted = TRUE, mode = "directed")
  igraph <- set_vertex_attr(igraph, "sex", 
                            value = ifelse(V(igraph)$name %in% LETTERS[1:12], "Male", "Female"))
  strength <- strength(igraph, mode = "in")
  igraph <- set_vertex_attr(igraph, "matings", value = strength)
  return(igraph)
}  

## FUNCTION 10: Assigns the rest of the node attributes to the randomized nodes using left_join; 
                # runs the glm for Prediction 3
# Input: a list of randomized igraph objects
# Steps: creates a dataframe combining the attributes from the different igraph objects
# then runs the glm using this new dataframe 
# Output: the coefficient from the glm for the relationship between strength and matings
func_random_model_p3 <- function(random_igraphs){
  sim_attr <- data.frame()
  for (i in 1:length(random_igraphs)){
    attr_i <- attr_observed %>% 
      filter(network == i) %>% 
      select(c("name", "size", "network", "treatment", "strength", "block"))
    
    new_attr <- as.data.frame(vertex_attr(random_igraphs[[i]])) %>% 
                left_join(attr_i, by = "name")
    new_attr <- new_attr %>% 
                filter(sex == "Female")
    sim_attr <- rbind(sim_attr, new_attr)
  }
  sim_model <- glm(matings ~ strength + size + network + treatment +
                   block, data = sim_attr, family = Gamma(link="log"))
  return(coef(sim_model)[2])
}
