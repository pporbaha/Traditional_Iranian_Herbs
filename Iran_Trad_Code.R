# Library
library(readxl)
library(tidyverse)
library(ggraph)
library(GGally)
library(igraph)
library(network)
library(minerva)
library(PhylogeneticEM)
library(AMR)
library(GOSemSim)
library(org.Hs.eg.db)
library(cowplot)
library(ggpubr)
library(kableExtra)
library(visNetwork)
library(webshot)
library(ChemmineR)
library(ggcorrplot)
library(CINNA)
library(rcdk)
library(fingerprint)
library(ChemmineR)
library(cowplot)
library(rcdk)

## Input of TCM herb dataset

df <- readxl::read_excel("our_data.xlsx") %>%
  dplyr::rename(drug=herb_id, target = ensymble_3, ingredient=`Ingredient name-DNP`, smile = `Canonical SMILE`) %>%
  mutate(drug = tolower(drug), ingredient = tolower(ingredient), cid=tolower(cid), target=toupper(target))

df <- df %>% dplyr::filter(properties_0 != "C/H")  # remove "C/H"

head(df)
str(df)

# remove digits from columns and make new columns
df$cleaned_properties <- gsub("[0-9]", "", df$Properties) 
df$cleaned_properties_0 <- gsub("[0-9]", "", df$properties_0)
df$cleaned_properties_1 <- gsub("[0-9]", "", df$properties_1)


# define two functions
count_unique <- function(x) length(unique(x))  # To measure number of unique values
count_na_percent <- function(x) round(sum(is.na(x))*100/length(x), digits = 5)  # To check percentage of NA values

df %>% dplyr::summarise_all(count_unique)
df %>% dplyr::summarise_all(count_na_percent)
df %>% dplyr::select(properties_0, properties_1) %>% table

term_freq <- df %>% dplyr::select(drug, Properties) %>% distinct 
unique(term_freq$drug) %>% length

tf_plot <- table(term_freq$Properties) %>% data.frame %>%  dplyr::rename(term = "Var1", freq = "Freq") %>% 
  mutate(term = as.character(term)) %>% arrange(desc(freq))

ggplot(tf_plot, aes(x = reorder(term, -freq), y = freq))+
  xlab("Properties")+
  theme_minimal() +
  geom_bar(stat = "identity", fill = "skyblue")+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),  # White background
    axis.text.x = element_text(size = 12, angle = 30, hjust = 0.5),
    legend.position = "none"
    )

# Target
term_freq <- df %>% dplyr::select(drug, target) %>% filter(!is.na(target)) %>% distinct 
unique(term_freq$drug) %>% length

tf_plot <- table(term_freq$target) %>% data.frame %>%  dplyr::rename(term = "Var1", freq = "Freq") %>% 
  mutate(term = as.character(term)) %>% arrange(desc(freq))

ggplot(tf_plot[1:20, ], aes(x = reorder(term, -freq), y = freq))+
  xlab("Targets")+
  geom_bar(stat = "identity", fill="pink")+
  theme(
    axis.text.x = element_text(size = 12, angle = 30),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none")

# Multipartite Network Analysis -------------------------------------------
## Reconstruct the UN-B graph drug and ingredient

# drug_ingredient <- df %>% dplyr::select(drug, ingredient) %>% distinct  ##? resulted 983 rows: the original code used cid instead of ingredient resulted 915 rows
drug_ingredient <- df %>% dplyr::select(drug, cid) %>% distinct

##Bipartite network
drug_ingredient_net <-  drug_ingredient %>% graph_from_data_frame(directed = FALSE)

number_drugs = df %>% distinct(drug) %>% nrow 
# number_ingredient = df %>% distinct(ingredient) %>% nrow
number_ingredient = df %>% distinct(cid) %>% nrow
total_number = number_drugs + number_ingredient

V(drug_ingredient_net)$type <- NA  # Pedram: Initialize type attribute for all vertices
V(drug_ingredient_net)$type[1:75] <- FALSE  # drugs 75 = df %>% distinct(drug) %>% nrow   or from summarize_all func
V(drug_ingredient_net)$type[76:852] <- TRUE # pubchem 993 = 75 + df %>% distinct(ingredient) %>% nrow

V(drug_ingredient_net)$color <- c("steel blue", "orange")[V(drug_ingredient_net)$type+1]
V(drug_ingredient_net)$shape <- c("square", "circle")[V(drug_ingredient_net)$type+1]

drug_ingredient_net <- giant_component_extract(drug_ingredient_net)[[1]]  ##? remove disconneted nodes or with little important connections to have better visualization.

### Visualization
plot.igraph(drug_ingredient_net, vertex.label = NA, vertex.size = (2 - V(drug_ingredient_net)$type) * 5, layout = layout_as_bipartite)

## Projection of drug similarity network (DSN) -----------------------------
drug_ingredient_incd <- as_biadjacency_matrix(drug_ingredient_net) # as_incidence_matrix() was deprecated and renamed to as_biadjacency_matrix() to create a more consistent API.
dim(drug_ingredient_incd) # result: 43 654 While it should be 75 799 (864 if ingredients) if there is not giant_component_extract 
# drug_ingredient_incd %>% data.frame %>% View

# Drugs
proj_drug_sim <- drug_ingredient_incd %*% t(drug_ingredient_incd)
dim(proj_drug_sim)
summary(as.vector(proj_drug_sim) > 0)
# proj_drug_sim %>% View()

(proj_drug_sim_net <- graph_from_adjacency_matrix(proj_drug_sim, mode = "upper", weighted = TRUE, diag = FALSE))
igraph::components(proj_drug_sim_net)$no

(E(proj_drug_sim_net)[weight > 0] %>% length) * 2 + (dim(proj_drug_sim)[1]) == summary(as.vector(proj_drug_sim) > 0)[3] # Check the number of edges

# Visualization
# Drugs
V(proj_drug_sim_net)$color <- "orange"
E(proj_drug_sim_net)$width <- 1
E(proj_drug_sim_net)$color <- "gray"
E(proj_drug_sim_net)[weight > median(E(proj_drug_sim_net)$weight)]$width <- 4
E(proj_drug_sim_net)[weight > median(E(proj_drug_sim_net)$weight)]$color <- "blue"

g5 <- ggnet2(proj_drug_sim_net, color = "red", size = 1)+
  theme(legend.position = "none")
g5

## Pedram just visualize more than 9 connections
# # Step 1: Filter edges with weight greater than 9
# edges_to_keep <- E(proj_drug_sim_net)[E(proj_drug_sim_net)$weight > 9]
# 
# # Step 2: Create a subgraph with the filtered edges
# filtered_net <- subgraph.edges(proj_drug_sim_net, edges_to_keep)
# 
# # Step 3: Visualize the filtered subgraph
# # plot(filtered_net, vertex.label = NA, vertex.size = 5, edge.width = E(filtered_net)$weight / 5, 
# #      edge.color = "gray", main = "Drug Similarity Network (Weight > 20)")
# g_drugs_just_more_weighted <- ggnet2(filtered_net, 
#                                      color = "blue",   # Color of the vertices
#                                      size = 3,        # Size of the vertices
#                                      edge.size = E(filtered_net)$weight / 5,  # Scale edge width
#                                      edge.color = "gray") + 
#   theme(legend.position = "none") +
#   ggtitle("Drug Similarity Network (Weight > 9)")
# g_drugs_just_more_weighted

### Projection of ingredient similarity network (ISN) -----------------------
# Ingredients
proj_ingd_sim <- t(drug_ingredient_incd)  %*% drug_ingredient_incd
dim(proj_ingd_sim)
summary(as.vector(proj_ingd_sim) > 0)

(proj_ingd_sim_net <- graph_from_adjacency_matrix(proj_ingd_sim, mode = "upper", weighted = TRUE, diag = FALSE))
igraph::components(proj_ingd_sim_net)$no

(E(proj_ingd_sim_net)[weight > 0] %>% length) * 2 + (dim(proj_ingd_sim)[1]) == summary(as.vector(proj_ingd_sim) > 0)[3] # Check the number of edges

### Visualization
# ingds
V(proj_ingd_sim_net)$color <- "blue"
E(proj_ingd_sim_net)$width <- 1
E(proj_ingd_sim_net)$color <- "gray"
E(proj_ingd_sim_net)[weight > median(E(proj_ingd_sim_net)$weight)]$width <- 4
E(proj_ingd_sim_net)[weight > median(E(proj_ingd_sim_net)$weight)]$color <- "orange"

g6 <- ggnet2(proj_drug_sim_net, color = "blue", size = 1)+
  theme(legend.position = "none")
g6

## Modularity analysis of DSN
## Drugs modularity analysis

#cluster_fast_greedy cluster_infomap	cluster_spinglass but cluster_fast_greedy resulted most modularity score in supp Table1
drug_cluster_net <- cluster_fast_greedy(proj_drug_sim_net, weights = E(proj_drug_sim_net)$weight)
modularity(drug_cluster_net)
drug_cluster_net <- cluster_infomap(proj_drug_sim_net)
modularity(drug_cluster_net)
# drug_cluster_net <- cluster_spinglass(proj_drug_sim_net, spins = 25)  # You can adjust spins to fit your problem.
# modularity(drug_cluster_net)

# The best algorithm with the best modularity score
drug_cluster_net <- cluster_fast_greedy(proj_drug_sim_net, weights = E(proj_drug_sim_net)$weight)
modularity(drug_cluster_net)
sizes(drug_cluster_net) %>% as.data.frame

mem_med_d <- data.frame("member_dsn" = unclass(membership(drug_cluster_net))) %>% rownames_to_column(var = "drug")
drug_ingredient <- drug_ingredient %>% left_join(mem_med_d, by= "drug") %>% filter(!is.na(member_dsn))
# head(drug_ingredient)

## Modularity analysis of ISN
## Ingredients modularity analysis
ingd_cluster_net <- cluster_fast_greedy(proj_ingd_sim_net, weights = E(proj_ingd_sim_net)$weight) 
modularity(ingd_cluster_net)

ingd_cluster_net <- cluster_infomap(proj_ingd_sim_net)
modularity(ingd_cluster_net)

# ingd_cluster_net <- cluster_spinglass(proj_ingd_sim_net, spins = 25)
# modularity(ingd_cluster_net)

# Best Algorithm according to modularity score 
ingd_cluster_net <- cluster_fast_greedy(proj_ingd_sim_net, weights = E(proj_ingd_sim_net)$weight) 
modularity(ingd_cluster_net)
sizes(ingd_cluster_net) %>% as.data.frame

mem_med_i <- data.frame("member_isn" = unclass(membership(ingd_cluster_net))) %>% rownames_to_column(var = "cid")

drug_ingredient <- drug_ingredient %>% left_join(mem_med_i, by= "cid") %>% filter(!is.na(member_isn))

# head(drug_ingredient)

### Visualization
ggplot(sizes(drug_cluster_net) %>% as.data.frame, 
       aes(x = reorder(Community.sizes, -Freq), y = Freq, fill = Community.sizes))+
  geom_bar(stat = "identity")+
  guides(fill = "none")+
  geom_text(aes(label = Freq), vjust = -0.5, color = "black", size = 5)+
  ylab("Frequency")+
  xlab("Cluster")+
  ggtitle('drug cluster net')+
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))  # Centered title with larger font


ggplot(sizes(ingd_cluster_net) %>% as.data.frame, 
       aes(x = reorder(Community.sizes, -Freq), y = Freq, fill = Community.sizes))+
  geom_bar(stat = "identity")+
  guides(fill = "none")+
  geom_text(aes(label = Freq), vjust = -0.5, color = "black", size = 4)+
  ylab("Frequency")+
  xlab("Cluster")+
  ggtitle('ingredient cluster net')+
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))  # Centered title with larger font

# ------------
### Computational Validation
### Drug's cluster
### 1. No. of common Meridian in each drug cluster


# # Function to calculate the median of average similarities between drugs in a cluster
# median_average_similarity_in_each_cluster <- function(chr1_chr2_mem, mem_filter, property) {
#   # Computes the median of the average similarities of the specified property 
#   #
#   # Args: 
#   #   chr1_chr2_mem: A dataframe to be converted to a list by cluster membership
#   #   mem_filter: A number to filter the membership cluster
#   #   property: The column name containing the property to compare
#   # Returns:
#   #   The median of the average of all similarities
#   
#   # Filter and split the data by drug
#   chr1_chr2_mem_list <- chr1_chr2_mem %>%
#     filter(member_dsn == mem_filter) %>%
#     split(.$drug)
#   
#   # If the cluster contains just one drug, return NA
#   if (length(chr1_chr2_mem_list) < 2) {
#     return(NA)
#   }
#   
#   total_comparisons <- c()
# 
#   # Calculate pairwise comparisons
#   for (i in seq_along(chr1_chr2_mem_list)) {
#     for (j in seq_along(chr1_chr2_mem_list)) {
#       if (i != j) {  # Avoid self-comparison
#         if (chr1_chr2_mem_list[[i]][[property]] == chr1_chr2_mem_list[[j]][[property]]) {
#           total_comparisons <- c(total_comparisons, 1)
#         }
#         else {
#           total_comparisons <- c(total_comparisons, 0)
#         }
#       }
#     }
#   }
# 
#   # Return the median of the average similarities
#   print(total_comparisons)
#   return(mean(total_comparisons, na.rm = TRUE))
# }
# 
# property_names <- c("Properties", "properties_0", "properties_1", "cleaned_properties_0", "cleaned_properties_1")
# for (property in property_names) {
#   
#   (tcm_merid <- df %>% dplyr::select(drug, !!sym(property)) %>% distinct(drug, !!sym(property)))
#   
#   tcm_merid_cl <- tcm_merid %>% inner_join(drug_ingredient, by = "drug") %>% dplyr::select(drug, !!sym(property), member_dsn) %>% distinct_all 
#   
#   
#   all_med_merid <- data.frame()
#   for(i in 1:6){
#     all_med_merid[i, 1] <- median_average_similarity_in_each_cluster(tcm_merid_cl, i, property)
#   }
#   # Profile of meridian in each cluster of DSN
#   cluster_size <- tcm_merid_cl %>% group_by(member_dsn) %>% dplyr::summarise(count = n())
#   
#   merid_mat <- tcm_merid_cl %>% group_by(member_dsn, !!sym(property)) %>% dplyr::summarise(n = n()) %>% inner_join(cluster_size, by = "member_dsn") %>% arrange(member_dsn, desc(n)) %>%  mutate(freq = signif(n / count, digits = 2)) %>% dplyr::select(-n, -count) %>% spread(key = !!sym(property), value = freq)
#   file_name = paste(property, "_mat.csv")
#   write_csv(merid_mat, file_name)
#   
#   #Random selection 
#   inters_mean <- c()
#   
#   for(k in 1:6){
#     total_comparisons <- c()
#     sample_num <- table(tcm_merid_cl$member_dsn) %>% sample(1)
#     tcm_merid_cl_list <- tcm_merid_cl %>% sample_n(sample_num, replace = FALSE) %>% split(.$drug)
#     
#     for (i in seq_along(tcm_merid_cl_list)){
#       for(j in seq_along(tcm_merid_cl_list)){
#         if (i != j) {  # Avoid self-comparison
#           if (tcm_merid_cl_list[[i]][[property]] == tcm_merid_cl_list[[j]][[property]]) {
#             total_comparisons <- c(total_comparisons, 1)
#           }
#           else {
#             total_comparisons <- c(total_comparisons, 0)
#           }
#         }
#       }
#     }
#     print(total_comparisons)
#     inters_mean[k] <- mean(total_comparisons, na.rm = TRUE)
#   }
#   
#   #### Visualization
#   intersect_merid <- data.frame(intersection = c(inters_mean,
#                                                  all_med_merid$V1),
#                                 compare = c(rep("random", times = length(inters_mean)),
#                                             rep("clusters", times = dim(all_med_merid)[1])))
#   
#   intersect_merid_clean <- na.omit(intersect_merid)
#   
#   print(property)
#   print(intersect_merid)
#   plot_title <- paste(property, "similarity plot")
#   
#   g1 <- ggplot(intersect_merid, aes(x = intersection, fill = compare, col = compare))+
#     geom_density(position="identity", alpha = 0.5, linewidth = 1)+ xlab(plot_title)+
#     # geom_vline(xintercept = c(median(inters_mean), median(all_med_merid$V1)),
#     #            color = c("darkblue", "red"), linetype = "dashed", size = 1) +
#     theme(legend.position = "top")
#   print(g1)
# }

#---------------
### property : both properties_0 and propersties_1
df_both <- df %>%
  pivot_longer(cols = starts_with("properties_"),
               names_to = "property_type",
               values_to = "both_properties")

(tcm_merid <- df_both %>% dplyr::select(drug, both_properties) %>% distinct(drug, both_properties))

tcm_merid_cl <- tcm_merid %>% inner_join(drug_ingredient, by = "drug") %>% dplyr::select(drug, both_properties, member_dsn) %>% distinct_all 

intersect_all <- function(chr1_chr2_mem, mem_filter){
  # Computes the intersection of multiple list and the average of intersections
  #
  # Args: 
  #   chr1_chr2_mem: A dataframe to converted to list by cluster membership No.
  #   mem_filter: A number to filter the membership cluster
  # Returns:
  #   The median of the average of all intersections
  #   
  
  chr1_chr2_mem_list <- chr1_chr2_mem %>% filter(member_dsn == mem_filter) %>% split(list(.$drug))
  inters_cl <- matrix(0, length(chr1_chr2_mem_list), length(chr1_chr2_mem_list))
  
  # If the cluster contains just one drug, return NA
  if (length(chr1_chr2_mem_list) < 2) {
    return(NA)
  }
  
  for (i in seq_along(chr1_chr2_mem_list)){
    for(j in seq_along(chr1_chr2_mem_list)){
      
      inters_cl[i, j] <- length(intersect(chr1_chr2_mem_list[[i]]$both_properties, chr1_chr2_mem_list[[j]]$both_properties))
    }
    diag(inters_cl) <- 0
    inters_cl <- inters_cl %>% data.frame %>% mutate(average = rowSums(.)/length(chr1_chr2_mem_list))
  }
  return(median(inters_cl$average))
}

all_med_merid <- data.frame()
for(i in 1:6){
  all_med_merid[i, 1] <- intersect_all(tcm_merid_cl, i)
}
# Profile of meridian in each cluster of DSN
cluster_size <- tcm_merid_cl %>% group_by(member_dsn) %>% dplyr::summarise(count = n())

merid_mat <- tcm_merid_cl %>% group_by(member_dsn, both_properties) %>% dplyr::summarise(n = n()) %>% inner_join(cluster_size, by = "member_dsn") %>% arrange(member_dsn, desc(n)) %>%  mutate(freq = signif(n / count, digits = 2)) %>% dplyr::select(-n, -count) %>% spread(key = both_properties, value = freq)

write_csv(merid_mat, "property_both_mat.csv")

#Random selection 
inters_mean <- c()
inters_list <- list()
for(k in 1:6){
  sample_num <- table(tcm_merid_cl$member_dsn) %>% sample(1)
  tcm_merid_cl_list <- tcm_merid_cl %>% sample_n(sample_num, replace = FALSE) %>% split(list(.$drug))
  
  inters <- matrix(0, length(tcm_merid_cl_list), length(tcm_merid_cl_list))
  for (i in seq_along(tcm_merid_cl_list)){
    for(j in seq_along(tcm_merid_cl_list)){
      
      inters[i, j] <- length(intersect(tcm_merid_cl_list[[i]]$both_properties, tcm_merid_cl_list[[j]]$both_properties))
    }
    diag(inters) <- 0
    inters <- inters %>% data.frame %>% mutate(average = rowSums(.)/length(tcm_merid_cl_list))
  }
  inters_list[k] <- list(inters)
  inters_mean[k] <- median(inters$average, na.rm=TRUE)
}

# Visualization
intersect_merid <- data.frame(intersection = c(inters_mean,
                                               all_med_merid$V1),
                              compare = c(rep("random", times = length(inters_mean)),
                                          rep("clusters", times = dim(all_med_merid)[1])))
print('Both property')

intersect_merid_clean <- na.omit(intersect_merid)

g1 <- ggplot(intersect_merid, aes(x = intersection, fill = compare, col = compare))+
  geom_density(position="identity", alpha = 0.5, linewidth = 1)+ xlab("Meridian similarity")+
  geom_vline(xintercept = c(median(inters_mean, na.rm=TRUE), median(all_med_merid$V1, na.rm=TRUE)),
             color = c("darkblue", "red"), linetype = "dashed", size = 1)+
  theme(legend.position = "top")  +xlab("Both properties similarity plot")
print(g1)
ggsave("properties_both_similarity_plot.jpg", plot = g1,
       width = 10, height = 4,
       dpi = 300)
print(intersect_merid)

# --------------
#### involving tempers property degree  #Hot4 = +4  Cold3 = -3 Moderate =0
df <- df %>%
  dplyr::mutate(properties_0_number = case_when(
    grepl("Hot", properties_0) ~ as.numeric(gsub("Hot", "", properties_0)),  # Hot -> +
    grepl("Cold", properties_0) ~ -as.numeric(gsub("Cold", "", properties_0)), # Cold -> -
    properties_0 == "Moderate" ~ 0  # Moderate -> 0
  ))

df <- df %>%
  dplyr::mutate(properties_1_number = case_when(
    grepl("Dry", properties_1) ~ as.numeric(gsub("Dry", "", properties_1)),   # Dry -> +
    grepl("Wet", properties_1) ~ -as.numeric(gsub("Wet", "", properties_1))   # Wet -> -
  ))


# Function to calculate the median of average similarities between drugs in a cluster
median_average_similarity__include_temper_intensity_in_each_cluster <- function(chr1_chr2_mem, mem_filter) {
  # Computes the median of the average similarities of the specified property 
  #
  # Args: 
  #   chr1_chr2_mem: A dataframe to be converted to a list by cluster membership
  #   mem_filter: A number to filter the membership cluster
  #   property: The column name containing the property to compare
  # Returns:
  #   The median of the average of all similarities
  
  # Filter and split the data by drug
  chr1_chr2_mem_list <- chr1_chr2_mem %>%
    filter(member_dsn == mem_filter) %>%
    split(.$drug)
  
  # If the cluster contains just one drug, return NA
  if (length(chr1_chr2_mem_list) < 2) {
    return(NA)
  }
  
  total_comparisons <- c()
  
  # Calculate pairwise comparisons
  for (i in seq_along(chr1_chr2_mem_list)) {
    for (j in seq_along(chr1_chr2_mem_list)) {
      if (i != j) {  # Avoid self-comparison
        # Calculate the absolute differences and append to total_comparisons
        difference <- abs(chr1_chr2_mem_list[[i]]$properties_0_number - chr1_chr2_mem_list[[j]]$properties_0_number) + 
          abs(chr1_chr2_mem_list[[i]]$properties_1_number - chr1_chr2_mem_list[[j]]$properties_1_number)
        # difference <- difference ** 2
        total_comparisons <- c(total_comparisons, difference)
      }
    }
  }
  
  
  # Return the median of the average similarities
  # print(total_comparisons)
  return(median(total_comparisons, na.rm = TRUE))
}

tcm_merid <- df %>% dplyr::select(drug, properties_0_number, properties_1_number) %>% distinct_all

tcm_merid_cl <- tcm_merid %>% inner_join(drug_ingredient, by = "drug") %>% dplyr::select(drug, properties_0_number, properties_1_number, member_dsn) %>% distinct_all 

all_med_merid <- data.frame()
for(i in 1:6){
  all_med_merid[i, 1] <- median_average_similarity__include_temper_intensity_in_each_cluster(tcm_merid_cl, i)
}
# Profile of meridian in each cluster of DSN
cluster_size <- tcm_merid_cl %>% group_by(member_dsn) %>% dplyr::summarise(count = n())

# merid_mat <- tcm_merid_cl %>% group_by(member_dsn, properties_0_number, properties_1_number) %>%
#   dplyr::summarise(n = n()) %>% inner_join(cluster_size, by = "member_dsn") %>% arrange(member_dsn, desc(n)) %>%
#   mutate(freq = signif(n / count, digits = 2)) %>% dplyr::select(-n, -count) %>%
#   spread(key = !!sym(property), value = freq)
# file_name = paste("Properties_numbers_mat.csv")
# write_csv(merid_mat, file_name)

#Random selection 
inters_mean <- c()

for(k in 1:6) {
  total_comparisons <- c()
  sample_num <- table(tcm_merid_cl$member_dsn) %>% sample(1)
  tcm_merid_cl_list <- tcm_merid_cl %>% sample_n(sample_num, replace = FALSE) %>% split(.$drug)
  
  for (i in seq_along(tcm_merid_cl_list)) {
    for (j in seq_along(tcm_merid_cl_list)) {
      if (i != j) {  # Avoid self-comparison
        # Calculate the absolute differences and append to total_comparisons
        difference <- abs(tcm_merid_cl_list[[i]]$properties_0_number - tcm_merid_cl_list[[j]]$properties_0_number) + 
          abs(tcm_merid_cl_list[[i]]$properties_1_number - tcm_merid_cl_list[[j]]$properties_1_number)
        # difference <- difference ** 2
        total_comparisons <- c(total_comparisons, difference)
      }
    }
  }
  
  # print(total_comparisons)
  if (is.null(total_comparisons)) {
    inters_mean[k] <- NA
  } else {
  inters_mean[k] <- median(total_comparisons, na.rm = TRUE)
  }
}

# Visualization
intersect_merid <- data.frame(intersection = c(inters_mean,
                                               all_med_merid$V1),
                              compare = c(rep("random", times = length(inters_mean)),
                                          rep("clusters", times = dim(all_med_merid)[1])))

intersect_merid_clean <- na.omit(intersect_merid)

print("Properties Numbers Difference plot")
plot_title <- paste( "Properties Numbers Difference plot")

g1 <- ggplot(intersect_merid, aes(x = intersection, fill = compare, col = compare))+
  geom_density(position="identity", alpha = 0.5, linewidth = 1)+ xlab(plot_title)+
   geom_vline(xintercept = c(median(inters_mean, na.rm=TRUE), median(all_med_merid$V1, na.rm=TRUE)),
              color = c("darkblue", "red"), linetype = "dashed", size = 1) +
  theme(legend.position = "top") +
  scale_x_continuous(limits = c(0, NA))  # Set x-axis limit starting from 0
  
ggsave("properties_difference_plot.jpg", plot = g1,
       width = 10, height = 4,
       dpi = 300)
print(g1)
print(intersect_merid)

# ------------------------------------------------------
## Ingredient's cluster
### 3. No. of common Targets in each ingredient cluster
tcm_target <- df %>% dplyr::select(cid, target) %>% distinct_all

tcm_target_cl <- tcm_target %>% inner_join(drug_ingredient, by = "cid") %>% dplyr::select(cid, target, member_isn) %>% distinct_all  

intersect_all <- function(chr1_chr2_mem, mem_filter){
  # Computes the intersection of multiple list and the average of intersections
  #
  # Args: 
  #   chr1_chr2_mem: A dataframe to converted to list by cluster membership No.
  #   mem_filter: A number to filter the membership cluster
  # Returns:
  #   The median of the average of all intersections
  #   
  
  chr1_chr2_mem_list <- chr1_chr2_mem %>% filter(member_isn == mem_filter) %>% split(list(.$cid))
  inters_cl <- matrix(0, length(chr1_chr2_mem_list), length(chr1_chr2_mem_list))
  
  for (i in seq_along(chr1_chr2_mem_list)){
    for(j in seq_along(chr1_chr2_mem_list)){
      
      inters_cl[i, j] <- length(intersect(chr1_chr2_mem_list[[i]]$target, chr1_chr2_mem_list[[j]]$target))
    }
    diag(inters_cl) <- 0
    inters_cl <- inters_cl %>% data.frame %>% mutate(average = rowSums(.)/length(chr1_chr2_mem_list))
  }
  return(median(inters_cl$average))
}

all_med_target <- data.frame()
for(i in 1:17){
  all_med_target[i, 1] <- intersect_all(tcm_target_cl, i)
}

table(tcm_target_cl$member_isn) %>% mean

#Random selection 
inters_mean <- c()
inters_list <- list()
for(k in 1:17){
  if (k%%5 ==0) {print(k)}
  sample_num <- table(tcm_target_cl$member_isn) %>% sample(1)
  tcm_target_cl_list <- tcm_target_cl %>% sample_n(sample_num, replace = FALSE) %>% split(list(.$cid))
  
  inters <- matrix(0, length(tcm_target_cl_list), length(tcm_target_cl_list))
  for (i in seq_along(tcm_target_cl_list)){
    for(j in seq_along(tcm_target_cl_list)){
      
      inters[i, j] <- length(intersect(tcm_target_cl_list[[i]]$target, tcm_target_cl_list[[j]]$target))
    }
    diag(inters) <- 0
    inters <- inters %>% data.frame %>% mutate(average = rowSums(.)/length(tcm_target_cl_list))
  }
  inters_list[k] <- list(inters)
  inters_mean[k] <- mean(inters$average)
}

# Visualization
intersect_target <- data.frame(intersection = c(inters_mean,
                                                all_med_target$V1),
                               compare = c(rep("random", times = length(inters_mean)),
                                           rep("clusters", times = dim(all_med_target)[1])))

g3 <- ggplot(intersect_target, aes(x = intersection, fill = compare, col = compare))+
  geom_density(position="identity", alpha = 0.5, linewidth = 1)+ 
  ylim(0, 4)+
  xlab("Targets similarity")+
  geom_vline(xintercept = c(median(filter(intersect_target, compare == "random")$intersection, na.rm = T), median(all_med_target$V1)),
             color = c("darkblue", "red"), linetype = "dashed", size = 1)+
  theme(legend.position = "top")+ xlim(0,1)

print(g3)

# -------------------------------------------------------
### 4. Similar SMILEs in each ingredient cluster

(tcm_smile <- df %>% dplyr::select(cid, smile) %>% distinct_all)

tcm_smile_cl <- tcm_smile %>% inner_join(drug_ingredient, by = "cid") %>% dplyr::select(cid, smile, member_isn) %>% distinct_all 

finger <- function(chr1_chr2_mem, mem_filter){
  # Computes the median of fingerprint similarity of multiple list of molecules
  #
  # Args: 
  #   chr1_chr2_mem: A dataframe to be filtered.
  #   mem_filter: A number to filter the membership cluster No.
  # Returns:
  #   The median of all Tanimoto simialrities
  #   
  chr1_chr2_mem_filter <- chr1_chr2_mem %>% filter(member_isn == mem_filter)
  
  mols <- parse.smiles(chr1_chr2_mem_filter[ , 2, drop = TRUE])
  fps <- lapply(mols, get.fingerprint, type='extended', size = 500)
  fp.sim <- fingerprint::fp.sim.matrix(fps, method = 'dice')
  diag(fp.sim) <- 0
  med_fing <- median(fp.sim)
  return(med_fing)
}

all_med_smile <- data.frame()
for(i in 1:17){
  all_med_smile[i, 1] <- finger(tcm_smile_cl, i)
}
table(tcm_smile_cl$member_isn)

#Random selection 
sim_mean <- c()

for(k in 1:20){
  if (k%%5 == 0) {print(k)}
  
  sample_num <- table(tcm_smile_cl$member_isn) %>% sample(1)
  rand_smiles <- tcm_smile_cl %>% sample_n(sample_num, replace = FALSE)
  mols <- parse.smiles(rand_smiles[ , 2, drop = TRUE])
  fps <- lapply(mols, get.fingerprint, type='extended', size = 500)
  fp.sim <- fingerprint::fp.sim.matrix(fps, method = 'dice')
  
  diag(fp.sim) <- 0
  sim_mean[k] <- median(fp.sim)
}



#Visualization
intersect_smile <- data.frame(intersection = c(sim_mean,
                                               all_med_smile$V1),
                              compare = c(rep("random", times = length(sim_mean)),
                                          rep("clusters", times = dim(all_med_smile)[1])))

g4 <- ggplot(intersect_smile, aes(x = intersection, fill = compare, col = compare))+
  geom_density(position="identity", alpha = 0.5, linewidth = 1)+ xlab("Dice similarity")+
  geom_vline(xintercept = c(median(sim_mean), median(all_med_smile$V1)),
             color = c("darkblue", "red"), linetype = "dashed", size = 1)+
  theme(legend.position = c(0.8, 0.7))
print(g4)

# -------
# All figures in One
plot_grid(g1, g3, g4, align = "hv", 
          nrow = 2, labels = LETTERS[1:3])


