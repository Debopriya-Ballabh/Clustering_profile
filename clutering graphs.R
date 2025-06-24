# Load libraries
BiocManager::install("phyloseq", force = TRUE)
library("phyloseq")
library("ggplot2")
library("vegan")

# Set working directory
setwd("C:/Users/Debopriya/Downloads/beta_diversity_NIWASM")
getwd()

# Step 1: Read the metadata
metadata <- read.csv("C:/Users/Debopriya/Downloads/beta_diversity_NIWASM/metadata_Mitali.csv", sep = ",", header = TRUE, row.names = 1)

# Step 2: Load phyloseq object
ps <- readRDS("C:/Users/Debopriya/Downloads/beta_diversity_NIWASM/ps.rds")
ps_tss <- transform_sample_counts(ps, function(x) x / sum(x))
sample_names(ps_tss) <- rownames(metadata)
sample_data(ps_tss) <- sample_data(metadata)
ps_tss@sam_data$Groups <- factor(ps_tss@sam_data$Groups, levels = unique(ps_tss@sam_data$Groups))
ps_tss@sam_data$Groups
str(sample_data(ps_tss))

metadata$Soil_type <- factor(metadata$Soil_type, levels = unique(metadata$Soil_type))
metadata$Soil_type

str(ps_tss@sam_data)


# Graph for bray curtis after clustering 

ps_nb <- subset_samples(physeq = ps_tss, Treatment == 'No_Biochar')
ps_nb

ps.nb <- transform_sample_counts(ps_nb,function(otu)otu/ sum(otu))
ord.PCoA.bray <- ordinate(ps.nb,method = "PCoA",distance="bray")

plot_ordination(ps.nb, ord.PCoA.bray, color = "Soil_type", shape="Treatment",title="No_biochar") + geom_point(size =2.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

unique_groups <- unique(sample_data(ps_combined)$Groups)
dark_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#1b4965", "#9a031e")

# Ensure we have enough colors for all groups; repeat colors if necessary
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)
# Step 1: Transform sample counts for the entire phyloseq object (without subsetting)
ps_combined <- transform_sample_counts(ps_tss, function(otu) otu / sum(otu))

library(ggplot2)
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

# Step 1: Perform clustering on Bray-Curtis distance matrix
distance_matrix <- distance(ps_combined, method = "bray")  # Use Bray-Curtis distance
hc <- hclust(as.dist(distance_matrix), method = "ward.D2")  # Hierarchical clustering
num_clusters <- 3  # Specify the number of clusters
clusters <- cutree(hc, k = num_clusters)  # Assign cluster numbers

# Step 2: Add clusters to the sample data
sample_data(ps_combined)$Cluster <- factor(clusters[rownames(sample_data(ps_combined))])

# Step 3: Define color palettes
cluster_colors <- c("#1b9e77", "#d95f02", "#7570b3")  # Adjust as needed
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)

# Step 4: Create PCoA ordination using Bray-Curtis distance
ord.PCoA.bray <- ordinate(ps_combined, method = "PCoA", distance = "bray")

# Step 5: Plot with cluster and group information
p <- plot_ordination(ps_combined, ord.PCoA.bray, color = "Cluster", shape = "Treatment", 
                     title = "PCoA Ordination with Clustering (Bray-Curtis)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # Apply cluster colors
  scale_color_manual(values = cluster_colors) +
  # Add group colors (optional for groups)
  scale_fill_manual(values = group_colors) +
  # Add sample labels
  geom_text_repel(aes(label = Plant_name), size = 3, max.overlaps = 20)

# Display the plot
print(p)

# graph for PcoA ordination with clustering for Jaccard
# Subset samples for Treatment == 'No_biochar'
ps_nb <- subset_samples(physeq = ps_tss, Treatment == 'No_Biochar')
ps_nb

# Step 1: Transform sample counts for subsetted data
ps.nb <- transform_sample_counts(ps_nb, function(otu) otu / sum(otu))

# Step 2: Perform PCoA ordination using Jaccard distance for 'No_biochar' samples
ord.PCoA.jaccard <- ordinate(ps.nb, method = "PCoA", distance = "jaccard")

# Plot PCoA for 'No_biochar' samples
plot_ordination(ps.nb, ord.PCoA.jaccard, color = "Soil_type", shape = "Treatment", title = "No_Biochar (PCoA - Jaccard)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Step 3: Transform the entire phyloseq object for combined analysis
ps_combined <- transform_sample_counts(ps, function(otu) otu / sum(otu))

# Perform clustering on Jaccard distance matrix
distance_matrix_jaccard <- distance(ps_combined, method = "jaccard")  # Use Jaccard distance
hc_jaccard <- hclust(as.dist(distance_matrix_jaccard), method = "ward.D2")  # Hierarchical clustering
num_clusters <- 3  # Specify number of clusters
clusters_jaccard <- cutree(hc_jaccard, k = num_clusters)  # Assign clusters

# Add clusters to sample data
sample_data(ps_combined)$Cluster <- factor(clusters_jaccard[rownames(sample_data(ps_combined))])

# Define color palettes for clusters and groups
cluster_colors <- c("#1b9e77", "#d95f02", "#7570b3")  # Adjust as needed
unique_groups <- unique(sample_data(ps_combined)$Groups)
dark_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#1b4965", "#9a031e")
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)

# Create PCoA ordination using Jaccard distance
ord.PCoA.jaccard_combined <- ordinate(ps_combined, method = "PCoA", distance = "jaccard")

# Plot PCoA with cluster and group information
p <- plot_ordination(ps_combined, ord.PCoA.jaccard_combined, color = "Cluster", shape = "Treatment", 
                     title = "PCoA Ordination with Clustering (Jaccard)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # Apply cluster colors
  scale_color_manual(values = cluster_colors) +
  # Add group colors (optional for groups)
  scale_fill_manual(values = group_colors) +
  # Add sample labels
  geom_text_repel(aes(label = SampleName), size = 3, max.overlaps = 20)

# Display the plot
print(p)

# NMDS Jaccard with clustering plot
# Subset samples for Treatment == 'No_biochar'
ps_nb <- subset_samples(physeq = ps, Treatment == 'No_biochar')
ps_nb

# Step 1: Transform sample counts for subsetted data
ps.nb <- transform_sample_counts(ps_nb, function(otu) otu / sum(otu))

# Step 2: Perform NMDS ordination using Jaccard distance for 'No_biochar' samples
ord.NMDS.jaccard <- ordinate(ps.nb, method = "NMDS", distance = "jaccard")

# Plot NMDS for 'No_biochar' samples
plot_ordination(ps.nb, ord.NMDS.jaccard, color = "Soil", shape = "Treatment", title = "No_biochar (NMDS - Jaccard)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Step 3: Transform the entire phyloseq object for combined analysis
ps_combined <- transform_sample_counts(ps, function(otu) otu / sum(otu))

# Perform clustering on Jaccard distance matrix
distance_matrix_jaccard <- distance(ps_combined, method = "jaccard")  # Use Jaccard distance
hc_jaccard <- hclust(as.dist(distance_matrix_jaccard), method = "ward.D2")  # Hierarchical clustering
num_clusters <- 3  # Specify number of clusters
clusters_jaccard <- cutree(hc_jaccard, k = num_clusters)  # Assign clusters

# Add clusters to sample data
sample_data(ps_combined)$Cluster <- factor(clusters_jaccard[rownames(sample_data(ps_combined))])

# Define color palettes for clusters and groups
cluster_colors <- c("#1b9e77", "#d95f02", "#7570b3")  # Adjust as needed
unique_groups <- unique(sample_data(ps_combined)$Groups)
dark_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#1b4965", "#9a031e")
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)

# Create NMDS ordination using Jaccard distance
ord.NMDS.jaccard_combined <- ordinate(ps_combined, method = "NMDS", distance = "jaccard")

# Plot NMDS with cluster and group information
p <- plot_ordination(ps_combined, ord.NMDS.jaccard_combined, color = "Cluster", shape = "Treatment", 
                     title = "NMDS Ordination with Clustering (Jaccard)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # Apply cluster colors
  scale_color_manual(values = cluster_colors) +
  # Add group colors (optional for groups)
  scale_fill_manual(values = group_colors) +
  # Add sample labels
  geom_text_repel(aes(label = Plant_name), size = 3, max.overlaps = 20)

# Display the plot
print(p)

# nmds ordination bray curtis after clustering 
# Subset samples for Treatment == 'No_biochar'
ps_nb <- subset_samples(physeq = ps, Treatment == 'No_Biochar')
ps_nb

# Step 1: Transform sample counts for subsetted data
ps.nb <- transform_sample_counts(ps_nb, function(otu) otu / sum(otu))

# Step 2: Perform NMDS ordination using Bray-Curtis distance for 'No_biochar' samples
ord.NMDS.bray <- ordinate(ps.nb, method = "NMDS", distance = "bray")

# Plot NMDS for 'No_biochar' samples
plot_ordination(ps.nb, ord.NMDS.bray, color = "Soil_type", shape = "Treatment", title = "No_Biochar (NMDS - Bray-Curtis)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Step 3: Transform the entire phyloseq object for combined analysis
ps_combined <- transform_sample_counts(ps, function(otu) otu / sum(otu))

# Perform clustering on Bray-Curtis distance matrix
distance_matrix_bray <- distance(ps_combined, method = "bray")  # Use Bray-Curtis distance
hc_bray <- hclust(as.dist(distance_matrix_bray), method = "ward.D2")  # Hierarchical clustering
num_clusters <- 3  # Specify number of clusters
clusters_bray <- cutree(hc_bray, k = num_clusters)  # Assign clusters

# Add clusters to sample data
sample_data(ps_combined)$Cluster <- factor(clusters_bray[rownames(sample_data(ps_combined))])

# Define color palettes for clusters and groups
cluster_colors <- c("#1b9e77", "#d95f02", "#7570b3")  # Adjust as needed
unique_groups <- unique(sample_data(ps_combined)$Groups)
dark_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#1b4965", "#9a031e")
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)

# Create NMDS ordination using Bray-Curtis distance
ord.NMDS.bray_combined <- ordinate(ps_combined, method = "NMDS", distance = "bray")

# Plot NMDS with cluster and group information
p <- plot_ordination(ps_combined, ord.NMDS.bray_combined, color = "Cluster", shape = "Treatment", 
                     title = "NMDS Ordination with Clustering (Bray-Curtis)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # Apply cluster colors
  scale_color_manual(values = cluster_colors) +
  # Add group colors (optional for groups)
  scale_fill_manual(values = group_colors) +
  # Add sample labels
  geom_text_repel(aes(label = Plant_name), size = 3, max.overlaps = 20)

# Display the plot
print(p)

# clustring using K means clustring 
 # nmds BC KMC
# Subset samples for Treatment == 'No_biochar'
ps_nb <- subset_samples(physeq = ps, Treatment == 'No_Biochar')
ps_nb

# Step 1: Transform sample counts for subsetted data
ps.nb <- transform_sample_counts(ps_nb, function(otu) otu / sum(otu))

# Step 2: Perform NMDS ordination using Bray-Curtis distance for 'No_biochar' samples
ord.NMDS.bray <- ordinate(ps.nb, method = "NMDS", distance = "bray")

# Plot NMDS for 'No_biochar' samples
plot_ordination(ps.nb, ord.NMDS.bray, color = "Soil_type", shape = "Treatment", title = "No_Biochar (NMDS - Bray-Curtis)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Step 3: Transform the entire phyloseq object for combined analysis
ps_combined <- transform_sample_counts(ps, function(otu) otu / sum(otu))

# Perform distance matrix calculation with Bray-Curtis
distance_matrix_bray <- distance(ps_combined, method = "bray")

# Convert distance matrix to a data frame for K-means clustering
distance_matrix_df <- as.matrix(distance_matrix_bray)
k <- 3  # Number of clusters
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(distance_matrix_df, centers = k)

# Assign K-means cluster numbers to sample data
sample_data(ps_combined)$Cluster <- factor(kmeans_result$cluster[rownames(sample_data(ps_combined))])

# Define color palettes for clusters and groups
cluster_colors <- c("#1b9e77", "#d95f02", "#7570b3")  # Adjust as needed
unique_groups <- unique(sample_data(ps_combined)$Groups)
dark_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#1b4965", "#9a031e")
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)

# Create NMDS ordination using Bray-Curtis distance
ord.NMDS.bray_combined <- ordinate(ps_combined, method = "NMDS", distance = "bray")

# Plot NMDS with K-means cluster and group information
p <- plot_ordination(ps_combined, ord.NMDS.bray_combined, color = "Cluster", shape = "Treatment", 
                     title = "NMDS Ordination with K-means Clustering (Bray-Curtis)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # Apply cluster colors
  scale_color_manual(values = cluster_colors) +
  # Add group colors (optional for groups)
  scale_fill_manual(values = group_colors) +
  # Add sample labels
  geom_text_repel(aes(label = Plant_name), size = 3, max.overlaps = 20)

# Display the plot
print(p)

# nmds jaccard kmc
# Subset samples for Treatment == 'No_biochar'
ps_nb <- subset_samples(physeq = ps, Treatment == 'No_Biochar')
ps_nb

# Step 1: Transform sample counts for subsetted data
ps.nb <- transform_sample_counts(ps_nb, function(otu) otu / sum(otu))

# Step 2: Perform NMDS ordination using Jaccard distance for 'No_biochar' samples
ord.NMDS.jaccard <- ordinate(ps.nb, method = "NMDS", distance = "jaccard")

# Plot NMDS for 'No_biochar' samples
plot_ordination(ps.nb, ord.NMDS.jaccard, color = "Soil_type", shape = "Treatment", title = "No_Biochar (NMDS - Jaccard)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Step 3: Transform the entire phyloseq object for combined analysis
ps_combined <- transform_sample_counts(ps, function(otu) otu / sum(otu))

# Perform distance matrix calculation with Jaccard
distance_matrix_jaccard <- distance(ps_combined, method = "jaccard")

# Convert distance matrix to a data frame for K-means clustering
distance_matrix_df <- as.matrix(distance_matrix_jaccard)
k <- 3  # Number of clusters
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(distance_matrix_df, centers = k)

# Assign K-means cluster numbers to sample data
sample_data(ps_combined)$Cluster <- factor(kmeans_result$cluster[rownames(sample_data(ps_combined))])

# Define color palettes for clusters and groups
cluster_colors <- c("#1b9e77", "#d95f02", "#7570b3")  # Adjust as needed
unique_groups <- unique(sample_data(ps_combined)$Groups)
dark_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#1b4965", "#9a031e")
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)

# Create NMDS ordination using Jaccard distance
ord.NMDS.jaccard_combined <- ordinate(ps_combined, method = "NMDS", distance = "jaccard")

# Plot NMDS with K-means cluster and group information
p <- plot_ordination(ps_combined, ord.NMDS.jaccard_combined, color = "Cluster", shape = "Treatment", 
                     title = "NMDS Ordination with K-means Clustering (Jaccard)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # Apply cluster colors
  scale_color_manual(values = cluster_colors) +
  # Add group colors (optional for groups)
  scale_fill_manual(values = group_colors) +
  # Add sample labels
  geom_text_repel(aes(label = Plant_name), size = 3, max.overlaps = 20)

# Display the plot
print(p)

# PCoA JACCARD KMC
# Subset samples for Treatment == 'No_biochar'
ps_nb <- subset_samples(physeq = ps, Treatment == 'No_Biochar')
ps_nb

# Step 1: Transform sample counts for subsetted data
ps.nb <- transform_sample_counts(ps_nb, function(otu) otu / sum(otu))

# Step 2: Perform PCoA ordination using Jaccard distance for 'No_biochar' samples
ord.PCoA.jaccard <- ordinate(ps.nb, method = "PCoA", distance = "jaccard")

# Plot PCoA for 'No_biochar' samples
plot_ordination(ps.nb, ord.PCoA.jaccard, color = "Soil", shape = "Treatment", title = "No_biochar (PCoA - Jaccard)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Step 3: Transform the entire phyloseq object for combined analysis
ps_combined <- transform_sample_counts(ps, function(otu) otu / sum(otu))

# Perform distance matrix calculation with Jaccard
distance_matrix_jaccard <- distance(ps_combined, method = "jaccard")

# Convert distance matrix to a data frame for K-means clustering
distance_matrix_df <- as.matrix(distance_matrix_jaccard)
k <- 3  # Number of clusters
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(distance_matrix_df, centers = k)

# Assign K-means cluster numbers to sample data
sample_data(ps_combined)$Cluster <- factor(kmeans_result$cluster[rownames(sample_data(ps_combined))])

# Define color palettes for clusters and groups
cluster_colors <- c("#1b9e77", "#d95f02", "#7570b3")  # Adjust as needed
unique_groups <- unique(sample_data(ps_combined)$Groups)
dark_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#1b4965", "#9a031e")
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)

# Create PCoA ordination using Jaccard distance
ord.PCoA.jaccard_combined <- ordinate(ps_combined, method = "PCoA", distance = "jaccard")

# Plot PCoA with K-means cluster and group information
p <- plot_ordination(ps_combined, ord.PCoA.jaccard_combined, color = "Cluster", shape = "Treatment", 
                     title = "PCoA Ordination with K-means Clustering (Jaccard)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # Apply cluster colors
  scale_color_manual(values = cluster_colors) +
  # Add group colors (optional for groups)
  scale_fill_manual(values = group_colors) +
  # Add sample labels
  geom_text_repel(aes(label = Plant_name), size = 3, max.overlaps = 20)

# Display the plot
print(p)

# PCoA BC KMC
# Subset samples for Treatment == 'No_biochar'
ps_nb <- subset_samples(physeq = ps, Treatment == 'No_Biochar')
ps_nb

# Step 1: Transform sample counts for subsetted data
ps.nb <- transform_sample_counts(ps_nb, function(otu) otu / sum(otu))

# Step 2: Perform PCoA ordination using Bray-Curtis distance for 'No_biochar' samples
ord.PCoA.bray <- ordinate(ps.nb, method = "PCoA", distance = "bray")

# Plot PCoA for 'No_biochar' samples
plot_ordination(ps.nb, ord.PCoA.bray, color = "Soil_type", shape = "Treatment", title = "No_biochar (PCoA - Bray-Curtis)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Step 3: Transform the entire phyloseq object for combined analysis
ps_combined <- transform_sample_counts(ps, function(otu) otu / sum(otu))

# Perform distance matrix calculation with Bray-Curtis
distance_matrix_bray <- distance(ps_combined, method = "bray")

# Convert distance matrix to a data frame for K-means clustering
distance_matrix_df <- as.matrix(distance_matrix_bray)
k <- 3  # Number of clusters
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(distance_matrix_df, centers = k)

# Assign K-means cluster numbers to sample data
sample_data(ps_combined)$Cluster <- factor(kmeans_result$cluster[rownames(sample_data(ps_combined))])

# Define color palettes for clusters and groups
cluster_colors <- c("#1b9e77", "#d95f02", "#7570b3")  # Adjust as needed
unique_groups <- unique(sample_data(ps_combined)$Groups)
dark_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#1b4965", "#9a031e")
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)

# Create PCoA ordination using Bray-Curtis distance
ord.PCoA.bray_combined <- ordinate(ps_combined, method = "PCoA", distance = "bray")

# Plot PCoA with K-means cluster and group information
p <- plot_ordination(ps_combined, ord.PCoA.bray_combined, color = "Cluster", shape = "Treatment", 
                     title = "PCoA Ordination with K-means Clustering (Bray-Curtis)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # Apply cluster colors
  scale_color_manual(values = cluster_colors) +
  # Add group colors (optional for groups)
  scale_fill_manual(values = group_colors) +
  # Add sample labels
  geom_text_repel(aes(label = Plant_name), size = 3, max.overlaps = 20)

# Display the plot
print(p)

#For Group Division
metadata$Groups <- factor(metadata$Groups, levels = unique(metadata$Groups))
metadata$Groups

str(ps@sam_data)


# Graph for bray curtis after clustering 

ps_nb <- subset_samples(physeq = ps, Subsetting == 'Subset1')
ps_nb

ps.nb <- transform_sample_counts(ps_nb,function(otu)otu/ sum(otu))
ord.PCoA.bray <- ordinate(ps.nb,method = "PCoA",distance="bray")

plot_ordination(ps.nb, ord.PCoA.bray, color = "Subsetting", shape="Groups",title="Grasses") + geom_point(size =2.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

unique_groups <- unique(sample_data(ps_combined)$Groups)
dark_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                 "#e6ab02", "#a6761d", "#666666", "#1b4965", "#9a031e")

# Ensure we have enough colors for all groups; repeat colors if necessary
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)
# Step 1: Transform sample counts for the entire phyloseq object (without subsetting)
ps_combined <- transform_sample_counts(ps, function(otu) otu / sum(otu))

library(ggplot2)
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

# Step 1: Perform clustering on Bray-Curtis distance matrix
distance_matrix <- distance(ps_combined, method = "bray")  # Use Bray-Curtis distance
hc <- hclust(as.dist(distance_matrix), method = "ward.D2")  # Hierarchical clustering
num_clusters <- 3  # Specify the number of clusters
clusters <- cutree(hc, k = num_clusters)  # Assign cluster numbers

# Step 2: Add clusters to the sample data
sample_data(ps_combined)$Cluster <- factor(clusters[rownames(sample_data(ps_combined))])

# Step 3: Define color palettes
cluster_colors <- c("#1b9e77", "#d95f02", "#7570b3")  # Adjust as needed
group_colors <- setNames(rep(dark_colors, length.out = length(unique_groups)), unique_groups)

# Step 4: Create PCoA ordination using Bray-Curtis distance
ord.PCoA.bray <- ordinate(ps_combined, method = "PCoA", distance = "bray")

# Step 5: Plot with cluster and group information
p <- plot_ordination(ps_combined, ord.PCoA.bray, color = "Groups", shape = "Subsetting", 
                     title = "PCoA Ordination with Clustering (Bray-Curtis)") +
  geom_point(size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # Apply cluster colors
  scale_color_manual(values = c("darkred","blue","darkgreen","violet", "red","orange", "#66a61e", "deeppink", "coral2", "deepskyblue4", "firebrick3")) +
  # Add group colors (optional for groups)
  scale_fill_manual(values = group_colors) +
  # Add sample labels
  geom_text_repel(aes(label = Plant_name), size = 3, max.overlaps = 20)

# Display the plot
print(p)








