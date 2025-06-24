library("phyloseq")
library("ggplot2")
library("vegan")
library("FNN")
library("igraph")
library("RColorBrewer")
library("patchwork") 
library("gridExtra")
library("ggrepel")

# Step 1: Read the metadata
metadata <- read.csv("C:/Users/Debopriya/Downloads/beta_diversity_NIWASM/clusterprofile/metadata_Mitali.csv", sep = ",", header = TRUE, row.names = 1)

# Step 2: Load phyloseq object
ps <- readRDS("C:/Users/Debopriya/Downloads/beta_diversity_NIWASM/clusterprofile/ps.rds")
ps_tss <- transform_sample_counts(ps, function(x) x / sum(x))
sample_names(ps_tss) <- rownames(metadata)
sample_data(ps_tss) <- sample_data(metadata)
ps_tss@sam_data$Groups <- factor(ps_tss@sam_data$Groups, levels = unique(ps_tss@sam_data$Groups))
ps_tss@sam_data$Groups
str(sample_data(ps_tss))

metadata$Soil <- factor(metadata$Soil, levels = unique(metadata$Soil))
metadata$Soil

str(ps_tss@sam_data)

# Step 2: PCoA using Bray-Curtis
ord_pcoa <- ordinate(ps_tss, method = "PCoA", distance = "bray")
coords <- as.data.frame(ord_pcoa$vectors)
coords$SampleID <- rownames(coords)


# Step 3: KNN + Louvain clustering
k <- 2  # Number of neighbors
knn_result <- get.knn(coords[, 1:2], k = k)  # Use first 2 axes for clustering
edges <- do.call(rbind, lapply(1:nrow(knn_result$nn.index), function(i) {
  cbind(i, knn_result$nn.index[i, ])
}))
g <- graph_from_edgelist(edges, directed = FALSE)
clusters <- cluster_louvain(g)

# Add clustering and treatment info to coords
coords$Cluster <- factor(clusters$membership)
coords$Treatment <- sample_data(ps_tss)$Treatment

# Step 4: Add sample labels from sample_data(ps_tss)
coords$Plant_name <- sample_data(ps_tss)$Plant_name  # <- Make sure this column exists

# Step 5: Plot with labels
ggplot(coords, aes(x = Axis.1, y = Axis.2, color = Cluster, shape = Treatment, label = Plant_name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(title = "PCoA Ordination with KNN Clustering (Bray-Curtis)",
       x = "PCoA 1", y = "PCoA 2") +
  scale_color_brewer(palette = "Set1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# Step 2: PCoA using Bray-Curtis
ord_pcoa <- ordinate(ps_tss, method = "PCoA", distance = "jaccard")
coords <- as.data.frame(ord_pcoa$vectors)
coords$SampleID <- rownames(coords)


# Step 3: KNN + Louvain clustering
k <- 2  # Number of neighbors
knn_result <- get.knn(coords[, 1:2], k = k)  # Use first 2 axes for clustering
edges <- do.call(rbind, lapply(1:nrow(knn_result$nn.index), function(i) {
  cbind(i, knn_result$nn.index[i, ])
}))
g <- graph_from_edgelist(edges, directed = FALSE)
clusters <- cluster_louvain(g)

# Add clustering and treatment info to coords
coords$Cluster <- factor(clusters$membership)
coords$Treatment <- sample_data(ps_tss)$Treatment

# Step 4: Add sample labels from sample_data(ps_tss)
coords$Plant_name <- sample_data(ps_tss)$Plant_name  # <- Make sure this column exists

# Step 5: Plot with labels
ggplot(coords, aes(x = Axis.1, y = Axis.2, color = Cluster, shape = Treatment, label = Plant_name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(title = "PCoA Ordination with KNN Clustering (Jaccard)",
       x = "PCoA 1", y = "PCoA 2") +
  scale_color_brewer(palette = "Set1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Step 2: NMDS using Bray-Curtis
ord_nmds <- ordinate(ps_tss, method = "NMDS", distance = "bray", trymax = 100)

# Extract sample coordinates
coords <- as.data.frame(scores(ord_nmds, display = "sites"))
coords$SampleID <- rownames(coords)

# Step 3: KNN + Louvain clustering
k <- 2  # Adjust the number of neighbors if needed
knn_result <- get.knn(coords[, 1:2], k = k)
edges <- do.call(rbind, lapply(1:nrow(knn_result$nn.index), function(i) {
  cbind(i, knn_result$nn.index[i, ])
}))
g <- graph_from_edgelist(edges, directed = FALSE)
clusters <- cluster_louvain(g)

# Step 4: Add sample metadata (Treatment & Plant_name)
coords$Cluster <- factor(clusters$membership)
coords$Treatment <- sample_data(ps_tss)$Treatment
coords$Plant_name <- sample_data(ps_tss)$Plant_name  # Ensure this column exists

# Step 5: Plot with clusters and labels
ggplot(coords, aes(x = NMDS1, y = NMDS2, color = Cluster, shape = Treatment, label = Plant_name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(title = "NMDS Ordination with KNN Clustering (Bray-Curtis)",
       x = "NMDS1", y = "NMDS2") +
  scale_color_brewer(palette = "Set1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Step 2: NMDS using jaccard
ord_nmds <- ordinate(ps_tss, method = "NMDS", distance = "jaccard", trymax = 100)

# Extract sample coordinates
coords <- as.data.frame(scores(ord_nmds, display = "sites"))
coords$SampleID <- rownames(coords)

# Step 3: KNN + Louvain clustering
k <- 2  # Adjust the number of neighbors if needed
knn_result <- get.knn(coords[, 1:2], k = k)
edges <- do.call(rbind, lapply(1:nrow(knn_result$nn.index), function(i) {
  cbind(i, knn_result$nn.index[i, ])
}))
g <- graph_from_edgelist(edges, directed = FALSE)
clusters <- cluster_louvain(g)

# Step 4: Add sample metadata (Treatment & Plant_name)
coords$Cluster <- factor(clusters$membership)
coords$Treatment <- sample_data(ps_tss)$Treatment
coords$Plant_name <- sample_data(ps_tss)$Plant_name  # Ensure this column exists

# Step 5: Plot with clusters and labels
ggplot(coords, aes(x = NMDS1, y = NMDS2, color = Cluster, shape = Treatment, label = Plant_name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(title = "NMDS Ordination with KNN Clustering (jaccard)",
       x = "NMDS1", y = "NMDS2") +
  scale_color_brewer(palette = "Set1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
