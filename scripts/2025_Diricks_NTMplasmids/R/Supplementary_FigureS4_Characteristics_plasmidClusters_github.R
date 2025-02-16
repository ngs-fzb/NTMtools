
library(ggplot2)
library(patchwork)
library(scales)  # For scales functions
library(dplyr)
library(tidyr)
library(reshape2)  # For melt function
#install.packages("gridExtra")
library(gridExtra)
library(cowplot)  # Load cowplot for adding labels to subplots

#Not same cluster means all comparisons between unclustered plasmids, plasmids that belong to different clusters and plasmids that are unclustered vs clustered

### A) Plasmid Size vs Mash cluster ###

# Read the data from the text file
data <- read.table("Plasmid_Lengths_mash.txt", header = TRUE, sep = "\t")
#View(data)

# Convert Plasmid_cluster to a factor with levels sorted numerically and "Unclustered" last
data <- data %>%
  mutate(Plasmid_cluster = factor(Plasmid_cluster, 
                                  levels = c(sort(as.numeric(levels(factor(data$Plasmid_cluster))[levels(factor(data$Plasmid_cluster)) != "U"]), na.last = TRUE), "U")))

# Calculate counts per Plasmid_cluster
counts <- data %>%
  group_by(Plasmid_cluster) %>%
  summarise(count = n())


# Plotting Plasmid Cluster vs Size with smaller dots and counts
plot1 <- ggplot(data, aes(x = Plasmid_cluster, y = Size / 1000)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(color = "darkgrey", size = 1, width = 0.2, height = 0) +  # Jitter to reduce overlap
  geom_hline(yintercept = 50, color = "red", linetype = "dashed") +
  geom_text(data = counts, aes(x = Plasmid_cluster, y = max(data$Size) / 1000 + 50, label = count), size = 3, vjust = 0) +  # Add text labels for counts
  labs(title = "Amount of plasmids per Cluster", x = "Plasmid Cluster (ID)", y = "Plasmid size (kbp)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, margin = margin(t = 20)),  # Center the plot title
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_text(margin = margin(t = 10), size = 11),  # Move the x-axis title down
    axis.title.y = element_text(margin = margin(r = 10), size = 11)  # Move the y-axis title left
  ) +
scale_y_continuous(
  limits = c(0, NA),  # Ensure y-axis starts at 0, with no upper limit defined (auto-calculated)
  breaks = function(limits) {
    c(0, scales::pretty_breaks(n = 10)(limits), 50)  # Ensure 0 and 0.05 are included
  }
)
plot1

### B) Mash_distance vs mash cluster ###

### Only do once to create correct file ###
File <- "Mash_distance"
matrix <- read.table (File, sep="\t", header=TRUE)
View(matrix)
rn <- matrix[,1]
matrix_g <- matrix [,-1]
rownames(matrix_g) <- rn
matrix_numeric <- as.matrix(matrix_g)
View(matrix_numeric)
# Convert the matrix into a long format with 3 columns
mash_long <- melt(as.matrix(matrix_numeric))

# Rename the columns for better clarity
colnames(mash_long) <- c("Plasmid1", "Plasmid2", "MashDistance")

# Remove rows where Plasmid1 == Plasmid2 (since these are distances of a plasmid to itself)
mash_long <- mash_long %>% filter(Plasmid1 != Plasmid2)

# View the result
#View(mash_long)

write.table(mash_long, file = "mash_long_output.txt", sep = "\t", row.names = FALSE, quote = FALSE)

data <- read.table("mash_long_output.txt", header = TRUE, sep = "\t")

# Filter the data for same cluster
same_cluster <- data %>%
  filter(plasmid_clust1 == plasmid_clust2 & plasmid_clust1 != "U")
#View(same_cluster)

# Filter the data for not same cluster
not_same_cluster <- data %>%
  filter(plasmid_clust1 != plasmid_clust2 | 
           (plasmid_clust1 == 'U' & plasmid_clust2 == 'U'))
#View(not_same_cluster)

# Add a new column to indicate 'Not same cluster'
not_same_cluster <- not_same_cluster %>%
  mutate(cluster = "O")
#O = Other = not same cluster
#View(not_same_cluster)

# Combine the data for plotting
plot_data <- same_cluster %>%
  mutate(cluster = as.character(plasmid_clust1)) %>%
  select(cluster, MashDistance) %>%
  bind_rows(not_same_cluster %>% select(cluster, MashDistance))

# Convert 'cluster' to numeric where possible and sort
numeric_clusters <- plot_data %>%
  filter(cluster != "O") %>%
  mutate(cluster = as.numeric(cluster)) %>%
  arrange(cluster) %>%
  mutate(cluster = as.character(cluster))

# Add 'Not same cluster' back
not_same_cluster <- plot_data %>%
  filter(cluster == "O")

# Combine numeric clusters and 'Not same cluster'
plot_data <- bind_rows(numeric_clusters, not_same_cluster)
#plot_data <- bind_rows(numeric_clusters)

# Set factor levels for 'cluster' to maintain the order
plot_data$cluster <- factor(plot_data$cluster, levels = unique(plot_data$cluster))

# Create the whisker plot
plot2 <- ggplot(plot_data, aes(x = cluster, y = MashDistance)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +  # Add horizontal line at y = 95%
  labs(x = "Plasmid Cluster (ID)", y = "mash distance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_text(margin = margin(t = 10), size = 11),  # Move the x-axis title down
    axis.title.y = element_text(margin = margin(r = 10), size = 11)  # Move the y-axis title left
  ) +
  scale_y_continuous(
    limits = c(0, NA),  # Ensure y-axis starts at 0, with no upper limit defined (auto-calculated)
    breaks = function(limits) {
      c(0, scales::pretty_breaks(n = 10)(limits), 0.05)  # Ensure 0 and 0.05 are included
    }
  )
  

### C) ANI vs Mash cluster###
data <- read.table("fastANI_196plasmids_fragLen500.txt", header = TRUE, sep = "\t")

# Filter the data for same cluster but different contig
same_cluster_diff_contig <- data %>%
  filter(plasmid_clust1 == plasmid_clust2 & Self_comparison == "Different_Contig" & plasmid_clust1 != "U")
#View(same_cluster_diff_contig)

# Filter the data for not same cluster
not_same_cluster <- data %>%
  filter(plasmid_clust1 != plasmid_clust2 | 
           (plasmid_clust1 == 'U' & plasmid_clust2 == 'U' & Self_comparison == "Different_Contig"))
#View(not_same_cluster)

# Add a new column to indicate 'Not same cluster'
not_same_cluster <- not_same_cluster %>%
  mutate(cluster = "Not same cluster")
#View(not_same_cluster)

# Combine the data for plotting
plot_data <- same_cluster_diff_contig %>%
  mutate(cluster = as.character(plasmid_clust1)) %>%
  select(cluster, ANI) %>%
  bind_rows(not_same_cluster %>% select(cluster, ANI))

# Convert 'cluster' to numeric where possible and sort
numeric_clusters <- plot_data %>%
  filter(cluster != "Not same cluster") %>%
  mutate(cluster = as.numeric(cluster)) %>%
  arrange(cluster) %>%
  mutate(cluster = as.character(cluster))

# Add 'Not same cluster' back
#not_same_cluster <- plot_data %>%
#  filter(cluster == "Not same cluster")

# Combine numeric clusters and 'Not same cluster'
#plot_data <- bind_rows(numeric_clusters, not_same_cluster)
plot_data <- bind_rows(numeric_clusters)

# Set factor levels for 'cluster' to maintain the order
plot_data$cluster <- factor(plot_data$cluster, levels = unique(plot_data$cluster))

# Create the whisker plot
# Create the whisker plot
plot3 <- ggplot(plot_data, aes(x = cluster, y = ANI)) +
  geom_boxplot() +
  #geom_hline(yintercept = 95, color = "red", linetype = "dashed") +  # Add horizontal line at y = 95%
  labs(x = "Plasmid Cluster (ID)", y = "ANI (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_text(margin = margin(t = 10), size = 11),  # Move the x-axis title down
    axis.title.y = element_text(margin = margin(r = 10), size = 11)  # Move the y-axis title left
  ) +
  scale_y_continuous(
    limits = c(0, NA),  # Ensure y-axis starts at 0, with no upper limit defined (auto-calculated)
    breaks = function(limits) {
      c(0, scales::pretty_breaks(n = 10)(limits))  # Ensure 0 and 0.05 are included
    }
  )

### D) ANI alignment fraction vs Mash cluster###
data <- read.table("fastANI_196plasmids_fragLen500.txt", header = TRUE, sep = "\t")

# Filter the data for same cluster but different contig
same_cluster_diff_contig <- data %>%
  filter(plasmid_clust1 == plasmid_clust2 & Self_comparison == "Different_Contig" & plasmid_clust1 != "U")
#View(same_cluster_diff_contig)

# Filter the data for not same cluster
not_same_cluster <- data %>%
  filter(plasmid_clust1 != plasmid_clust2 | 
           (plasmid_clust1 == 'U' & plasmid_clust2 == 'U' & Self_comparison == "Different_Contig"))
#View(not_same_cluster)

# Add a new column to indicate 'Not same cluster'
not_same_cluster <- not_same_cluster %>%
  mutate(cluster = "Not same cluster")
#View(not_same_cluster)

# Combine the data for plotting
plot_data <- same_cluster_diff_contig %>%
  mutate(cluster = as.character(plasmid_clust1)) %>%
  select(cluster, Alignment_fraction_perc) %>%
  bind_rows(not_same_cluster %>% select(cluster, Alignment_fraction_perc))

# Convert 'cluster' to numeric where possible and sort
numeric_clusters <- plot_data %>%
  filter(cluster != "Not same cluster") %>%
  mutate(cluster = as.numeric(cluster)) %>%
  arrange(cluster) %>%
  mutate(cluster = as.character(cluster))

# Add 'Not same cluster' back
#not_same_cluster <- plot_data %>%
#  filter(cluster == "Not same cluster")

# Combine numeric clusters and 'Not same cluster'
#plot_data <- bind_rows(numeric_clusters, not_same_cluster)
plot_data <- bind_rows(numeric_clusters)

# Set factor levels for 'cluster' to maintain the order
plot_data$cluster <- factor(plot_data$cluster, levels = unique(plot_data$cluster))

# Create the whisker plot
# Create the whisker plot
plot4 <- ggplot(plot_data, aes(x = cluster, y = Alignment_fraction_perc)) +
  geom_boxplot() +
  #geom_hline(yintercept = 95, color = "red", linetype = "dashed") +  # Add horizontal line at y = 95%
  labs(x = "Plasmid Cluster (ID)", y = "Alignment fraction (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_text(margin = margin(t = 10), size = 11),  # Move the x-axis title down
    axis.title.y = element_text(margin = margin(r = 10), size = 11)  # Move the y-axis title left
  ) +
  scale_y_continuous(
    limits = c(0, NA),  # Ensure y-axis starts at 0, with no upper limit defined (auto-calculated)
    breaks = function(limits) {
      c(0, scales::pretty_breaks(n = 10)(limits))  # Ensure 0 and 0.05 are included
    }
  )

# Create the labeled plots
plot1_labeled <- ggdraw() + draw_plot(plot1) + draw_label("A", x = 0.03, y = 0.95, size = 10)
plot2_labeled <- ggdraw() + draw_plot(plot2) + draw_label("B", x = 0.03, y = 0.95, size = 10)
plot3_labeled <- ggdraw() + draw_plot(plot3) + draw_label("C", x = 0.03, y = 0.95, size = 10)
plot4_labeled <- ggdraw() + draw_plot(plot4) + draw_label("D", x = 0.03, y = 0.95, size = 10)

# Use gridExtra to combine labeled plots
combined_plot <- grid.arrange(plot1_labeled, plot2_labeled, plot3_labeled, plot4_labeled, ncol = 2)  # Specify the number of columns

# Print the combined plot
print(combined_plot)

