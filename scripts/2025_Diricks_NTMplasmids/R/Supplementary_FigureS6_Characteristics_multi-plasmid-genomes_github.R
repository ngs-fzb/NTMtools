
library(ggplot2)
library(patchwork)
library(scales)  # For scales functions
library(dplyr)
library(tidyr)
library(reshape2)  # For melt function
#install.packages("gridExtra")
library(gridExtra)
library(cowplot)  # Load cowplot for adding labels to subplots
library(RColorBrewer)

#Not same cluster means all comparisons between unclustered plasmids, plasmids that belong to different clusters and plasmids that are unclustered vs clustered

### A) Multi-isolate genome vs Size ###

# Read the data from the text file
data <- read.table("Plasmid_Lengths_mash.txt", header = TRUE, sep = "\t")
View(data)

# Filter data for multi-plasmid isolates
filtered_data <- subset(data, multi_plasmid_isolate == "Yes")

# Convert organism_ID to a factor with levels sorted numerically
filtered_data <- filtered_data %>%
  mutate(organism_ID = factor(organism_ID, levels = sort(unique(organism_ID))))

# Calculate counts per multi-plasmid isolate
counts <- filtered_data %>%
  group_by(organism_ID) %>%
  summarise(count = n())

# Choose a color palette based on the number of unique Plasmid_cluster values
n_clusters <- length(unique(filtered_data$Plasmid_cluster))
color_palette <- brewer.pal(n = min(n_clusters, 12), name = "Set3")  # 'Set3' can have up to 12 colors

max_y <- max(filtered_data$Size / 1000, na.rm = TRUE)  # Size divided by 1000 to convert to kbp

# Create a summary for cluster labels with the maximum size
cluster_labels_top <- filtered_data %>%
  group_by(organism_ID) %>%
  summarise(
    label_text = paste(Plasmid_cluster, collapse = "\n"),  # Stack cluster labels
    count_clusters = n()  # Count the number of clusters
  )

# Define fixed y-axis positions for the first cluster and subsequent clusters
base_y_position <- 700  # Fixed y-value for the first cluster
spacing <- 10  # Spacing between cluster labels
cluster_labels_top <- cluster_labels_top %>%
  mutate(y_position = base_y_position - (count_clusters - 1) * spacing)  # Calculate y positions


# Create the plot with cluster labels at the top, stacked vertically
plot1 <- ggplot(filtered_data, aes(x = organism_ID, y = Size / 1000)) +  
  #geom_point(color = "darkgrey", size = 1) +  # Plot individual points in dark grey
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(color = "darkgrey", size = 1, width = 0.2, height = 0) +  # Jitter to reduce overlap
    #geom_violin(fill = "lightgrey", color = "black", adjust = 1.5) +  # Violin plot with light grey fill; very bad
  geom_hline(yintercept = 50, color = "red", linetype = "dashed") +  # Horizontal reference line at 50
  geom_text(data = counts, aes(x = organism_ID, y = max_y + 100, label = count), size = 3, vjust = -0.5) +  # Position counts above max_y
  geom_text(data = cluster_labels_top, aes(x = organism_ID, y = y_position, label = label_text),  # Use fixed y-position based on count
            size = 3, vjust = 0, color = "black", lineheight = 0.85) +  # Adjust vertical position
  #labs(title = "Amount of plasmids per genome\n\n plasmid cluster (ID)", x = "Multi-plasmid genome (organism ID)", y = "Size (kbp)") +  
  labs(title = "Amount of plasmids per genome\n\n plasmid cluster (ID)", x = "Strain (organism ID)", y = "Size (kbp)") + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, margin = margin(t = 20)),  # Center plot title
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels
    axis.title.x = element_text(margin = margin(t = 10), size = 11),  # Adjust x-axis title position
    axis.title.y = element_text(margin = margin(r = 10), size = 11)  # Adjust y-axis title position
  ) +
  scale_y_continuous(
    limits = c(0, NA),  # Start y-axis at 0, upper limit auto-calculated
    breaks = function(limits) {
      c(0, scales::pretty_breaks(n = 10)(limits), 50)  # Include custom breaks and reference line at 50
    }
  )

# Plot the result
plot1



### B) Mash distance vs multi-plasmid isolate ###
# Read the data from the text file
data <- read.table("mash_long_output.txt", header = TRUE, sep = "\t")
#Mash_long_output file generated in "SupplementaryFigure_PlasmidClusterCharacteristics"
View(data)

# Filter the data
filtered_data <- subset(data, 
                        Multi_plasmid_isolate == "Yes" & 
                          organism_ID1 == organism_ID2 & 
                          Plasmid1 != Plasmid2)

# Ensure organism_ID1 is numeric
filtered_data$organism_ID1 <- as.numeric(filtered_data$organism_ID1)

# Sort organism_ID1
filtered_data <- filtered_data[order(filtered_data$organism_ID1),]

# Create the boxplot
plot2 <- ggplot(filtered_data, aes(x = factor(organism_ID1), y = MashDistance)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +  # Add horizontal line at y = 0.5
  #labs(x = "Multi-plasmid genome (organism ID)", y = "Mash distance") +
  labs(x = "Strain (organism ID)", y = "Mash distance") +
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
plot2

### ANI and alignment fraction does not make sense, because all ANI values below 70% or something like that are not outputted by fastANI).
### Within clusters this was not a problem because there all ANI values were above 70%

#Combined plot function was not working --> save plot 1 and 2 seperately, drag and drop in inkscape together and adjust

