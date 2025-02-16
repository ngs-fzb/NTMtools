#Written by Margo Diricks, with help of chatGPT.


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)
library(grid)

# Define the column classes, where "character" is set for the column with the problematic value
col_classes <- c("character", "character", "numeric", "character", "character", "NULL", "NULL", "NULL")

# Read the data from the text file
data <- read.table("AMRfinder_results_mash_WithBra100.txt", header = TRUE, sep = "\t", colClasses = col_classes, quote = "\"")
#This table contains all results that are not "internal_stop" and that have >30% identity and >70% coverage
View(data)
# The plasmid names are preceeded by the cluster number and Gene by the type (S = stress, A = AMR, V = virulence)
#QAC = Quaternary ammonia compounds
# Ensure that columns are named correctly
# Assuming the columns are: Plasmid, Gene, Identity, Type, Group
colnames(data) <- c("Plasmid", "Gene", "Identity", "Type", "Class")

# Keep the original order of plasmids
data$Plasmid <- factor(data$Plasmid, levels = unique(data$Plasmid))




###Version 1###
#Add categories 
#data$category <- gsub("\\-.*","",data$Gene)
#data$Gene_2 <- gsub(".\\-","",data$Gene)
data$Gene <- with(data, factor(Gene, levels = Gene[order(Class)]))

data$Gene_unique <- paste(data$Gene, data$Class, sep = "_")
data$Gene_unique <- with(data, factor(Gene_unique, levels = unique(Gene_unique[order(Class)])))

#Heatmap with facets and ggplot
ggplot(data)+
  geom_tile(aes(y=Plasmid,x=Gene_unique,fill=Identity))+
  #facet_grid(.~category,drop=TRUE, scales="free", space="free")+
  facet_grid(.~data$Type,drop=TRUE, scales="free", space="free")+
  labs(x = "Gene_Class")+
  scale_fill_gradient(name = "Identity (%)")+  # Change the legend title+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_text(vjust = 5),
        panel.spacing = unit(0.5, "lines")  # Increase space between horizontal rows) 
  )
  
