library(DESeq2)
library(tidyverse)
library(openxlsx)
library(svglite)
library(EnhancedVolcano)
library(ggplot2)
library(ggbreak)
library(svglite)
library(biomaRt)
library(dplyr)

#Read data
WD <- setwd("//isdsyn1-1.isd.med.uni-muenchen.de/BD-Bernhagen/Simon Ebert/CD74_Adrian")
source("Adrina CD74_function2.R")






# Define different conditions
conditions_list <- list(
  list(experiment = "Th0 vs. resting",
       Cell_type = "CD4_Memory",
       timepoint = "5d",
       condition = c("Th0", "Resting"),
       Reference = "Resting"),
  
  list(experiment = "Th0 vs. resting",
       Cell_type = "CD4_Memory",
       timepoint = "16h",
       condition = c("Th0", "Resting"),
       Reference = "Resting"),
  
  list(experiment = "Th0 vs. resting",
       Cell_type = "CD4_Naive",
       timepoint = "16h",
       condition = c("Th0", "Resting"),
       Reference = "Resting"),
  
  list(experiment = "Th0 vs. resting",
       Cell_type = "CD4_Naive",
       timepoint = "5d",
       condition = c("Th0", "Resting"),
       Reference = "Resting"),
  
  list(experiment = "Th1 vs. Th0",
       Cell_type = "CD4_Memory",
       timepoint = "16h",
       condition = c("Th1", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "Th1 vs. Th0",
       Cell_type = "CD4_Memory",
       timepoint = "5d",
       condition = c("Th1", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "Th1 vs. Th0",
       Cell_type = "CD4_Naive",
       timepoint = "16h",
       condition = c("Th1", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "Th1 vs. Th0",
       Cell_type = "CD4_Naive",
       timepoint = "5d",
       condition = c("Th1", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "Th2 vs. Th0",
       Cell_type = "CD4_Naive",
       timepoint = "16h",
       condition = c("Th2", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "Th2 vs. Th0",
       Cell_type = "CD4_Naive",
       timepoint = "5d",
       condition = c("Th2", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "Th2 vs. Th0",
       Cell_type = "CD4_Memory",
       timepoint = "16h",
       condition = c("Th2", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "Th2 vs. Th0",
       Cell_type = "CD4_Memory",
       timepoint = "5d",
       condition = c("Th2", "Th0"),
       Reference = "Th0"),
  
  
  list(experiment = "Th17 vs. Th0",
       Cell_type = "CD4_Memory",
       timepoint = "16h",
       condition = c("Th17", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "Th17 vs. Th0",
       Cell_type = "CD4_Memory",
       timepoint = "5d",
       condition = c("Th17", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "Th17 vs. Th0",
       Cell_type = "CD4_Naive",
       timepoint = "16h",
       condition = c("Th17", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "Th17 vs. Th0",
       Cell_type = "CD4_Naive",
       timepoint = "5d",
       condition = c("Th17", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "iTreg vs. Th0",
       Cell_type = "CD4_Naive",
       timepoint = "16h",
       condition = c("iTreg", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "iTreg vs. Th0",
       Cell_type = "CD4_Naive",
       timepoint = "5d",
       condition = c("iTreg", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "iTreg vs. Th0",
       Cell_type = "CD4_Memory",
       timepoint = "16h",
       condition = c("iTreg", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "iTreg vs. Th0",
       Cell_type = "CD4_Memory",
       timepoint = "5d",
       condition = c("iTreg", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "IFNB vs. Th0",
       Cell_type = "CD4_Memory",
       timepoint = "16h",
       condition = c("IFNB", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "IFNB vs. Th0",
       Cell_type = "CD4_Memory",
       timepoint = "5d",
       condition = c("IFNB", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "IFNB vs. Th0",
       Cell_type = "CD4_Naive",
       timepoint = "16h",
       condition = c("IFNB", "Th0"),
       Reference = "Th0"),
  
  list(experiment = "IFNB vs. Th0",
       Cell_type = "CD4_Naive",
       timepoint = "5d",
       condition = c("IFNB", "Th0"),
       Reference = "Th0")
  
)

# Run the analysis for each condition
for (cond in conditions_list) {
  performRNASeqAnalysis(cond$experiment, cond$Cell_type, cond$timepoint, cond$condition, cond$Reference)
}


library(pander)

# Set the path to the folder containing SVG files
folder_path <- "//isdsyn1-1.isd.med.uni-muenchen.de/BD-Bernhagen/Simon Ebert/CD74_Adrian/"

# List all SVG files in the folder that have "volcano" in their name
svg_files <- list.files(folder_path, pattern = "volcano", full.names = TRUE)

# Check if any files were found
if (length(svg_files) == 0) {
  stop("No SVG files with 'volcano' in their name found.")
}

# Create a character vector with HTML image tags
html_images <- sapply(svg_files, function(file) {
  sprintf('<img src="%s" width="400">', file)
})

# Combine HTML image tags into a single HTML document
html_content <- paste(html_images, collapse = "\n")

# Save the HTML document to a file
html_file <- "combined_volcano_plots.html"
cat(sprintf('<html>\n<body>\n%s\n</body>\n</html>', html_content), file = html_file)

# Open the HTML file in a web browser
browseURL(html_file)

# Set the path to the folder containing SVG files
folder_path <- "//isdsyn1-1.isd.med.uni-muenchen.de/BD-Bernhagen/Simon Ebert/CD74_Adrian/"

# List all SVG files in the folder that have "volcano" in their name
svg_files <- list.files(folder_path, pattern = "dot", full.names = TRUE)

# Check if any files were found
if (length(svg_files) == 0) {
  stop("No SVG files with 'volcano' in their name found.")
}

# Create a character vector with HTML image tags
html_images <- sapply(svg_files, function(file) {
  sprintf('<img src="%s" width="400">', file)
})

# Combine HTML image tags into a single HTML document
html_content <- paste(html_images, collapse = "\n")

# Save the HTML document to a file
html_file <- "combined_dot_plots.html"
cat(sprintf('<html>\n<body>\n%s\n</body>\n</html>', html_content), file = html_file)

# Open the HTML file in a web browser
browseURL(html_file)

# Set the path to the folder containing SVG files
folder_path <- "//isdsyn1-1.isd.med.uni-muenchen.de/BD-Bernhagen/Simon Ebert/CD74_Adrian/"

# List all SVG files in the folder that have "volcano" in their name
svg_files <- list.files(folder_path, pattern = "PCA", full.names = TRUE)

# Check if any files were found
if (length(svg_files) == 0) {
  stop("No SVG files with 'volcano' in their name found.")
}

# Create a character vector with HTML image tags
html_images <- sapply(svg_files, function(file) {
  sprintf('<img src="%s" width="400">', file)
})

# Combine HTML image tags into a single HTML document
html_content <- paste(html_images, collapse = "\n")

# Save the HTML document to a file
html_file <- "combined_PCA_plots.html"
cat(sprintf('<html>\n<body>\n%s\n</body>\n</html>', html_content), file = html_file)

# Open the HTML file in a web browser
browseURL(html_file)

####Done###






