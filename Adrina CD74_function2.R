# script to perform differential gene expression analysis using DESeq2 package




# Step 1: preparing count data ----------------


performRNASeqAnalysis <- function(experiment, Cell_type, timepoint, condition, Reference) {

# read in counts data
WD <- setwd("//isdsyn1-1.isd.med.uni-muenchen.de/BD-Bernhagen/Simon Ebert/CD74_Adrian")

counts_data <- read.table("counts.txt", header = TRUE, sep="\t")

str(counts_data)


# Check the structure again
str(counts_data)

mode(counts_data)



head(counts_data)


# read in sample info
colData <- project_info<- read.table("meta.txt", header = TRUE, sep="\t")
row.names(colData)<-project_info$sample_id



# Keep samples with cell_type == Cell_type
project_info <- project_info[project_info$cell_type == Cell_type, ]

# Keep samples with stimulation_time == timepoint
project_info <- project_info[project_info$stimulation_time == timepoint, ]

# Keep samples with cytokine_condition == condition
project_info <- project_info[project_info$cytokine_condition %in% condition, ]





# Subset the count matrix based on sample names in project_info
selected_samples <- project_info$sample_id
counts_data <- counts_data[, selected_samples]
rownames(project_info) <- project_info$sample_id

# Verify the updated count matrix
head(counts_data)



colData <- project_info

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = colData,
                       design = ~ cytokine_condition)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]




dds

# set the factor level
dds$Condition <- relevel(dds$cytokine_condition, ref = Reference)



# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)


Cell_type <- gsub("_", "-", Cell_type)

# PCA plot
plotPCA(rlog(dds), intgroup = "cytokine_condition")
ggsave(paste0(experiment, timepoint,Cell_type, "_PCA.svg"), width = 7, height = 8)

# MA plot
plotMA(res)
res
res_data <- as.data.frame(res)

res_data$Ensemble_ID <- rownames(res_data)








summary(res)




# Load the biomaRt library


# Specify the BioMart database
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", path = "/biomart/martservice")

# Specify the dataset you want to use
dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)



# Use biomaRt to get gene symbols
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = res_data$Ensemble_ID,
                      mart = dataset)

# Merge the original dataframe with the gene symbols
result_dataframe <- merge(res_data, gene_symbols, by.x = "Ensemble_ID", by.y = "ensembl_gene_id")
result_dataframe <- result_dataframe %>%
  rename(SYMBOL = external_gene_name)

# Print the result
print(result_dataframe)


write.xlsx(result_dataframe, paste0(experiment,timepoint,Cell_type, "_DEGs.xlsx"))


####







#Plot volcano and GO

df <- result_dataframe
head(df)

gene_column <- "SYMBOL" 

# Remove genes starting with "Gm" or ending with "Rik" because no interest
filtered_df <- df[!(grepl("^Gm", df[[gene_column]]) | grepl("Rik$", df[[gene_column]])), ]

# Now, filtered_df contains the DataFrame with genes removed
filtered_df

#y_limits <- c(0, 10)
#x_limits <- c(-10, 10)

p1<-EnhancedVolcano(filtered_df,
                   lab = filtered_df$SYMBOL,       
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = paste("DeSeq2 results of ", Cell_type, "T-Cells", experiment,"after", timepoint),
                    subtitle= "",
                    pCutoff = 0.05,
                    selectLab = c("CD74", "CXCR4","MIF"), #Lable the genes of interest
                    FCcutoff = 1.5,
                    max.overlaps = 10,
                    cutoffLineType = "blank",
                    axisLabSize = 12,
                    titleLabSize = 12,
                   boxedLabels = TRUE,
                    subtitleLabSize = 12,
                    captionLabSize = 12,
                    legendLabSize = 12,
                    ylab = bquote(~-Log[10] ~ (p.adjust)),
                    #xlim = c(-3, 3),
                    xlab=paste("Log2fold change (", experiment, ")"),
                    pointSize = 2,
                    labSize = 5,
                    caption =  "Log2fold change > |1.5|, p.adjust < 0.05",
                    drawConnectors = TRUE,
                    widthConnectors = 0.5)
#+
#  coord_cartesian(ylim = y_limits)



show(p1)
ggsave(paste0(experiment,timepoint, Cell_type, "_volcano.svg"), width = 9, height = 10)

result_dataframe
result_dataframe <- subset(result_dataframe, padj < 0.05)



# Select genes you want to include in the dot plot
selected_genes <- c("CD74", "CXCR4", "MIF")

# Filter dataframe for selected genes
selected_data <- result_dataframe %>%
  filter(SYMBOL %in% selected_genes)

# Check if selected_data is not empty before creating the dot plot
if (!is_empty(selected_data)) {
  
  # Create a dot plot using ggplot2
  if (nrow(selected_data) == 1) {
    # If only one gene is found, use black color for dot and p value without gradient
    dotplot <- ggplot(selected_data, aes(x = SYMBOL, y = log2FoldChange)) +
      geom_point(aes(color = as.factor(sprintf("%.3g", padj))), size = 10) +  # Use formatted p.adjust for color
      labs(title = paste("Dot Plot of", Cell_type, "T-Cells", "after", timepoint, "in", experiment),
           x = "Genes of interest",
           y = paste("log2fold change (", experiment, ")"),
           color = "p.adjust") +
      theme(legend.title = element_text(size = 12)) +  # Change the legend caption to "p.adjust"
      scale_color_manual(values = "black", guide = "legend") +  # Set color to black
      ylim(c(-4, 3))
    
    # Save the dot plot with width 3 and height 4
    ggsave(paste0(experiment, timepoint, Cell_type, "_dotplot.svg"), plot = dotplot, width = 6, height = 5)
  } else {
    # If more than one gene is found, use gradient for dot colors
    dotplot <- ggplot(selected_data, aes(x = SYMBOL, y = log2FoldChange, color = padj)) +
      geom_point(size = 10) +
      scale_color_gradient(low = "red", high = "blue", guide = "colourbar") +
      labs(title = paste("Dot Plot of", Cell_type, "T-Cells", "after", timepoint, "in", experiment),
           x = "Genes of interest",
           y = paste("log2fold change (", experiment, ")"),
           color = "p.adjust") +
      theme(legend.title = element_text(size = 12)) +  # Change the legend caption to "p.adjust"
      ylim(c(-4, 3))
    
    # Save the dot plot with default width and height
    ggsave(paste0(experiment, timepoint, Cell_type, "_dotplot.svg"), plot = dotplot, width = 6, height = 5)
  }
  
  # Show the dot plot
  print(dotplot)
} else {
  cat("No selected genes found. Skipping dot plot creation.\n")
}


rm(list = ls())
cat("RNA-Seq analysis completed. Workspace cleared.\n")

}


