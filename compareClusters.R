library(ggplot2)
library(tidyverse)
library(UpSetR)
library(ComplexUpset)
library(ComplexHeatmap)
library(VennDiagram)
library(pheatmap)
library(RColorBrewer)

# We want overlap between enrichment genes and other trubetskoy enrichments...

csv_data1 <- read.csv("Ora/Enriched/Trubetskoy2022broadcoding_significant_rows_sorted.csv", stringsAsFactors = FALSE)
csv_data2 <- read.csv("Ora/Enriched/GOBPID_significant_rows_sorted.csv", stringsAsFactors = FALSE)
# Combine 'alg' and 'cl' columns into a new column 'algcl'


split_gene_overlap <- function(csv_data){
  csv_data$algcl <- paste(csv_data$alg, csv_data$cl, sep = "")

  gene_df_long <- csv_data %>%
    separate_rows(overlapGenes, sep = ", ") %>%
    filter(!is.na(overlapGenes))


  df <- data.frame(
    gene = gene_df_long$overlapGenes,
    algcl = gene_df_long$algcl
  )

  df <- df[order(df$algcl),]

  return(df)
}

calculate_overlap <- function(set1, set2) {
  length(intersect(set1, set2))
}

same_comparison <- function(df) {
  # Create a matrix to store the overlap between each pair of algcl permutations
  unique_algcl <- unique(df$algcl)
  overlap_matrix <- matrix(NA, nrow = length(unique_algcl), ncol = length(unique_algcl))
  rownames(overlap_matrix) <- colnames(overlap_matrix) <- unique_algcl


  # Fill in the overlap matrix using the calculate_overlap function
  for (i in seq_along(unique_algcl)) {
    for (j in seq_along(unique_algcl)) {
      genes_i <- df$gene[df$algcl == unique_algcl[i]]
      genes_j <- df$gene[df$algcl == unique_algcl[j]]
      if (i == j) {
        overlap_matrix[i, j] <- 0
      } else {
        overlapped <- calculate_overlap(genes_i, genes_j)
        overlap_matrix[i, j] <- overlapped
      }
    }
  }
  return(overlap_matrix)
}

# Function to calculate gene overlap between two algcl values
# calculate_overlap <- function(algcl1, algcl2, df1, df2) {
#   genes1 <- unlist(strsplit(df1$gene[df1$algcl == algcl1], " "))
#   genes2 <- unlist(strsplit(df2$gene[df2$algcl == algcl2], " "))
#   length(intersect(genes1, genes2))
# }

compare_diff <- function(df1, df2){
  # Extract unique algcl values from both dataframes
  unique_algcl1 <- unique(df1$algcl)
  unique_algcl2 <- unique(df2$algcl)

  # Initialize an empty matrix to store the overlap
  overlap_matrix <- matrix(NA, nrow = length(unique_algcl1), ncol = length(unique_algcl2))
  rownames(overlap_matrix) <- unique_algcl1
  colnames(overlap_matrix) <- unique_algcl2

  # Fill in the overlap matrix using the calculate_overlap function
  for (i in seq_along(unique_algcl1)) {
    for (j in seq_along(unique_algcl2)) {
      overlap_matrix[i, j] <- ifelse(i == j, 0, calculate_overlap(unique_algcl1[i], unique_algcl2[j], df1, df2))
    }
  }
  # Identify rows and columns with no overlaps
  rows_to_drop <- which(rowSums(overlap_matrix == 0, na.rm = TRUE) == ncol(overlap_matrix))
  cols_to_drop <- which(colSums(overlap_matrix == 0, na.rm = TRUE) == nrow(overlap_matrix))

  # Drop rows and columns with no overlaps
  overlap_matrix_filtered <- overlap_matrix[-rows_to_drop, -cols_to_drop]

  return(overlap_matrix_filtered)

}


df1 <- split_gene_overlap(csv_data1)
df2 <- split_gene_overlap(csv_data2)
print("Read CSV")
#verlap_matrix <- compare_diff(df1, df2)

overlap_matrix <- same_comparison(df1)

print("Matrix created")
# Create a heatmap using pheatmap with color palette
color_palette <- colorRampPalette(c("white", "blue", "darkblue"))(n = 100)
pheatmap(
  overlap_matrix,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color = color_palette,
  main = "Overlap Heatmap for Trubetskoy and GOPBPID Cluster Comparison",
  legend = TRUE,
  #filename = "UpsetPlot/TrubetskoyGOBPIDComparision.png",
  angle_col = 45,
  width = 1000,  # Adjust width of the image
  height = 600
)
print("Done")


# Function to calculate gene overlap between two algcl values
calculate_overlap_genes <- function(algcl1, algcl2, df1, df2) {
  genes1 <- unlist(strsplit(df1$gene[df1$algcl == algcl1], " "))
  genes2 <- unlist(strsplit(df2$gene[df2$algcl == algcl2], " "))
  intersect_genes <- intersect(genes1, genes2)
  return(ifelse(length(intersect_genes) > 0, paste(intersect_genes, collapse = ", "), NA))
}

# Create a table of overlapped genes
overlapped_genes_table <- matrix(NA, nrow = nrow(overlap_matrix), ncol = ncol(overlap_matrix))
rownames(overlapped_genes_table) <- rownames(overlap_matrix)
colnames(overlapped_genes_table) <- colnames(overlap_matrix)

for (i in seq_along(rownames(overlap_matrix))) {
  for (j in seq_along(colnames(overlap_matrix))) {
    overlap_value <- overlap_matrix[i, j]
    if (overlap_value > 0) {
      algcl1 <- rownames(overlap_matrix)[i]
      algcl2 <- colnames(overlap_matrix)[j]
      overlapped_genes <- calculate_overlap_genes(algcl1, algcl2, df1, df2)
      overlapped_genes_table[i, j] <- overlapped_genes
    } else {
      overlapped_genes_table[i, j] <- NA
    }
  }
}


overlapped_genes_table <- overlapped_genes_table[!apply(is.na(overlapped_genes_table), 1, all), ]
overlapped_genes_table <- overlapped_genes_table[, !apply(is.na(overlapped_genes_table), 2, all)]


# Print the table of overlapped genes
print(overlapped_genes_table)





# Function to calculate gene overlap between two algcl values
calculate_overlap_genes <- function(algcl1, algcl2, df1, df2) {
  genes1 <- unlist(strsplit(df1$gene[df1$algcl == algcl1], " "))
  genes2 <- unlist(strsplit(df2$gene[df2$algcl == algcl2], " "))
  intersect_genes <- intersect(genes1, genes2)
  return(intersect_genes)
}





# Create a data frame to store unique genes with counts
gene_frequency_df <- data.frame(Gene = character(), Frequency = integer(), stringsAsFactors = FALSE)

for (i in seq_along(rownames(overlap_matrix))) {
  for (j in seq_along(colnames(overlap_matrix))) {
    overlap_value <- overlap_matrix[i, j]
    if (overlap_value > 0) {
      algcl1 <- rownames(overlap_matrix)[i]
      algcl2 <- colnames(overlap_matrix)[j]
      overlapped_genes <- calculate_overlap_genes(algcl1, algcl2, df1, df2)

      # Update gene frequency table
      gene_frequency_df <- rbind(gene_frequency_df, data.frame(Gene = overlapped_genes, Frequency = 1))
    }
  }
}

# Sum the frequencies for each gene
gene_frequency_df <- aggregate(Frequency ~ Gene, gene_frequency_df, sum)

# Print the table of genes with frequencies
print(gene_frequency_df)


# I need to generate a heatmap now trubetskoy genes in other enrichments,

# I maybe want to re-do with comparing all cluster genes, not just enriched ones

#(1) TrubetskoyFull and Priortised,
# Add in Priortised, just highlight priortised,
#
# higher resolution for enriched clusters, what terms are associated and how many genes, just summary for spectral,
#
# disgard infomap too many lone clusters, sgG1, sgG2. Want a summary sheet for each one what are the enriched genes, what the enriched GO terms. See enriched go terms for spin glass one
# Re-run spin glasses and spectral, just with synGO subset. See if there enough left in syngo for enrichment, is syngo too biased.
#
# flawed assume you only know syngo and re-build network with just syngo stuff, cluster using same algorithms and genetic data, what do you gain / lose for syngo. what is present in both