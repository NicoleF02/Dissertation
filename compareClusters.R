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
        overlap_matrix[i, j] <- calculate_overlap(genes_i, genes_j)
      }
    }
  }
  return(overlap_matrix)
}

# Function to calculate gene overlap between two algcl values
calculate_overlap <- function(algcl1, algcl2, df1, df2) {
  genes1 <- unlist(strsplit(df1$gene[df1$algcl == algcl1], " "))
  genes2 <- unlist(strsplit(df2$gene[df2$algcl == algcl2], " "))
  length(intersect(genes1, genes2))
}

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


# df1 <- split_gene_overlap(csv_data1)
# df2 <- split_gene_overlap(csv_data2)
# print("Read CSV")
# overlap_matrix <- compare_diff(df1, df2)

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
  filename = "UpsetPlot/TrubetskoyGOBPIDComparision.png"
)
print("Done")


# I need to generate a heatmap now trubetskoy genes in other enrichments,

# I maybe want to re-do with comparing all cluster genes, not just enriched ones