library(ggplot2)
library(tidyverse)


data <- read.csv("Ora/Enriched/Trubetskoy2022broadcoding_significant_rows_sorted.csv", stringsAsFactors = FALSE)

data <- data %>%
  separate_rows(overlapGenes, sep = ",\\s*")

data$overlapGenes <- str_trim(data$overlapGenes)

# We want to do an upset graph and venn diagrams

# We want to do an upset plot of all enriched values in the structure and their overlays