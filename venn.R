# library(DBI)
#
# mydb <- dbConnect(RSQLite::SQLite(), "synaptic.proteome_SR_20210408.db")
#
# dbListTables(mydb)
library(data.table)
library(RSQLite)
library(DBI)
library(VennDiagram)

ilePathFullDB <- "Files"

dbFile <- file.path("Files/Full_DB_Rat_Aut22.txt")
trubetskoyFilePath <- file.path("Files/Trubetskoy_2022_magma.txt")

synapseDB <- read.table(file=dbFile, sep="\t", header=TRUE)

# Read the data from the text file
data <- readLines(trubetskoyFilePath)

# Initialize a list to store the rows of Table 2
table2 <- list()

# Determine the maximum number of columns
max_cols <- 0

# Iterate through each line in Table 1
for (line in data) {
  # Split the line based on the tab character "\t" to separate the row header
  parts <- strsplit(line, "\t")[[1]]

  # The first part is the row header for Table 2
  row_header <- parts[1]

  # The rest of the parts are space-separated values
  values <- strsplit(parts[-1], " ")[[1]]

  # Update the maximum number of columns if needed
  max_cols <- max(max_cols, length(values))

  # Append the row to Table 2
  table2[[row_header]] <- values
}

# Fill in missing values with NA to create a data frame
trubetskoyFile <- as.data.frame(
  lapply(table2, function(x) {
    if (length(x) < max_cols) {
      c(x, rep(NA, max_cols - length(x)))
    } else {
      x
    }
  })
)

# Rename the columns with the row headers
colnames(trubetskoyFile) <- names(table2)


db <- dbConnect(SQLite(), dbname = "Database/DS_10283_3877/synaptic.proteome_SR_20210408.db.sqlite")

query <- "SELECT distinct HumanEntrez
FROM FullGeneFullDisease
WHERE HDOID = \"DOID:5419\""

dbData <- dbGetQuery(db, query)

# dbData
# synapseDB
# trubetskoyFile

# we trubetskoyFile and dbData venn diagrams.
grid.newpage()

# -1 as it has headers \\double check
dbTotal <- nrow(dbData) - 1

# Only do the coding ones, broad coding and prioritised coding.
# Compare both of them seperately with the 3 current data sets, to get two venn diagrams from them.
# As soon as they are compared, find which pre,post will have more schizophrenia genes for it. Another set with in syn go and not in syn go. Another venn diagram with the database

find_totals <- function(header_name, table_data){
  # Find the index of the specified header in the column names
  header_index <- which(names(table_data) == header_name)
  # Check if the header exists
  if (length(header_index) > 0) {
    # Calculate the total number of non-NA items for the specified header, -1 to remove header in calc
    total_items_for_header <- sum(!is.na(table_data[-1, header_index]))
    return(total_items_for_header)
  } else {
    cat("Header", header_name, "not found in the table.\n")
  }
}

matching_values_total <- function(header_name, table_data, dbData){
  total <- 0
  for (i in 1:nrow(dbData)){
    value <- dbData[i,0]
    if (any(which(names(table_data) == header_name) == value, na.rm = TRUE)){
      total <- total + 1
      cat("matching value is ",value)
    }
  }
  return(total)
}


# Just join gene and disease
# Extract list and filter it for non null, make sure its just for human.
Trubetskoy_2022_broad_coding <- find_totals("Trubetskoy_2022_broad_coding", trubetskoyFile)
Trubetskoy_2022_prioritised_coding <- find_totals("Trubetskoy_2022_prioritised_coding", trubetskoyFile)


Trubetskoy_broad_db_cross = matching_values_total("Trubetskoy_2022_broad_coding", trubetskoyFile, dbData)
print(Trubetskoy_broad_db_cross)

draw.pairwise.venn(area1 = Trubetskoy_2022_broad_coding_total, area2 = dbTotal, cross.area = Trubetskoy_broad_db_cross)