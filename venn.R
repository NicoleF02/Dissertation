# library(DBI)
#
# mydb <- dbConnect(RSQLite::SQLite(), "synaptic.proteome_SR_20210408.db")
#
# dbListTables(mydb)
library(data.table)
library(RSQLite)
library(DBI)
library(VennDiagram)
install.packages("ggplot2")
install.packages("dplyr")
library(ggplot2)
library(dplyr)

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

# Only do the coding ones, broad coding and prioritised coding.
# Compare both of them seperately with the 3 current data sets, to get two venn diagrams from them.
# As soon as they are compared, find which pre,post will have more schizophrenia genes for it. Another set with in syn go and not in syn go. Another venn diagram with the database

extractColumnAsList <- function(data_frame, column_header) {
  if (column_header %in% names(data_frame)) {
    column_data <- data_frame[[column_header]]
    column_data <- as.integer(column_data)
    return(column_data)
  } else {
    cat("Column header not found in the table.\n")
    return(NULL)
  }
}

runningTotalMatching <- function(list1, list2){
  total <- 0
  for (value in list1){
    if (any(list2 == value, na.rm = TRUE)){
      total <- total + 1
    }

  }
  return(total)
}

matchingValues <- function(list1, list2){
  matching_list <- vector("numeric")
  for (value in list1){
    if (any(list2 == value, na.rm = TRUE)){
      matching_list <- c(matching_list,value)
    }
  }
  return(matching_list)
}

matchingValuesSynapse <- function(matchedList, synapseDB){
  key <- "HUMAN.ENTREZ.ID"
  matching_rows <- synapseDB[synapseDB[, key] %in% matchedList, ]
  return(matching_rows)
}

vennSynapse <- function(vennSynapse){
  psdpres <- sum(subset(vennSynapse, psd == 1 & pres == 1)$psd)
  pressyn <- sum(subset(vennSynapse, syn == 1 & pres == 1)$syn)
  psdsyn <- sum(subset(vennSynapse, psd == 1 & syn == 1)$psd)
  psdpressyn <- sum(subset(vennSynapse, psd == 1 & pres == 1 & syn == 1)$psd)
  psd <- sum(vennSynapse$psd)
  pres <- sum(vennSynapse$pres)
  syn <- sum(vennSynapse$syn)


  grid.newpage()
  draw.triple.venn(area1 = psd, area2 = pres, area3= syn,
                   n12 = psdpres, n23= pressyn, n13 = psdsyn, n123=psdpressyn,
                   fill = c('yellow','brown','blue'),
                   category = c("PSD","PRES","SYN"))

}

filteredSynGo <- function(totalDB, dbSynList){
  key <- "HumanEntrez"
  matching_rows <- totalDB[totalDB[, key] %in% dbSynList, ]
  return(matching_rows)
}

# Just join gene and disease
# Extract list and filter it for non null, make sure its just for human.
Trubetskoy_2022_broad_coding <- na.omit(extractColumnAsList(trubetskoyFile, "Trubetskoy_2022_broad_coding"))
Trubetskoy_2022_prioritised_coding <- na.omit(extractColumnAsList(trubetskoyFile, "Trubetskoy_2022_prioritised_coding"))
HumanEntrezDB <- na.omit(extractColumnAsList(dbData, "HumanEntrez"))

dbTotal <- length(HumanEntrezDB)

Trubetskoy_2022_broad_coding_total <- length(Trubetskoy_2022_broad_coding)

Trubetskoy_2022_prioritised_coding_total <- length(Trubetskoy_2022_prioritised_coding)

Trubetskoy_broad_db_cross <- matchingValues(Trubetskoy_2022_broad_coding, HumanEntrezDB)
Trubetskoy_broad_db_cross_total <- runningTotalMatching(Trubetskoy_2022_broad_coding, HumanEntrezDB)


Trubetskoy_prioritised_db_cross <- matchingValues(Trubetskoy_2022_prioritised_coding, HumanEntrezDB)
Trubetskoy_prioritised_db_cross_total <- runningTotalMatching(Trubetskoy_prioritised_db_cross, HumanEntrezDB)


# This has the entire database as a flat file, we will use this to check for overlapping genes in general
dbFile <- file.path("Files/Full_DB_Rat_Aut22.txt")
synapseDB <- read.table(file=dbFile, sep="\t", header=TRUE)

synapseList <- na.omit(extractColumnAsList(synapseDB, "HUMAN.ENTREZ.ID"))

synapsePriortisedCompleteOverlap <- runningTotalMatching(Trubetskoy_broad_db_cross, synapseList)


tableSynapsePriortisedOverlap <- matchingValuesSynapse(Trubetskoy_2022_prioritised_coding, synapseDB)
tableSynapseBroadOverlap <- matchingValuesSynapse(Trubetskoy_2022_broad_coding, synapseDB)

#vennSynapse(tableSynapsePriortisedOverlap)

# Want to compare Full_DB_Rat and full database

query <- "SELECT distinct HumanEntrez
FROM FullGeneFullDisease"

dbDataFull <- extractColumnAsList(dbGetQuery(db, query), "HumanEntrez")


# We don't want distinct HumanEntrez as SynGO can be accosciated with more than one
# query <- "SELECT HumanEntrez, SynGO, B.Des
# FROM Gene A
# LEFT JOIN GOGene B
# ON A.
# WHERE HumanEntrez is NOT NULL and SynGO is NOT NULL"
#
# dbDataSynGo <- dbGetQuery(db, query)
# dbDataSynGoList <- extractColumnAsList(dbDataSynGo, "HumanEntrez")
# synMatchBroad <- filteredSynGo(dbDataSynGo, Trubetskoy_2022_broad_coding)

# We want which ones have synGO associations, so yeah, maybe also do most common syngo associations?

#crossSynapseDBRAT <- matchingValuesSynapse(dbDataFull, synapseDB)

query <- "SELECT a.GeneID, b.HumanEntrez, a.GOID
FROM GOGene a
LEFT JOIN Gene b on a.GeneID=b.ID
WHERE GOID is NOT NULL"

dbDataGo <- dbGetQuery(db, query)
dbDataGoList <- extractColumnAsList(dbDataSynGo, "HumanEntrez")
# broadGOCross <- matchingValues(dbDataGoList, Trubetskoy_2022_broad_coding)
# # TODO possible bar on what GO associations
# draw.pairwise.venn(area1 = Trubetskoy_2022_broad_coding_total, area2 = length(dbDataGoList),
#                    cross.area = length(broadGOCross), fill = c('red','blue'),
#                    category = c("Trubetskoy_2022_broad_coding","GO Genes"))


prioristiedGOCross <- filteredSynGo(dbDataGoList, Trubetskoy_2022_prioritised_coding)



data_count <- prioristiedGOCross %>%
  group_by(GOID) %>%
  summarize(count = n())

ggplot(data_count, aes(x = GOID, y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Bar Chart of GOID Occurrences", x = "GOID", y = "Count")





# draw.pairwise.venn(area1 = Trubetskoy_2022_prioritised_coding_total, area2 = length(dbDataGoList),
#                    cross.area = length(prioristiedGOCross), fill = c('red','blue'),
#                    category = c("Trubetskoy_2022_priortised_coding","GO Genes"))



#vennSynapse(tableSynapsePriortisedOverlap)




#
# draw.pairwise.venn(area1 = length(dbDataFull), area2 = length(synapseList),
#                    cross.area = length(crossSynapseDBRAT), fill = c('red','blue'),
#                    category = c("Full_DB_Rat","Synaptic DB"))




# draw.pairwise.venn(area1 = Trubetskoy_2022_prioritised_coding_total, area2 = length(synapseList),
#                    cross.area = synapsePriortisedCompleteOverlap, fill = c('red','blue'),
#                    category = c("Trubetskoy_2022_prioritised_coding_total","Synaptic DB Total"))




# synapseBroadCompleteOverlap <- runningTotalMatching(Trubetskoy_2022_broad_coding, synapseList)
# draw.pairwise.venn(area1 = Trubetskoy_2022_broad_coding_total, area2 = length(synapseList),
#                    cross.area = synapseBroadCompleteOverlap, fill = c('red','blue'),
#                    category = c("Trubetskoy_2022_broad_coding_total","Synaptic DB Total"))

# draw.pairwise.venn(area1 = Trubetskoy_2022_broad_coding_total, area2 = dbTotal,
#                    cross.area = Trubetskoy_broad_db_cross_total, fill = c('red','blue'),
#                    category = c("Trubetskoy_2022_broad_coding_total","Synaptic DB"))


# draw.pairwise.venn(area1 = Trubetskoy_2022_prioritised_coding_total, area2 = dbTotal,
#                    cross.area = Trubetskoy_prioritised_db_cross_total, fill = c('red','blue'),
#                    category = c("Trubetskoy_2022_prioritised_coding_total","Synaptic DB"))

