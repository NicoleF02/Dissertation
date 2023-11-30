library(data.table)
library(RSQLite)
library(DBI)
library(VennDiagram)
library(ggplot2)
library(dplyr)
library(synaptome.db)


ReadTrubetskoyFile <- function(){
  # Read the data from the text file
  trubetskoyFilePath <- file.path("Files/Trubetskoy_2022_magma.txt")

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

  return(trubetskoyFile)
}


ReadDatabase <- function (query){
  db <- dbConnect(SQLite(), dbname = "Database/DS_10283_3877/synaptic.proteome_SR_20210408.db.sqlite")
  dbData <- dbGetQuery(db, query)

  return(dbData)
}


ExtractColumnAsList <- function(data_frame, column_header) {
  if (column_header %in% names(data_frame)) {
    column_data <- data_frame[[column_header]]
    column_data <- as.integer(column_data)
    return(column_data)
  } else {
    cat("Column header not found in the table.\n")
    return(NULL)
  }
}


RunningTotalMatching <- function(list1, list2){
  total <- 0
  for (value in list1){
    if (any(list2 == value, na.rm = TRUE)){
      total <- total + 1
    }

  }
  return(total)
}


MatchingValues <- function(list1, list2){
  matching_list <- vector("numeric")
  for (value in list1){
    if (any(list2 == value, na.rm = TRUE)){
      matching_list <- c(matching_list,value)
    }
  }
  return(matching_list)
}


MatchingValuesSynapse <- function(matchedList, synapseDB){
  key <- "HUMAN.ENTREZ.ID"
  matching_rows <- synapseDB[synapseDB[, key] %in% matchedList, ]
  return(matching_rows)
}


VennSynapse <- function(vennSynapse){
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


FilteredSynGo <- function(totalDB, dbSynList){
  key <- "EntrezID_Human"
  matching_rows <- totalDB[totalDB[, key] %in% dbSynList, ]
  return(matching_rows)
}


Barchart <- function(dataFrame, titleName, xLabel){
  data_count <- dataFrame %>%
    group_by(SynGO) %>%
    summarize(count = n())

  ggplot(data_count, aes(x = GOID, y = count)) +
    geom_bar(stat = "identity") +
    labs(title = titleName, x = xLabel, y = "Count")
}



dbFile <- file.path("Files/Full_DB_Rat_Aut22.txt")


synapseDB <- read.table(file=dbFile, sep="\t", header=TRUE)

grid.newpage()

query <- "SELECT distinct HumanEntrez
FROM Gene"

dbData <- ReadDatabase(query=query)
HumanEntrezDB <- na.omit(ExtractColumnAsList(dbData, "HumanEntrez"))
dbTotal <- length(HumanEntrezDB)

# Just join gene and disease


trubetskoyFile <- ReadTrubetskoyFile()
# Extract list and filter it for non null, make sure its just for human.
Trubetskoy_2022_broad_coding <- na.omit(ExtractColumnAsList(trubetskoyFile, "Trubetskoy_2022_broad_coding"))
Trubetskoy_2022_prioritised_coding <- na.omit(ExtractColumnAsList(trubetskoyFile, "Trubetskoy_2022_prioritised_coding"))

Trubetskoy_2022_broad_coding_total <- length(Trubetskoy_2022_broad_coding)
Trubetskoy_2022_prioritised_coding_total <- length(Trubetskoy_2022_prioritised_coding)

Trubetskoy_broad_db_cross <- MatchingValues(Trubetskoy_2022_broad_coding, HumanEntrezDB)
Trubetskoy_broad_db_cross_total <- RunningTotalMatching(Trubetskoy_2022_broad_coding, HumanEntrezDB)

Trubetskoy_prioritised_db_cross <- MatchingValues(Trubetskoy_2022_prioritised_coding, HumanEntrezDB)
Trubetskoy_prioritised_db_cross_total <- RunningTotalMatching(Trubetskoy_prioritised_db_cross, HumanEntrezDB)


# This has the entire database as a flat file, we will use this to check for overlapping genes in general
dbFile <- file.path("Files/Full_DB_Rat_Aut22.txt")
synapseDB <- read.table(file=dbFile, sep="\t", header=TRUE)

synapseList <- na.omit(ExtractColumnAsList(synapseDB, "HUMAN.ENTREZ.ID"))

ratDBDiff <- synapseList[!(synapseList %in% HumanEntrezDB)]


missingValues <- synapseDB[synapseDB$HUMAN.ENTREZ.ID %in% ratDBDiff, ]

listSynaptombeDB <- findGenesByEntrez(ratDBDiff)

temp <- listSynaptombeDB$HumanEntrez

syntombeDiff <- ratDBDiff[!(ratDBDiff %in% temp)]

file_name <- "Files/FlatFileExtraGenes.txt"
write.table(missingValues, file = file_name, sep= "\t", quote=FALSE, row.names = FALSE)



synapsePriortisedCompleteOverlap <- RunningTotalMatching(Trubetskoy_broad_db_cross, synapseList)


tableSynapsePriortisedOverlap <- MatchingValuesSynapse(Trubetskoy_2022_prioritised_coding, synapseDB)
tableSynapseBroadOverlap <- MatchingValuesSynapse(Trubetskoy_2022_broad_coding, synapseDB)

SynGOFile <- "Files/SynGO.txt"
SynGO <- read.table(file=SynGOFile, sep="\t", header=TRUE)



broadSynGOCross <- FilteredSynGo(SynGO, Trubetskoy_2022_prioritised_coding)
#unique(data$Product)


draw.pairwise.venn(area1 = length(unique(SynGO$EntrezID_Human)), area2 = length(Trubetskoy_2022_prioritised_coding),
                   cross.area = length(unique(broadSynGOCross$EntrezID_Human)), fill = c('red','blue'),
                   category = c("Synaptic SynGO","Trubetskoy_2022_prioritised_coding"))




file_name <- "Files/Trubetskoy_2022_broad_coding.txt"

# Write the header to the file

data_frame <- data.frame("EntrezID_Human" = Trubetskoy_2022_broad_coding, "Trubetskoy_Broad_coding" = TRUE)

# Append the numbers to the file
write.table(data_frame, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)

file_name <- "Files/Trubetskoy_2022_priortised_coding.txt"

# Write the header to the file
data_frame <- data.frame("EntrezID_Human" = Trubetskoy_2022_prioritised_coding, "Trubetskoy_Priortised_coding" = TRUE)

# Append the numbers to the file
write.table(data_frame, file = file_name, sep= "\t", quote=TRUE, row.names = FALSE)

file_name <- "Files/Trubetskoy_broad_db_cross.txt"
data_frame <- data.frame("EntrezID_Human" = Trubetskoy_broad_db_cross, "Trubetskoy_broad_db_cross" = TRUE)
write.table(data_frame, file = file_name, sep= "\t", quote=TRUE, row.names = FALSE)

#vennSynapse(tableSynapsePriortisedOverlap)
file_name <- "Files/Trubetskoy_prioritised_db_cross.txt"
data_frame <- data.frame("EntrezID_Human" = Trubetskoy_prioritised_db_cross, "Trubetskoy_prioritised_db_cross" = TRUE)
write.table(data_frame, file = file_name, sep= "\t", quote=TRUE, row.names = FALSE)

# Want to compare Full_DB_Rat and full database

# query <- "SELECT HumanEntrez
# FROM Gene
# WHERE SynGO is not NULL"
#
# dbDataFullSynGO <- dbGetQuery(db, query)
#
# # this is a list and not what we want, we want all columns
# # My to-do is all venns having more than one circle
# broadSynGOCross <- filteredSynGo(dbDataFullSynGO, Trubetskoy_2022_broad_coding)

# dataFrame <- broadSynGOCross
#
#
# #dataFrame, titleName, xLabel, groubByValue
# data_count <- dataFrame %>%
#   group_by(SynGO) %>%
#   summarize(count = n())
#
# ggplot(data_count, aes(x = GOID, y = count)) +
#   geom_bar(stat = "identity") +
#   labs(title = titleName, x = xLabel, y = "Count")



# dbDataFullSynGO <- extractColumnAsList(dbGetQuery(db, query), "HumanEntrez")
#
#




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

# query <- "SELECT a.GeneID, b.HumanEntrez, a.GOID
# FROM GOGene a
# LEFT JOIN Gene b on a.GeneID=b.ID
# WHERE GOID is NOT NULL"
#
# dbDataGo <- dbGetQuery(db, query)
# dbDataGoList <- extractColumnAsList(dbDataSynGo, "HumanEntrez")
# # broadGOCross <- matchingValues(dbDataGoList, Trubetskoy_2022_broad_coding)
# # # TODO possible bar on what GO associations
# # draw.pairwise.venn(area1 = Trubetskoy_2022_broad_coding_total, area2 = length(dbDataGoList),
# #                    cross.area = length(broadGOCross), fill = c('red','blue'),
# #                    category = c("Trubetskoy_2022_broad_coding","GO Genes"))
#
#
# prioristiedGOCross <- filteredSynGo(dbDataGo, Trubetskoy_2022_prioritised_coding)











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

