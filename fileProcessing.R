library(data.table)
library(RSQLite)
library(DBI)
library(ggplot2)
library(dplyr)

db <- dbConnect(SQLite(), dbname = "Database/DS_10283_3877/synaptic.proteome_SR_20210408.db.sqlite")


query <- "SELECT
    G.HumanEntrez,
    COUNT(*) AS DuplicateCount
FROM
    DiseaseGene DG
JOIN
    PaperGene PG ON DG.GeneID = PG.GeneID
JOIN
    Gene G ON DG.GeneID = G.ID
WHERE
    DG.HDOID = \"DOID:5419\"
GROUP BY
    G.HumanEntrez
HAVING
    COUNT(*) > 1
    AND COUNT(DISTINCT PG.PaperPMID) > 1"

dbData <- dbGetQuery(db, query)


df <- data.frame(
  "EntrezID_Human" = dbData$HumanEntrez,
  "NoGenes" = dbData$DuplicateCount
)

# Specify the file name
file_name <- "Files/SchizphreniaPriortisedDB.txt"

# Write the data frame to the file
write.table(df, file = file_name, sep = "\t", quote = TRUE, row.names = FALSE)
