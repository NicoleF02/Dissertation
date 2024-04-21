# Synaptic Protein-Protein Interactions for insights into Schizophrenia
This dissertation has generated a great many files. I have created so many data processing scripts so this will explain what does what.

I will also explain the files and the contents. Thanks for reading, and goodluck with dealing with this code, sorry :(.

### Useful files

FinalDissertationNicoleFerguson.pdf: The actual report!! Read this!!

network.R: This is the main file that generates the networks and performs the clustering. It is the most important file in the project.

ORAAnalysis.R: This is the file that generates the ORA analysis and saves the csvs.

ORAProcessing1.py: This is the file that processes the ORA csvs. It drops the columns that have less than 0.05 pval or padj.

ORAProcessing2.py: This is the file that refines the previously generated csvs from ORAProcessing1.py. It combines all 
enriched data into one mega csv, only for clusters with enriched trubetskoy genes.

wordCloud.py: This is the file that generates all of the word clouds for the BioNAR pipeline.

venn.R: This is the file that generates the Venn Diagrams for the BioNAR pipeline. Lots of them are commented out but just uncomment.

spreadsheetSplitter.py: This is the file that splits the massive ORAProcess2 csv rows into columns. It was used for a dissertation meeting.



### Folders

Bridgeness: This directory contains all of the bridgeness data generated from BioNAR. The .png are the correct ones, the 
pdfs were generated with faulty data.

Consensus, SynGO, Ora: These directories contain the data generated from the BioNAR pipeline. The data is in the form of 
.csv files. It includes the algorithms, their enriched clusters and their enriched terms. The Ora folder is for the full
network.

Database: This contains the database from the universities file share.

Files: Miscelinous data files that include things like the flat file of the database, trubetskoy broad and prioritised 
coding and an list of other schizophrenic associations from the database.

FullDBClustered: Some older ORA ran stuff on the full database.

HeatmapEnrichment: This contains the heatmap figures for all networks and for broad and prioritised, full and reduced.

OldCSVGenerators: This contains the old scripts that generated the csv files for the BioNAR pipeline.

PostsynapticNetwork: the gml files used that represent the PSD networks.

tempCSV: This contains the csv files that were generated from the BioNAR pipeline, older and useless placeholders.

UpsetPlot: An abandoned project to create upsets plots.

VennDiagrams: VennDiagrams generated that were used to see where the project would go. 

venv: virtual environment for my dataspell IDE.

WordClouds: Wordclouds generated from the BioNAR pipeline for all networks and for disease, biological functions, 
molecular functions, and SynGO annotations.


### CSV Generators

csvProcessing.py: Read in the robustness csv generated from bridgeness and split it into different algorithm csvs

fileProcessing.py: Older checker used to try and identify what was wrong with my network early in project. 

fileProcessing.R: As above. Made use of R with SQLite to figure out a bug early on. Functionally useless now.

fileProcessing2.py: Used to get a quick summary of algorithm, no. trubetskoy, pval and padj.

FileProcessing3.py: Used to find maximum cluster numbers quickly.

quickGraph.py: used to generate CommuntiespValue.png, functionally useless.


clusterCompareModularity.py: Calculates modularity for a network. Useless as BioNAR has built it, wanted to double check.

clusterComparison.py: Printed the genes that overlapped in a frequency table.

ClusterGeneVisualise.py: Generates the Heatmaps. Actually Useful.

compareClusters.R: A different heatmap generator. Actually generated data for overlapped genes.

compareClustersRandIndex.py: Useful. This does the ARI and jaccard combination comparison!

compareEnrichedClusters.py: Useful. Does jaccard similarity for enriched clusters.

csvRowDelimmer.py: Older decrepted script that was used to delimit the csv files.

dataProcessing.R: Old script that used to double check something and bug find.

dissertation.R: Old script that was used to generate a simple venn diagram. useless.

enrichment.R: Intial ORA clustering script.

enrichmentReduced.R: This was used to bug hunt for fitsigmoid.

GetUniprotID.py: Unused, and didn't work, just for BLAST things.

network.R: Actually generates the networks and performs the clustering.

OraAnalysis.R: Used to generate ORA analysis and save the csvs.

ORAProcessing.py: Used to process ORA csv, drops the columns have less than 0.05 pval or padj.

ORAProcessing2.py: Used to refine the previously generated csvs from ORAProcessing.py. It combines all enriched data into
one mega csv, only for clusters with enriched trubetskoy genes.

practiceUpset.R: unused, useless.

readJsonUniProt.py: unused, useless.

spreadsheetSplitter.py: Used to split the massive ORAProcess2 csv rows into columns, used for a dissertation meeting.

venn.R: Generates the Venn Diagrams for the BioNAR pipeline. Lots of them are commented out but just uncomment.

wordCloud.py: generates all of the word clouds for the BioNAR pipeline.

### Random Files

BLASTUniProt.list: A list of the proteins that were used in the BLAST search.

CommuntiespValue.png: Unused plot, just showed as no. communties increased avg. pval decreased.

diss.ipynb: The jupyter notebook that was orginally going to be used instead of so many python files. Jupyter did not play
nice with my pc or my laptop and it required logging in every time or manually setting up the kernel. I gave up on it.

go-basic.obo: The Gene Ontology file that was used in the BioNAR pipeline.

README.md: Hey thats here! :)
