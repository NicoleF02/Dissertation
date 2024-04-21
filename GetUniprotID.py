from Bio import Entrez
from Bio import SeqIO
import pandas as pd

# This code is inspired by the Bioinformatics 1 course, a lot of credit to Ian Simpson
def get_uniprot(entrez_ids):
    Entrez.email = "s1958199@ed.ac.uk"
    entrez_ids_str = ",".join(entrez_ids)
    handle = Entrez.efetch(db="gene", id=entrez_ids_str, retmode="xml")
    records = Entrez.read(handle)
    uniprot_accessions = set()
    for record in records:
        try:
            uniprot_id = record['Entrezgene_comments'][0]['Gene-commentary_products'][0]['Gene-commentary_products'][0]['Gene-commentary_accession']
            uniprot_accessions.add(uniprot_id)
        except KeyError:
            pass

    return list(uniprot_accessions)




# Get entrez_ids from pandas table

filepath = "Ora/Enriched/Trubetskoy2022broadcoding_significant_rows_sorted.csv"

df = pd.read_csv(filepath)

# Get all values from overlapGenes column, convert into list of entrez_ids
entrez_ids = df['overlapGenes'].values.tolist()
entrez_ids = [item for sublist in [x.split(",") for x in entrez_ids] for item in sublist]
entrez_ids = list(set(entrez_ids))
print(entrez_ids)


print(get_uniprot(entrez_ids))