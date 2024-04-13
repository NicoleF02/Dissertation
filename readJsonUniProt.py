# read in the json from uniprot and get the ids

import json
import pandas as pd

# Read in the json file
with open('idmapping_2024_04_11.json') as f:
    data = json.load(f)
    # get the primaryAccession
    primaryAccession = data['primaryAccession']
    print(primaryAccession)