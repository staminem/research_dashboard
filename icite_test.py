import pandas as pd
import requests
from esearch_test import record

pmids = record['IdList']
icite_dict = []
for x in pmids:
    response = requests.get(
        "/".join([
            "https://icite.od.nih.gov/api",
            "pubs",
            x,
            ]),
            )
    pub = response.json()
    icite_dict.append(pub)
df_icite_dict = pd.DataFrame.from_dict(icite_dict)
icite_list = []
for x in df_icite_dict['relative_citation_ratio']:
    if x >0:
        icite_list.append(x)
author_rcr = round(sum(icite_list), 2)
print(author_rcr)