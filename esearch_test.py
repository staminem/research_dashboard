import pandas as pd
import requests
from Bio import Entrez
from entrez_esummary_test import combined_dictionaries
Entrez.email = "girsberger.stefan@gmail.com"

def first_author_rcr(combined_dictionaries):
    for idx in combined_dictionaries.index:
        author = combined_dictionaries.loc[idx, 'Authors']
        first_author = author[0]
        handle = Entrez.esearch(db="pubmed", term=first_author + "[author]", retmax=200)
        record = Entrez.read(handle)
        handle.close()
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
        try:
            for x in df_icite_dict['relative_citation_ratio']:
                if x >0:
                    icite_list.append(x)
        except KeyError:
            icite_list.append(0)
        combined_dictionaries.loc[idx, 'First Author RCR'] = round(sum(icite_list), 2)
    return(combined_dictionaries)

first_author_rcr = first_author_rcr(combined_dictionaries)

first_author_rcr.to_excel("test_author_rcr.xlsx")
