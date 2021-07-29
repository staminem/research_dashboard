from Bio import Medline
from Bio import Entrez
import pandas as pd
import numpy as np

pmids = ["29189291","22809664","33343912","28648412"]

def create_efetch_file():
    Entrez.email = "girsberger.stefan@gmail.com"
    handle = Entrez.efetch(db='pubmed', id=pmids, rettype='medline', retmode='text')
    records = Medline.parse(handle)
    record_list = []
    for record in records:
        record_list.append(record)
    df = pd.DataFrame(record_list)
    df_renamed = df.rename(columns={'PMID':'Pubmed ID',
                   'OWN':'Owner',
                   'STAT':'Status',
                   'DCOM':'Date Completed',
                   'LR':'Date last revised',
                   'IS':'ISSN',
                   'VI':'Volume',
                   'IP':'Issue',
                   'DP':'Date of Publication',
                   'TI':'Title',
                   'PG':'Pagination',
                   'LID':'Location Identifier',
                   'AB':'Abstract',
                   'FAU':'Full Author',
                   'AU':'Author',
                   'AD':'Affiliation',
                   'LA':'Language',
                   'PT':'Publication Type',
                   'PL':'Place of Publication',
                   'TA':'Journal Title Abbreviation',
                   'JT':'Journal Title',
                   'JID':'NLM Unique ID',
                   'RN':'Registry Number/EC Number',
                   'CIN':'Comment in',
                   'MH':'MeSH Terms',
                   'EDAT':'Entrez Date',
                   'MHDA':'MeSH Date',
                   'CRDT':'Create Date',
                   'PHST':'Publication History Status',
                   'AID':'Article Identifier',
                   'PST':'Publication Status',
                   'SO':'Source',
                   'CI':'Copyright Information',
                   'DEP':'Date of Electronic Publication',
                   'SB':'Subset',
                   'AUID':'Author Identifier',
                   'PMC':'PubMed Central Identifier',
                   'OTO':'Other Term Owner',
                   'OT':'Other Term',
                   'COIS':'Conflict of Interest Statement'
                   })
    efetch_file = df_renamed
    return efetch_file

efetch_file = create_efetch_file()

efetch_file.to_excel('test_efetch.xlsx')