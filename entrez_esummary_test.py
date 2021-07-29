from Bio import Entrez
from Bio import Medline
from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import requests
import numpy as np
import datetime
from impact_factor import ImpactFactor
from geograpy import places
from geopy.geocoders import Nominatim
from geopy.geocoders import Nominatim
from geopy import distance
from nltk.tokenize import word_tokenize
import geonamescache
import pandas as pd
import geocoder


user_email = "girsberger.stefan@gmail.com"
user_tool = "test_for_my_dissertation"
pmids = ["29189291","22809664","33343912","28648412","605691","33822812","32930074","33151926"]

def create_icite_dict(pmids):
    icite_dict = []
    for x in pmids:
        response = requests.get("/".join(["https://icite.od.nih.gov/api",
                                      "pubs",
                                      x]))
        icite_json = response.json()
        icite_dict.append(icite_json)
    df_1 = pd.DataFrame(icite_dict)
    return df_1

def create_esummary_file(pmids):
    Entrez.email = "girsberger.stefan@gmail.com"
    handle = Entrez.esummary(db="pubmed",
                           id= ",".join(pmids),
                           retmode="text")
    records = Entrez.parse(handle)
    df_2 = pd.DataFrame.from_dict(records)
    return df_2

def create_efetch_file(pmids, user_email):
    Entrez.email = user_email
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
                   'FAU':'Full Authors',
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
    df_3 = df_renamed
    return df_3

def combine_dictionaries(pmids, user_email):
    df_1 = create_icite_dict(pmids)
    df_2 = create_esummary_file(pmids)
    df_3 = create_efetch_file(pmids, user_email)
    esummary_file_combined = pd.concat([df_2,df_1,df_3['Affiliation'],df_3['MeSH Terms'],df_3['Abstract'], df_3['Full Authors']], axis=1)
    esummary_file_combined_renamed = esummary_file_combined.rename(columns={'Id' :'Pubmed ID',
                                                                        'Source':'Journal Abbreviated',
                                                                        'Full Authors':'Authors',
                                                                        'LangList':'Languages',
                                                                        'PubTypeList':'Publication Type',
                                                                        'FullJournalName':'Journal',
                                                                        'year':'Year',
                                                                        'relative_citation_ratio':'Relative Citation Ratio',
                                                                        'citation_count':'Citation count',
                                                                        'citations_per_year':'Citations per year',
                                                                        'expected_citations_per_year':'Expected citations per year',
                                                                        'field_citation_rate':'Field citation rate',
                                                                        })
    esummary_file_combined_renamed_dropped = esummary_file_combined_renamed.drop(columns=['Citation count', 'AuthorList', 'Citations per year', 'Expected citations per year', 'Field citation rate', 'Item','PubDate','EPubDate','LastAuthor','Volume','Issue','Pages','NlmUniqueID','ISSN','ESSN','RecordStatus','PubStatus','ArticleIds','History','References','HasAbstract','PmcRefCount','ELocationID','SO','pmid','title','authors','journal','is_research_article','nih_percentile','human','animal','molecular_cellular','apt','is_clinical','provisional','x_coord','y_coord','cited_by_clin','cited_by','references','doi'], axis=1)
    return esummary_file_combined_renamed_dropped

def calculate_distance(combined_dictionaries, user_tool):
    g = geocoder.ip('me')
    gc = geonamescache.GeonamesCache()
    geolocator = Nominatim(user_agent=user_tool)
    for idx in combined_dictionaries.index:
        location_text = combined_dictionaries.loc[idx, 'Affiliation']
        if not isinstance(location_text, float):
            location_tokens = word_tokenize(location_text[-1])
        else:
            location_tokens = ""
        location_tokens_words = [word for word in location_tokens if word.isalpha()]
        filtered_list = []
        for word in location_tokens_words:
            if word not in filtered_list:
                filtered_list.append(word)
        city_list = gc.get_cities()
        city_list_df = pd.DataFrame(city_list)
        city_data = []
        for x in filtered_list:
            city_found = [x for x in city_list_df.iloc[1] if(x in filtered_list)]
            city_data.append(city_found)
        filtered_list_cities = []
        for city in city_data:
            if city not in filtered_list_cities:
                filtered_list_cities.append(city)
        if filtered_list_cities != [] and filtered_list_cities != [[]]:
            for city in filtered_list_cities:
                location = geolocator.geocode(city)
                break
            me = g.latlng
            target_location = (location.latitude, location.longitude)
            distance_km = distance.distance(me, target_location).km
            combined_dictionaries.loc[idx, 'Distance to me (km)'] = round(distance_km, 0)
        else:
            country_list = gc.get_countries()
            country_list_df = pd.DataFrame(country_list)
            country_data = []
            for x in filtered_list:
                country_found = [x for x in country_list_df.iloc[1] if (x in filtered_list)]
                country_data.append(country_found)
            filtered_list_countries = []
            for country in country_data:
                if country not in filtered_list_countries:
                    filtered_list_countries.append(country)
            if filtered_list_countries != []:
                if filtered_list_countries != [[]]:
                    for country in filtered_list_countries:
                        location = geolocator.geocode(country)
                        break
                    me = g.latlng
                    target_location = (location.latitude, location.longitude)
                    distance_km = distance.distance(me, target_location).km
                    combined_dictionaries.loc[idx, 'Distance to me (km)'] = round(distance_km, 0)
            else:
                combined_dictionaries.loc[idx, 'Distance to me (km)'] = "Not found"
    return combined_dictionaries

def impact_factor(combined_dictionaries):
    IF = ImpactFactor()
    for idx in combined_dictionaries.index:
        journal_text = combined_dictionaries.loc[idx, 'Journal']
        impact_factor = IF.search(journal_text)
        if not impact_factor == None:
            impact_factor_df = pd.DataFrame.from_dict(impact_factor)
            combined_dictionaries.loc[idx, 'Impact Factor'] = impact_factor_df['factor'].iloc[-1]
        else:
            combined_dictionaries.loc[idx, 'Impact Factor'] = "None or not found"
    return combined_dictionaries

def article_age(combined_dictionaries):
    now = datetime.datetime.now()
    for idx in combined_dictionaries.index:
        article_year = combined_dictionaries.loc[idx, 'Year']
        article_age = now.year - article_year
        combined_dictionaries.loc[idx, 'Article age'] = article_age
    return combined_dictionaries

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
                if x >0 and not None:
                    icite_list.append(x)
        except KeyError:
            icite_list.append(0)
        combined_dictionaries.loc[idx, 'First Author RCR'] = round(sum(icite_list), 2)
    return(combined_dictionaries)

def last_author_rcr(combined_dictionaries):
    for idx in combined_dictionaries.index:
        author = combined_dictionaries.loc[idx, 'Authors']
        last_author = author[-1]
        handle = Entrez.esearch(db="pubmed", term=last_author + "[author]", retmax=200)
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
                if x >0 and not None:
                    icite_list.append(x)
        except KeyError:
            icite_list.append(0)
        combined_dictionaries.loc[idx, 'Last Author RCR'] = round(sum(icite_list), 2)
    return(combined_dictionaries)


combined_dictionaries = combine_dictionaries(pmids, user_email)
combined_dictionaries = calculate_distance(combined_dictionaries, user_tool)
combined_dictionaries = impact_factor(combined_dictionaries)
combined_dictionaries = article_age(combined_dictionaries)
combined_dictionaries = first_author_rcr(combined_dictionaries)
combined_dictionaries = last_author_rcr(combined_dictionaries)

"""rules_RCR = {
    1: (combined_dataframes['Relative Citation Ratio'] > 0) & (combined_dataframes['Relative Citation Ratio'] < 0.6),
    2: (combined_dataframes['Relative Citation Ratio'] > 0.5) & (combined_dataframes['Relative Citation Ratio'] < 1.0),
    3: (combined_dataframes['Relative Citation Ratio'] > 0.9) & (combined_dataframes['Relative Citation Ratio'] < 1.6),
    4: (combined_dataframes['Relative Citation Ratio'] > 1.5) & (combined_dataframes['Relative Citation Ratio'] < 2.0),
    5: (combined_dataframes['Relative Citation Ratio'] > 1.9),
}
combined_dataframes['RCR Score'] = np.select(rules_RCR.values(), rules_RCR.keys())

rules_recent = {
    1: (now.year - combined_dataframes['Year'] > 20),
    2: (now.year - combined_dataframes['Year'] > 14) & (now.year - combined_dataframes['Year'] < 21),
    3: (now.year - combined_dataframes['Year'] > 9) & (now.year - combined_dataframes['Year'] < 15),
    4: (now.year - combined_dataframes['Year'] > 4) & (now.year - combined_dataframes['Year'] < 10),
    5: (now.year - combined_dataframes['Year'] < 5),
}

combined_dataframes['Recent Score'] = np.select(rules_recent.values(), rules_recent.keys())"""

combined_dictionaries = combined_dictionaries[['Pubmed ID',
                                               'Title',
                                               'Journal Abbreviated',
                                               'Journal',
                                               'Impact Factor',
                                               'Authors',
                                               'First Author RCR',
                                               'Last Author RCR',
                                               'Affiliation',
                                               'Languages',
                                               'Publication Type',
                                               'Year',
                                               'Article age',
                                               'Relative Citation Ratio',
                                               'Distance to me (km)',
                                               'MeSH Terms',
                                               'DOI',
                                               'Abstract']]


"""normalised_for_score = pd.concat([combined_dataframes['Impact Factor'], combined_dataframes['Year'], combined_dataframes['Relative Citation Ratio'], combined_dataframes['Distance to me (km)']])
normalised_for_score_df = MinMaxScaler().fit_transform(np.array(normalised_for_score).reshape(0,5))"""

combined_dictionaries.to_excel('test_esummary.xlsx')