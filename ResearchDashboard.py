from Bio import Entrez
from Bio import Medline
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import datetime
from impact_factor import ImpactFactor
from geopy.geocoders import Nominatim
from geopy import distance
from nltk.tokenize import word_tokenize
import geonamescache
import geocoder
from fuzzywuzzy import fuzz
import plotly.graph_objects as go

def create_icite_dict(pmids):
    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.1)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    icite_dict = []
    for x in pmids:
        response = session.get("/".join(["https://icite.od.nih.gov/api",
                                      "pubs",
                                      x]))
        icite_json = response.json()
        icite_dict.append(icite_json)
    df_1 = pd.DataFrame(icite_dict)
    return df_1

def create_esummary_file(pmids, user_email):
    Entrez.email = user_email
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
    if 'MH' not in df:
        df.insert(23, 'MH', "Not found", allow_duplicates=False)
    if 'AD' not in df:
        df.insert(14, 'AD', "Not found", allow_duplicates=False)
    if 'AB' not in df:
        df.insert(11, 'AB', "Not found", allow_duplicates=False)
    if 'FAU' not in df:
        df.insert(12, 'FAU', "Not found", allow_duplicates=False)
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
    df_3 = df_renamed[['Affiliation',
                       'MeSH Terms',
                       'Abstract',
                       'Full Authors']].fillna(0)
    return df_3

def combine_dataframes(pmids, user_email):
    df_1 = create_icite_dict(pmids)
    df_2 = create_esummary_file(pmids, user_email)
    df_3 = create_efetch_file(pmids, user_email)
    esummary_file_combined = pd.concat([df_2,df_1,df_3], axis=1)
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
    combined_dataframes = esummary_file_combined_renamed.drop(columns=['Citation count', 'AuthorList', 'Citations per year', 'Expected citations per year', 'Field citation rate', 'Item','PubDate','EPubDate','LastAuthor','Volume','Issue','Pages','NlmUniqueID','ISSN','ESSN','RecordStatus','PubStatus','ArticleIds','History','References','HasAbstract','PmcRefCount','ELocationID','SO','pmid','title','authors','journal','is_research_article','nih_percentile','human','animal','molecular_cellular','apt','is_clinical','provisional','x_coord','y_coord','cited_by_clin','cited_by','references','doi'], axis=1)
    return combined_dataframes

def calculate_distance(combined_dataframes, user_tool):
    g = geocoder.ip('me')
    gc = geonamescache.GeonamesCache()
    geolocator = Nominatim(user_agent=user_tool)
    for idx in combined_dataframes.index:
        location_text = combined_dataframes.loc[idx, 'Affiliation']
        if not isinstance(location_text, float) and location_text != "Not found" and location_text != 0:
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
            combined_dataframes.loc[idx, 'Distance to me (km)'] = round(distance_km, 0)
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
                    combined_dataframes.loc[idx, 'Distance to me (km)'] = round(distance_km, 0)
            else:
                combined_dataframes.loc[idx, 'Distance to me (km)'] = None
    return combined_dataframes

def impact_factor(combined_dataframes):
    IF = ImpactFactor()
    for idx in combined_dataframes.index:
        journal_text = combined_dataframes.loc[idx, 'Journal']
        impact_factor = IF.search(journal_text)
        if not impact_factor == None:
            impact_factor_df = pd.DataFrame.from_dict(impact_factor)
            combined_dataframes.loc[idx, 'Impact Factor'] = impact_factor_df['factor'].iloc[-1]
        else:
            combined_dataframes.loc[idx, 'Impact Factor'] = None
    return combined_dataframes

def article_age(combined_dataframes):
    now = datetime.datetime.now()
    for idx in combined_dataframes.index:
        article_year = combined_dataframes.loc[idx, 'Year']
        article_age = now.year - article_year
        combined_dataframes.loc[idx, 'Article age (years)'] = article_age
    return combined_dataframes

def first_author_rcr(combined_dataframes, user_email):
    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.1)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    for idx in combined_dataframes.index:
        author = combined_dataframes.loc[idx, 'Authors']
        first_author = author[0]
        Entrez.email = user_email
        handle = Entrez.esearch(db="pubmed", term=first_author + "[author]", retmax=500)
        record = Entrez.read(handle)
        pmids = record['IdList']
        icite_dict = []
        for x in pmids:
            response = session.get(
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
                if x is not None and x > 0:
                    icite_list.append(x)
        except KeyError:
            icite_list.append(0)
        combined_dataframes.loc[idx, 'First Author RCR'] = round(sum(icite_list), 2)
    return combined_dataframes

def last_author_rcr(combined_dataframes, user_email):
    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.1)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    for idx in combined_dataframes.index:
        author = combined_dataframes.loc[idx, 'Authors']
        last_author = author[-1]
        Entrez.email = user_email
        handle = Entrez.esearch(db="pubmed", term=last_author + "[author]", retmax=500)
        record = Entrez.read(handle)
        pmids = record['IdList']
        icite_dict = []
        for x in pmids:
            response = session.get(
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
                if x is not None and x >0:
                    icite_list.append(x)
        except KeyError:
            icite_list.append(0)
        combined_dataframes.loc[idx, 'Last Author RCR'] = round(sum(icite_list), 2)
    return combined_dataframes

def keyword_match(combined_dataframes, keyword_ratio, keywords):
    for idx in combined_dataframes.index:
        mesh_terms = combined_dataframes.loc[idx, "MeSH Terms"]
        if not isinstance(mesh_terms, float) and not mesh_terms == 0:
            mesh_tokens = word_tokenize(" ".join(mesh_terms))
            mesh_tokens_punct = []
            for mesh_token in mesh_tokens:
                if len(mesh_token) > 1:
                    mesh_word = mesh_token.strip("/")
                    mesh_word = mesh_word.split("/")
                    mesh_word = word_tokenize(",".join(mesh_word))
                    for x in mesh_word:
                        if len(x) > 1:
                            mesh_tokens_punct.append(x)
            mesh_words_filtered = []
            for word in mesh_tokens_punct:
                if word not in mesh_words_filtered:
                    mesh_words_filtered.append(word)
            if not isinstance(keywords, float) and not len(keywords) == 0:
                key_words = word_tokenize(" ".join(keywords))
                if not len(key_words) == 0:
                    keyword_matches = []
                    for mesh_word in mesh_words_filtered:
                        for key_word in key_words:
                            keyword_match = fuzz.ratio(mesh_word.lower(), key_word.lower())
                            if keyword_match >= keyword_ratio:
                                keyword_matches.append(keyword_match)
                    combined_dataframes.loc[idx, "Keyword match ratio (%)"] = round((len(keyword_matches)/len(key_words) * 100),1)
                else:
                    combined_dataframes.loc[idx, "Keyword match ratio (%)"] = 0
            else:
                combined_dataframes.loc[idx, "Keyword match ratio (%)"] = 0
        else:
            combined_dataframes.loc[idx, "Keyword match ratio (%)"] = 0
    return combined_dataframes

def radar_plot(combined_dataframes, pmids):
    normalised_for_score = pd.concat([
    combined_dataframes['Impact Factor'],
    combined_dataframes['Article age (years)'],
    combined_dataframes['Relative Citation Ratio'],
    combined_dataframes['Distance to me (km)'],
    combined_dataframes['First Author RCR'],
    combined_dataframes['Last Author RCR'],
    combined_dataframes['Keyword match ratio (%)']
    ], axis=1)
    normalised_for_score = normalised_for_score.replace(to_replace=["None or not found", "Not found", 0], value=0.0000000001)
    vars_normal = ['Impact Factor','Relative Citation Ratio','First Author RCR','Last Author RCR','Keyword match ratio (%)']
    vars_inv = ['Article age (years)', 'Distance to me (km)']
    normalised_for_score[vars_normal] = (normalised_for_score[vars_normal] - normalised_for_score[vars_normal].min())/(normalised_for_score[vars_normal].max() - normalised_for_score[vars_normal].min())
    normalised_for_score[vars_inv] = (normalised_for_score[vars_inv].max() - normalised_for_score[vars_inv])/(normalised_for_score[vars_inv].max() - normalised_for_score[vars_inv].min())
    categories = ['Impact Factor',
              'Article age (years)',
              'Relative Citation Ratio',
              'Distance to me (km)',
              'First Author RCR',
              'Last Author RCR',
              'Keyword match ratio (%)']
    fig = go.Figure()
    for i in range(len(pmids)):
        fig.add_trace(go.Scatterpolar(
            r=[
                normalised_for_score['Impact Factor'].iloc[i],
                normalised_for_score['Article age (years)'].iloc[i],
                normalised_for_score['Relative Citation Ratio'].iloc[i],
                normalised_for_score['Distance to me (km)'].iloc[i],
                normalised_for_score['First Author RCR'].iloc[i],
                normalised_for_score['Last Author RCR'].iloc[i],
                normalised_for_score['Keyword match ratio (%)'].iloc[i]
                ], theta=categories, fill='toself', name=str(" ".join(word_tokenize(combined_dataframes['Title'].iloc[i])[0:4]) + "...")))
        fig.update_layout(polar=dict(radialaxis=dict(visible=True,range=[0,1])),showlegend=True)
    fig.write_html("/Users/stefangirsberger/PycharmProjects/ResearchDashboard/www/static/dashboard_plot.html", full_html=False, include_plotlyjs='cdn')

def bar_plot(combined_dataframes, pmids):
    normalised_for_score = pd.concat([
    combined_dataframes['Impact Factor'],
    combined_dataframes['Article age (years)'],
    combined_dataframes['Relative Citation Ratio'],
    combined_dataframes['Distance to me (km)'],
    combined_dataframes['First Author RCR'],
    combined_dataframes['Last Author RCR'],
    combined_dataframes['Keyword match ratio (%)']
    ], axis=1)
    normalised_for_score = normalised_for_score.replace(to_replace=["None or not found", "Not found", 0], value=0.0000000001)
    vars_normal = ['Impact Factor','Relative Citation Ratio','First Author RCR','Last Author RCR','Keyword match ratio (%)']
    vars_inv = ['Article age (years)', 'Distance to me (km)']
    normalised_for_score[vars_normal] = (normalised_for_score[vars_normal] - normalised_for_score[vars_normal].min())/(normalised_for_score[vars_normal].max() - normalised_for_score[vars_normal].min())
    normalised_for_score[vars_inv] = (normalised_for_score[vars_inv].max() - normalised_for_score[vars_inv])/(normalised_for_score[vars_inv].max() - normalised_for_score[vars_inv].min())
    categories = ['Impact Factor',
              'Article age (years)',
              'Relative Citation Ratio',
              'Distance to me (km)',
              'First Author RCR',
              'Last Author RCR',
              'Keyword match ratio (%)']
    fig = go.Figure()
    for i in range(len(pmids)):
        fig.add_trace(go.Bar(
        x=categories,
        y=[normalised_for_score['Impact Factor'].iloc[i],
                normalised_for_score['Article age (years)'].iloc[i],
                normalised_for_score['Relative Citation Ratio'].iloc[i],
                normalised_for_score['Distance to me (km)'].iloc[i],
                normalised_for_score['First Author RCR'].iloc[i],
                normalised_for_score['Last Author RCR'].iloc[i],
                normalised_for_score['Keyword match ratio (%)'].iloc[i]
                ],
        name=str(" ".join(word_tokenize(combined_dataframes['Title'].iloc[i])[0:4]) + "...")))
        fig.update_layout(barmode='group', xaxis_tickangle=-45)
    fig.write_html("/Users/stefangirsberger/PycharmProjects/ResearchDashboard/www/static/bar_plot.html", full_html=False, include_plotlyjs='cdn')