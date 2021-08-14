from flask import Flask, render_template, request

#29189291 33343912 28648412 605691 33822812 32930074 33151926

#local anaesthetics epidural analgesia acetaminophen pain Covid-19 SARS-cov-2

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

def datatable(user_email, user_tool, pmids, keywords):
    from ResearchDashboard import create_icite_dict, create_esummary_file, create_efetch_file, create_fulltextlinks, combine_dataframes, calculate_distance, impact_factor, article_age, first_author_rcr, last_author_rcr, keyword_match, bar_plot, radar_plot
    keyword_ratio = 80
    create_icite_dict(pmids)
    create_esummary_file(pmids, user_email)
    create_efetch_file(pmids, user_email)
    combined_dataframes = combine_dataframes(pmids, user_email)
    create_fulltextlinks(combined_dataframes, user_email)
    calculate_distance(combined_dataframes, user_tool)
    impact_factor(combined_dataframes)
    article_age(combined_dataframes)
    first_author_rcr(combined_dataframes, user_email)
    last_author_rcr(combined_dataframes, user_email)
    keyword_match(combined_dataframes, keyword_ratio, keywords)
    combined_dataframes = combined_dataframes[['Pubmed ID',
                                                   'Title',
                                                   'Journal',
                                                   'Authors',
                                                   'Affiliation',
                                                   'Languages',
                                                   'Publication Type',
                                                   'Year',
                                                   'Article age (years)',
                                                   'Impact Factor',
                                                   'Relative Citation Ratio',
                                                    'First Author RCR',
                                                    'Last Author RCR',
                                                    'Keyword match ratio (%)',
                                                    'Distance to me (km)',
                                                    'MeSH Terms',
                                                    'DOI',
                                                    'Full text link',
                                                    'Abstract']]
    for idx in combined_dataframes.index:
        title = combined_dataframes.loc[idx, 'Title']
        if not title is None:
            combined_dataframes.loc[idx, 'Title'] = title.lstrip("[").rstrip("].")
    for idx in combined_dataframes.index:
        authors = combined_dataframes.loc[idx, 'Authors']
        if not authors is None and len(authors) > 1:
            combined_dataframes.loc[idx, 'Authors'] = ", ".join(authors)
    for idx in combined_dataframes.index:
        affiliation = combined_dataframes.loc[idx, 'Affiliation']
        if isinstance(affiliation, list) and len(affiliation) > 1:
            combined_dataframes.loc[idx, 'Affiliation'] = ", ".join(affiliation)
        else:
            combined_dataframes.loc[idx, 'Affiliation'] = str(affiliation).lstrip("[").rstrip("]")
    for idx in combined_dataframes.index:
        languages = combined_dataframes.loc[idx, 'Languages']
        if not languages is None and len(languages) > 1:
            combined_dataframes.loc[idx, 'Languages'] = ", ".join(languages)
        else:
            combined_dataframes.loc[idx, 'Languages'] = str(languages).lstrip("['").rstrip("']")
    for idx in combined_dataframes.index:
        publication_type = combined_dataframes.loc[idx, 'Publication Type']
        if not publication_type is None and len(publication_type) > 1:
            combined_dataframes.loc[idx, 'Publication Type'] = ", ".join(publication_type)
        else:
            combined_dataframes.loc[idx, 'Publication Type'] = str(publication_type).lstrip("['").rstrip("']")
    for idx in combined_dataframes.index:
        mesh_terms = combined_dataframes.loc[idx, 'MeSH Terms']
        combined_dataframes.loc[idx, 'MeSH Terms'] = str(mesh_terms).lstrip("[").rstrip("]")
    for idx in combined_dataframes.index:
        doi = combined_dataframes.loc[idx, 'DOI']
        if not doi is None:
            combined_dataframes.loc[idx, 'DOI'] = str(doi).lstrip("[").rstrip("]")
    bar_plot(combined_dataframes, pmids)
    radar_plot(combined_dataframes, pmids)
    for x in combined_dataframes:
        data = combined_dataframes.drop(columns=['Abstract','DOI','Full text link']).iloc
    for x in combined_dataframes['Full text link']:
        fulltext_data = combined_dataframes['Full text link'].iloc
    for x in combined_dataframes['DOI']:
        doi_data = combined_dataframes['DOI'].iloc
    for x in combined_dataframes['Abstract']:
        abstract_data = combined_dataframes['Abstract'].iloc
    return render_template('datatable.html', data=data, abstract_data=abstract_data, fulltext_data=fulltext_data, doi_data=doi_data, num_rows=len(pmids))

@app.route('/process_input',methods=['POST'])
def process_input():
    user_email = request.form['user_email']
    user_tool = request.form['user_tool']
    pmids = request.form['pmids'].split(" ")
    keywords = request.form['keywords'].split(" ")
    return datatable(user_email, user_tool, pmids, keywords)

if __name__ == "__main__":
    app.run(debug=True)