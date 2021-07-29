from ResearchDashboard import create_icite_dict, create_esummary_file, create_efetch_file, combine_dataframes, calculate_distance, impact_factor, article_age, first_author_rcr, last_author_rcr, keyword_match, bar_plot, radar_plot



user_email = "girsberger.stefan@gmail.com"
user_tool = "test_for_my_dissertation"
pmids = ["33343912","32930074"]
#["29189291","33343912","28648412","605691","33822812","32930074","33151926"]
keywords = ['local anaesthetics', 'epidural', 'analgesia', 'acetaminophen', 'pain', 'Covid-19', 'SARS-cov-2']
keyword_ratio = 80

create_icite_dict(pmids)
create_esummary_file(pmids)
create_efetch_file(pmids, user_email)
combined_dataframes = combine_dataframes(pmids, user_email)
calculate_distance(combined_dataframes, user_tool)
impact_factor(combined_dataframes)
article_age(combined_dataframes)
first_author_rcr(combined_dataframes)
last_author_rcr(combined_dataframes)
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
                                           'Abstract']].fillna(0)
bar_plot(combined_dataframes, pmids)
radar_plot(combined_dataframes, pmids)
"""combined_dataframes.to_excel("test_researchdashboard.xlsx")"""