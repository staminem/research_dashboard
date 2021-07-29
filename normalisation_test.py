from researchdashboard_test import combined_dataframes
import pandas as pd
import plotly.graph_objects as go
from researchdashboard_test import pmids
from nltk import word_tokenize


def radar_plot():
    normalised_for_score = pd.concat([
    combined_dataframes['Impact Factor'],
    combined_dataframes['Article age (years)'],
    combined_dataframes['Relative Citation Ratio'],
    combined_dataframes['Distance to me (km)'],
    combined_dataframes['First Author RCR'],
    combined_dataframes['Last Author RCR'],
    combined_dataframes['Keyword match ratio (%)']
    ], axis=1)
    normalised_for_score = normalised_for_score.replace(to_replace=["None or not found", "Not found"], value=0)
    normalised_for_score = ((normalised_for_score-normalised_for_score.min())/(normalised_for_score.max()-normalised_for_score.min()))*1000
    normalised_for_score['Distance to me (km)'] = 1000 - normalised_for_score['Distance to me (km)']
    normalised_for_score['Article age (years)'] = 1000 - normalised_for_score['Article age (years)']
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
                ], theta=categories, fill='toself', name=str(" ".join(word_tokenize(combined_dataframes['Title'].iloc[i])[0:3]) + "...")))
        fig.update_layout(polar=dict(radialaxis=dict(visible=True,range=[0,1000])),showlegend=True)
    fig.write_html("/Users/stefangirsberger/PycharmProjects/pythonProject/www/static/dashboard_plot.html", full_html=False, include_plotlyjs='cdn')

def bar_plot():
    normalised_for_score = pd.concat([
    combined_dataframes['Impact Factor'],
    combined_dataframes['Article age (years)'],
    combined_dataframes['Relative Citation Ratio'],
    combined_dataframes['Distance to me (km)'],
    combined_dataframes['First Author RCR'],
    combined_dataframes['Last Author RCR'],
    combined_dataframes['Keyword match ratio (%)']
    ], axis=1)
    normalised_for_score = normalised_for_score.replace(to_replace=["None or not found", "Not found"], value=0)
    normalised_for_score = ((normalised_for_score-normalised_for_score.min())/(normalised_for_score.max()-normalised_for_score.min()))*1000
    normalised_for_score['Distance to me (km)'] = 1000 - normalised_for_score['Distance to me (km)']
    normalised_for_score['Article age (years)'] = 1000 - normalised_for_score['Article age (years)']
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
        name=str(" ".join(word_tokenize(combined_dataframes['Title'].iloc[i])[0:3]) + "...")))
        fig.update_layout(barmode='group', xaxis_tickangle=-45)
    fig.write_html("/Users/stefangirsberger/PycharmProjects/pythonProject/www/static/bar_plot.html", full_html=False, include_plotlyjs='cdn')

radar_plot()
bar_plot()