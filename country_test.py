import geonamescache
from geopy.geocoders import Nominatim
gc = geonamescache.GeonamesCache()
import geocoder
from entrez_esummary_test import combined_dictionaries, user_tool
from nltk.tokenize import word_tokenize
import pandas as pd
from geopy import distance

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

print(combined_dictionaries['Distance to me (km)'])