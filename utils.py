import requests
import pandas as pd
import numpy as np

def get_osm_objects(country='CZ', key='natural', value='peak', extract_tags={'name':'??','ele':''}):
    '''Get pandas dataframe with OSM objects of given type in a country.

    :param country: 'CZ', 'DE', ...
    :param key: e.g. 'natural' or 'amenity'
    :param value: e.g. 'peak' or 'pub'
    :param extract_tags: list of tags that will be extracted to the dataframe
    :return: pandas dataframe


    '''
    overpass_url = "http://overpass-api.de/api/interpreter"
    overpass_query = """
    [out:json];
    area["ISO3166-1"="%s"][admin_level=2];
    (node["%s"="%s"](area);
    );
    out center;
    """ % (country, key, value)

    response = requests.get(overpass_url, params={'data': overpass_query})

    data = response.json()

    df = pd.DataFrame(data['elements'])

    if 'tags' in df:
        for tag, missing in extract_tags.items():
            df[tag] = df['tags'].apply(lambda x: x.get(tag, missing))
        df.drop('tags', inplace=True, axis=1)

    return df