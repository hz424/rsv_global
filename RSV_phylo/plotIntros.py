# %%
import pandas as pd
import plotly.express as px

import plotly.graph_objects as go
from plotly.subplots import make_subplots


wdir = '/Users/ahenschel/Dropbox/Flu_Madikay/Scripts/WorldMap'
option = 1

combinations = (['A', 'Introduction'],
                ['A', 'Export'],
                ['B', 'Introduction'],
                ['B', 'Export'])
for strain, direction in combinations:
    # Load the spreadsheet data
    if direction =='Introduction':
        countFile = f'{wdir}/meta_{strain}_Clusters.csv'
    else:
        countFile = f'{wdir}/exports_{strain}.csv'

    data = pd.read_csv(countFile)  # Make sure this file exists and has the right format
    data['origin'] = data['origin'].replace('England', 'United Kingdom')

    if direction != 'Export':
        origs = pd.DataFrame(data.origin.value_counts())
        origs = origs.reset_index()
    else:
        origs = data

    # 2nd option: color code by # introductions from a spec. country
    if option==2:
        origs = pd.DataFrame([(c, len(df), list(df.origin)[0]) for c, df in data.groupby('cluster')])
        origs = origs.dropna()
        origs = pd.DataFrame(origs.value_counts(2))
        origs = origs.reset_index()
        origs.columns = [ 'origin', 'count']
    # Create the choropleth map
    fig = px.choropleth(
        data_frame=origs,
        locations="origin",
        locationmode="country names",
        color="count",
        color_continuous_scale="Viridis",  
        title=f"{direction}s RSV {strain}",
        width=1200,   # width in pixels
        height=800    # height in pixels
    )
    fig.write_image(f'{wdir}/rsv{strain}world_{direction}.svg')


# %%
