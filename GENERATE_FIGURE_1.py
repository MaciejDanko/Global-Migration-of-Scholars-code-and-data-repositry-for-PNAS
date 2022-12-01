
#%%
##########
# Script to produce figure 1 (world-map with emigration-rates as chloropeth) 
# from the paper "Global Migration of Scholars: Trends and Patterns
# Revealed by Bibliometric Data"
##########

path = "./"
startyear = 2013
endyear = 2017

#%%
##########
# install dependencies:

# %conda install pandas
# %conda install plotly

# to export as pdf we need python-kaleido
# (according to https://plotly.com/python/static-image-export/ )
# %conda install -c conda-forge python-kaleido
#%%
# from distutils.command.config import config
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

#%%

dfmipop = pd.read_csv("AGGREGATED_DATA.CSV")
df_100_countries = pd.read_csv("COUNTRIES_LIST_100.CSV")
list_countrycodes_100 = df_100_countries["ISO 3"].unique()
dfmipop = dfmipop[dfmipop["countrycode"].isin(list_countrycodes_100)]


#%%
print("Check whether the dataset still contains 100 countries: ", dfmipop["countrycode"].nunique(), " unique countries.")


#%%
# aggregate by country over the the years:
dfmigrouped = dfmipop[(dfmipop["year"] >= startyear ) & (dfmipop["year"] <= endyear )].groupby("countrycode").agg(
        {
            "number_of_outmigrations":"mean",
            "padded_population_of_researchers":"mean",
            "year":"mean",
            "countryname":"first"
        })

# calculate rates:
dfmigrouped["outmigrationrate2"] = dfmigrouped["number_of_outmigrations"]/dfmigrouped["padded_population_of_researchers"]
dfmigrouped["Emigration rate"] = dfmigrouped["outmigrationrate2"] * 1000
dfmigrouped["ER"] = dfmigrouped["outmigrationrate2"] * 1000

worldwide_total_researchers_in_period = dfmigrouped["padded_population_of_researchers"].sum()
print(f"{worldwide_total_researchers_in_period=}")

#%%
# create the plot with plotly:
hover_data={
            'countryname':True,
            'outmigrationrate2':True, 
            'padded_population_of_researchers':True,
            'ER':True,
        }

conf = {"scrollZoom": False}
fig = px.choropleth(dfmigrouped.reset_index(),
              locations='countrycode',
              locationmode='ISO-3',
              color='ER',
              range_color=[.0,100],
              color_continuous_scale = [[0.0, 'rgb(170,220,249,0.5)'],
                      [.04, 'rgb(140,180,215,0.5)'],
                      [.12, 'rgb(185,225,115,0.5)'],
                      [.30, 'rgb(51,160,44,0.5)'],
                      [1.00, 'rgb(227,5,5,0.5)']
                     ],
              projection="robinson",
              scope='world',
              hover_data =hover_data,
              #title=f'Emigration rate over 2013 - 2017 <br><sup>100 countries</sup>',#<br><sup>{len(dfmigrouped)=} countries</sup>',
              )

fig.update_traces(
    hovertemplate="<br>".join([
        "%{customdata[0]}",
        "Scholars: %{customdata[2]}",
        "Emigration rate: %{customdata[1]:.3f}",
    ])
)

# choose angle and center so we see everything except antarctica:
fig.update_geos(
    visible=False,showcountries=True,
    center=dict(lon=7, lat=12),
    lataxis_range=[-74,62], lonaxis_range=[-142, 167]
)
# remove the border around the world:
fig.update_traces(marker_line_width=0)

fig.write_html(path +  f"FIGURES/Figure_1_{startyear}-{endyear}.html",config=conf)
# %%
# create pdf:
fig.update_layout(width =900, height=400, font_size=16, 
                  font_family="Sans-Serif",
                  margin_l=0, margin_t=2, margin_b=1, margin_r=0)
fig.write_image(path +  f"FIGURES/Figure_1_{startyear}-{endyear}.pdf")
# %%
