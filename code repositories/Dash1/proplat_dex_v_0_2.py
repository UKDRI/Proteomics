# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, dcc, html ,dash_table
import plotly.express as px
import pandas as pd
import math
import numpy as np
from dash.dependencies import Input, Output, State

external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]
app = Dash(__name__,external_stylesheets=external_stylesheets)

colors = {
    'background': '#FFFFFF',
    'text': '#002C6A'
}

# assume you have a "long-form" data frame
# see https://plotly.com/python/px-arguments/ for more options
df = pd.DataFrame({
    "Fruit": ["Apples", "Oranges", "Bananas", "Apples", "Oranges", "Bananas"],
    "Amount": [4, 1, 2, 2, 4, 5],
    "City": ["SF", "SF", "SF", "Montreal", "Montreal", "Montreal"]
})

fig = px.bar(df, x="Fruit", y="Amount", color="City", barmode="group")



fig.update_layout(
    plot_bgcolor=colors['background'],
    paper_bgcolor=colors['background'],
    font_color=colors['text']
)

def generate_table(dataframe, max_rows=10):
    return html.Table([
        html.Thead(
            html.Tr([html.Th(col) for col in dataframe.columns])
        ),
        html.Tbody([
            html.Tr([
                html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
            ]) for i in range(min(len(dataframe), max_rows))
        ])
    ])

def maketable(dataf):
    tab = html.Div([dash_table.DataTable(
                id='adding-rows-table',
                editable=True,
                data=dataf.to_dict('rows'),
                columns=[{'name': i, 'id': i} for i in dataf.columns])])
    return tab


#
# DB load and connection
#


# Example python program to read data from a PostgreSQL table

# and load into a pandas DataFrame

import psycopg2

import pandas as pd

from sqlalchemy import create_engine

 

# Create an engine instance

alchemyEngine   = create_engine('postgresql+psycopg2://postgres:119583rescue@localhost:5432/postgres', pool_recycle=3600);

 

# Connect to PostgreSQL server

dbConnection    = alchemyEngine.connect();

 

dffrsql       = pd.read_sql("select * from public.\"LimmaTable\" WHERE \"ID\" = 'P16858'" , dbConnection);


#
# Volcano plot loading
#
dffrsql.iloc[:,5] = np.negative(np.log10(dffrsql.iloc[:,5]))


figvolc = px.scatter(dffrsql, x = "logFC",y="P.Value")
figvolc.update_layout(xaxis_range=[-4,4])
figvolc.update_layout(yaxis_range=[0,3])
figvolc.update_layout(plot_bgcolor="#FFFFFF")
figvolc.add_hline(y=1.3013, line_width=2, line_dash='dash')
figvolc.add_vline(x=1, line_width=1.5, line_dash='dash')
figvolc.add_vline(x=-1, line_width=1.5, line_dash='dash')
figvolc.add_vline(x=0, line_width=1.5)

# App Layout


# Close the database connection

dbConnection.close();

# Input button search

@app.callback(
    
    #Output("table", component_property="table")
    [Output("dash_neu", "data"),
    Output("volcano", component_property="figure")],
    #],
    [Input("input2", "value")])
    
def update(input2):
    print(input2)
    if input2 is not None:
        dbConnection    = alchemyEngine.connect();
        print("Connected")
        statment= f""" select * from public.\"LimmaTable\" WHERE \"ID\" = '{input2}'"""
        print(statment)
        dffrsql = pd.read_sql(statment, dbConnection);
        dbConnection.close();
        print(dffrsql.head())
        dffrsql.iloc[:,5] = np.negative(np.log10(dffrsql.iloc[:,5]))
        figvolc = px.scatter(dffrsql, x = "logFC",y="P.Value")
        figvolc.update_layout(xaxis_range=[-5,5])
        figvolc.update_layout(yaxis_range=[0,6])
        figvolc.update_layout(plot_bgcolor="#FFFFFF")
        figvolc.add_hline(y=1.3013, line_width=2, line_dash='dash')
        figvolc.add_vline(x=1, line_width=1.5, line_dash='dash')
        figvolc.add_vline(x=-1, line_width=1.5, line_dash='dash')
        figvolc.add_vline(x=0, line_width=1.5)
        data = dffrsql.to_dict('records')
        
        
    return data, figvolc
    
#def update_output(input2):
#    return u'Input 2 {}'.format(input2)

# Read data from PostgreSQL database table and load into a DataFrame instance

pd.set_option('display.expand_frame_repr', False);


app.layout = html.Div(style={'backgroundColor': colors['background'],'textAlign': 'center'}, children=[
    html.Img(src=app.get_asset_url('dataexplorer_title.jpg')),
    
    html.H1(
        children='UKDRI Results Query Tool',
        style={
            'textAlign': 'center',
            'color': colors['text']
        }
    ),
    
    html.Div(
    [
        html.I("Protein search:"),
        html.Br(),
        dcc.Input(id="input2", type="text", placeholder="", debounce=True),
        html.Div(id="output"),
    ]),

    html.Div(style={'textAlign': 'center'}, children =[#maketable(dffrsql),
    #dash_table.DataTable(dffrsql.to_dict('records'), [{"name": i, "id": i} for i in dffrsql.columns]),
    #generate_table(dffrsql),
    
    dash_table.DataTable(id="dash_neu", columns=[{'name': i, 'id': i} for i in dffrsql.columns], data=dffrsql.to_dict('records')),
    html.Br(),
    html.Br(),
    ]),

    html.Div(children='Volcano Plot', style={
        'textAlign': 'center',
        'color': colors['text']
    }),
    
    dcc.Graph(
        id='volcano',
        figure=figvolc
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)