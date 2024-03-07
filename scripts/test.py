import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(
    
    html.Div(className='row', children=[
        
        html.Div(children=[
                html.Label(['Dataset:'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='dropdown_dataset_filter',
                    options=['BG']+['All'],
                    value='All',
                    searchable=False,
                    clearable=False,
                ),
                html.Label('Slider', style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Slider(
                    min=0,
                    max=9,
                    marks={i: 'Label {}'.format(i) if i == 1 else str(i) for i in range(1, 6)},
                    value=5,
                ),
            ], style=dict(width='33.33%')),

        
        html.Div(children=[
                html.Label(['Class filter'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='dropdown_class_filter',
                    options=cldf.class_label.unique().tolist()+['All'],
                    value='All',
                    searchable=True,
                    clearable=False
                ),
                html.Label('Slider', style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Slider(
                    min=0,
                    max=9,
                    marks={i: 'Label {}'.format(i) if i == 1 else str(i) for i in range(1, 6)},
                    value=5,
                ),
            ],style=dict(width='33.33%')),

        html.Div(children=[
                html.Label(['Subclass filter'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='dropdown_subclass_filter',
                    options=cldf.subclass_label.unique().tolist()+['All'],
                    value='All',
                    clearable=False,
                    searchable=True
                ),
                html.Label('Slider', style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Slider(
                    min=0,
                    max=9,
                    marks={i: 'Label {}'.format(i) if i == 1 else str(i) for i in range(1, 6)},
                    value=5,
                ),
            ],style=dict(width='33.33%')),
    ],style=dict(display='flex')),
)

if __name__ == '__main__':
    app.run_server(host='127.0.0.1')