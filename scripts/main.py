from dash import Dash, html, dcc, callback, Output, Input, dash_table
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
#from dash import dash_table as dt
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
import pandas as pd
from dash import dash_table as dt
from PIL import Image
import PIL
import requests
from io import BytesIO
import matplotlib
import urllib
import io
import urllib.request
import re
import os
import dash
import pandas as pd
import seaborn as sns
import networkx
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from itertools import compress
from networkx.drawing.nx_agraph import graphviz_layout
import seaborn as sns
import json
import s3fs
#initial setup ---------------------------------------------------




test = False
if test == True:
    exec(open("./scripts/set_up_graph.py").read())#
    data_file = 's3://mh-allen-gt-wb/app/Allen_GT_WB/data/SmartSeq_all_annotation_hmapped.csv'
    cldf = pd.read_csv('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/AIT21_updated_cldf_for_BG_with_parent.csv')
    clus = pd.read_table('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/WB_colorpal - clusters 230815.tsv')
    sub = pd.read_table('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/WB_colorpal - subclasses 230815.tsv')
    clas = pd.read_table('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/WB_colorpal - classes 230815.tsv')
    sup = pd.read_table('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/WB_colorpal - supertypes 230815.tsv')
    MER = pd.read_feather('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/app_MERFISH_data.feather')
elif test == False:
    #os.chdir('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/WB_hierarchy/')
    exec(open("./scripts/set_up_graph.py").read())#
    data_file = 's3://mh-allen-gt-wb/app/Allen_GT_WB/data/SmartSeq_all_annotation_hmapped.csv'
    cldf = pd.read_csv('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/AIT21_updated_cldf_for_BG_with_parent.csv')
    clus = pd.read_table('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/WB_colorpal - clusters 230815.tsv')
    sub = pd.read_table('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/WB_colorpal - subclasses 230815.tsv')
    clas = pd.read_table('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/WB_colorpal - classes 230815.tsv')
    sup = pd.read_table('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/WB_colorpal - supertypes 230815.tsv')
    MER = pd.read_feather('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/app_MERFISH_data.feather')


#### read data and set up colors
clus.columns = clus.columns.str.replace('cluster_label','cluster_id_label')
clus = clus[['cluster_id_label','cluster_color']]
sub = sub[['subclass_id_label','subclass_color']]
clas = clas[['class_id_label','class_color']]
sup = sup[['supertype_id_label','supertype_color']]

#add colors
cldf = cldf.merge(right = clus[['cluster_id_label','cluster_color']], on = 'cluster_id_label', how = 'left')
cldf = cldf.merge(right = sub[['subclass_id_label','subclass_color']], on = 'subclass_id_label', how = 'left')
cldf = cldf.merge(right = sup[['supertype_id_label','supertype_color']], on = 'supertype_id_label', how = 'left')
cldf = cldf.merge(right = clas[['class_id_label','class_color']], on = 'class_id_label', how = 'left')





#### read data and set up colors
root = pd.DataFrame({'id' : ['0_WB'],'node_var' : [root_value]})



#read data


counts_data = pd.read_csv(data_file)

#table for filtering
df = cldf[['cl','cluster_id_label', 'subclass_id_label', 'class_id_label','anatomical_annotation']] 
df = df[df['subclass_id_label'] != "0_WB"]
df['id'] = df['cl']


#cleanup - if id_label columns don't exist, use regular label columns
#counts_data.columns = counts_data.columns.str.replace('_WB_map','')
if not any(counts_data.columns.str.contains('cluster_id_label')):
    counts_data['cluster_id_label'] = counts_data['cluster_label']
if not any(counts_data.columns.str.contains('class_id_label')):
    counts_data['class_id_label'] = counts_data['class_label']
if not any(counts_data.columns.str.contains('subclass_id_label')):
    counts_data['subclass_id_label'] = counts_data['subclass_label']





cldf.index = cldf['cl']


counts_data['LabTracksID'] = counts_data['external_donor_name'].astype('str')
counts_data['LabTracksID'] = counts_data['LabTracksID'].astype('str')
counts_data = counts_data.fillna('')
counts_data["unique_condition"] = counts_data[["LabTracksID",'full_genotype','injection_materials',"facs_population_plan",'Region','roi']].agg(" ".join, axis=1)
counts_data = counts_data.merge(right = cldf[['cluster_id_label','cluster_color','subclass_color','supertype_color','class_color']], on = 'cluster_id_label')
#temporary
#counts_data['unique_condition'] = counts_data['library_prep']


#may not need this later- for click data display
styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}



# app-------------------------------------------------------------
app = Dash(__name__)

app.layout = html.Div([
    html.H1(children='WB data viewer', style={'textAlign':'center'}),
    
# input-----------------------------------------------------------
    html.Div(className='row', children=[
        
        html.Div([ 
            html.Label(['Filter'], style={'font-weight': 'bold', "text-align": "center"}),
            html.Button("clear selection", id="clear"),
            dbc.Container([
                dash_table.DataTable(
                	page_action='none',
                	style_table={'height': '300px', 'overflowY': 'auto'},

                    columns=[
                    {'name': 'Anatomy', 'id': 'anatomical_annotation', 'type': 'text'},
                    {'name': 'Class', 'id': 'class_id_label', 'type': 'text'},
                    {'name': 'Subclass', 'id': 'subclass_id_label', 'type': 'text'},
                    {'name': 'Cluster', 'id': 'cluster_id_label', 'type': 'text'},
                    {'name': 'cl', 'id': 'cl', 'type': 'text'}],

                    data=df.to_dict('records'),
                    filter_action='native',
                    style_data={
                    'width': '150px', 'minWidth': '150px', 'maxWidth': '150px',
                    'overflow': 'hidden',
                    'textOverflow': 'ellipsis'},
                    id='tbl')	,
                dbc.Alert(id='tbl_out')
                ])

            ]),
        
        html.Div(children=[
                html.Label(['Taxonomy subset:'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='dropdown_taxonomy_filter',
                    options=['TH','BG']+['All'],value='All',
                    searchable=False,clearable=False,
                ),
                html.Label(['EC ID'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='dropdown_dataset_filter',
                    options=counts_data.unique_condition.unique().tolist()+['All'],
                    value='All',
                    searchable=True,
                    clearable=False
                ),

                
            ], style=dict(width='33.33%')),


        html.Div(children=[
                html.Label(['Resolution'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='resolution',
                    options=['class_id_label','subclass_id_label','cluster_id_label']+['All'],
                    value='subclass_id_label',
                    clearable=False,
                    searchable=True
                ),

            ],style=dict(width='33.33%')),
    ],style=dict(display='flex')),

    # Flex container - taxonomy graphs
    html.Div([
#
        # Graph container
        html.Div([

            dcc.Graph(id = 'graph1'),

        ], style={'width': '33%', 'display': 'inline-block'}),

        html.Div([

            dcc.Graph(id = 'graph2'),

        ], style={'width': '33%', 'display': 'inline-block'}),

        html.Div([

            dcc.Graph(id = 'graph3'),

        ], style={'width': '33%', 'display': 'inline-block'}),

        #Table container
        html.Div([
            dbc.Container([
                dash_table.DataTable(id = 'tbl_sub')
                ])


        ], style={'width': '49%', 'display': 'inline-block'}),

    ], style={'display': 'flex'}),
    html.Div([
        html.Label('MERFISH Section', style={'font-weight': 'bold', "text-align": "center"}),
        dcc.Slider(
            id = 'mer_slider',
            min=0,
            max=len(MER.section.unique()),
            step = 1,
            value=5,
            ),]),
#
        # MERFISH container
        html.Div([dcc.Graph(id = 'graph_mer'),], style={'width': '40%','height': '40%', 'display': 'inline-block'}),
        html.Div([dcc.Graph(id = 'graph_mer2'),], style={'width': '40%','height': '40%', 'display': 'inline-block'}),

        #display click data

        html.Div([
            dcc.Markdown("""
                **Click Data**

                Click on points in the graph.
            """),
            html.Pre(id='click-data', style=styles['pre']),
        ], className='three columns'),




    ],)

#display click data
@callback(
    Output('click-data', 'children'),
    Input('graph3', 'clickData'))
def display_click_data(clickData):
    return json.dumps(clickData, indent=2)




#clear table
@app.callback(
    Output('dropdown_taxonomy_filter', 'value'),
    Output("tbl", "derived_virtual_row_ids"),
    Output("tbl", "selected_row_ids"),
    Output("tbl", "active_cell"),
    Input("clear", "n_clicks"))
def clear(n_clicks):

    return "All",None, None,"All"

@callback(
    Output('graph1', 'figure'),
    Output('graph2','figure'),
    Input('dropdown_taxonomy_filter', 'value'),
    Input('dropdown_dataset_filter', 'value'),
    Input('tbl', 'derived_virtual_row_ids'),
    Input('tbl', 'selected_row_ids'),
    Input('tbl', 'active_cell'),
    Input('resolution', 'value'))
#@callback(Output('tbl_sub', 'children'), Input('tbl', 'active_cell'))
def update_graphs(taxonomy_filter,dataset_filter, row_ids, selected_row_ids,active_cell,resolution, cldf = cldf, df = df, counts_data = counts_data ):
    #subset to dataset of interest
    #get counts and percents by unique condition, currently library_prep
    
    # if dataset_filter on, get counts by unique condition, else, get counts for every expt
    #if dataset_filter == "All":
    #	counts_data['unique_condition'] = "All"
    counts_data_cluster = counts_data
    counts_data_subclass = counts_data
    cldf_cluster = compute_counts_return_cldf(groupBy = "cluster_id_label", df = counts_data, cldf = cldf, dataset_filter = dataset_filter)
    cldf_select = compute_counts_return_cldf(groupBy = resolution, df = counts_data, cldf = cldf, dataset_filter = dataset_filter)
    #subset graph
    selected = active_cell if active_cell else "All"
    if selected == "All":
        sub_df_graph = cldf_cluster
    else:
        groupBy = selected['column_id']
        active_row_id = active_cell['row_id'] if active_cell else None
        selected_value = pd.DataFrame(df).loc[selected['row_id'], groupBy]
        sub_df_graph = cldf_cluster[(cldf_cluster[groupBy] == selected_value) | (cldf_cluster['source'] == "WB")]

    if ((taxonomy_filter != "All") ):
    	classes = sub_df_graph[sub_df_graph[taxonomy_filter]==True].class_id_label.unique().tolist()+["0_WB"]
    	sub_df_graph = sub_df_graph[(sub_df_graph[taxonomy_filter]==True) | ( (sub_df_graph['class_id_label'].isin(classes)) & (sub_df_graph['source'] == "WB")  ) ]
    if ( (selected == "All") & (taxonomy_filter == "All") & (dataset_filter == "All")):
    	graph1 = orig_full_graph
    	graph2 = orig_sub_graph
    else:
        clusters =sub_df_graph.cluster_id_label.unique().tolist()
        graph1 = build_plotly_taxonomy_full_graph(graph=G,data = cldf_cluster,clusters = clusters)
        graph2 = build_plotly_taxonomy_sub_graph(graph=G,data = cldf_cluster,clusters = clusters)
    #graph3 = build_plotly_bar(data = cldf_select,groupBy = resolution)
    return graph1, graph2


@callback(
    Output('graph3','figure'),
    Input('dropdown_taxonomy_filter', 'value'),
    Input('dropdown_dataset_filter', 'value'),
    Input('tbl', 'derived_virtual_row_ids'),
    Input('tbl', 'selected_row_ids'),
    Input('graph3', 'clickData'),
    Input('resolution', 'value'))
#@callback(Output('tbl_sub', 'children'), Input('tbl', 'active_cell'))
def make_bar_chart(taxonomy_filter,dataset_filter, row_ids, selected_row_ids,clickData,resolution, cldf = cldf, df = df, counts_data = counts_data ):
    #subset to dataset of interest
    cldf_select = compute_counts_return_cldf(groupBy = resolution, df = counts_data, cldf = cldf, dataset_filter = dataset_filter)
    try:
        clickData
        print(clickData)
        graph3 = build_plotly_bar(data = cldf_select,groupBy = resolution, clickData = clickData)
    except NameError:
        graph3 = build_plotly_bar_init(data = cldf_select,groupBy = resolution)
    else:
        graph3 = build_plotly_bar_init(data = cldf_select,groupBy = resolution)
    return graph3






@callback(
    Output('graph_mer', 'figure'),
    Input('mer_slider', 'value'),
    Input('dropdown_taxonomy_filter', 'value'))




def update_merfish_1(mer_slider, taxonomy_filter,
    cldf = cldf, MER = MER):
    if ((taxonomy_filter != "All") ):
        cldf = cldf[(cldf[taxonomy_filter]==True)    ]
        clusters =list(set(cldf.cluster_id_label.unique().tolist()) - set(["0_WB"]))
    else:
        clusters = MER.cluster_id_label.unique().tolist()
    #subset for clusters and section
    MER_sub = MER[(MER['section'] == mer_slider) & MER['cluster_id_label'].isin(clusters)]
    layout = go.Layout(
    title = 'MERFISH',
    xaxis = go.XAxis(
        title = '',
        visible=False,
        showticklabels=False),
    yaxis = go.YAxis(
        title = '',
        visible=False,
        showticklabels=False),
    )
    fig = go.Figure(data=go.Scatter(x=MER_sub['x'], y=MER_sub['y'], mode='markers', marker_color = MER_sub['class_color'],
        text=MER_sub['class_id_label']), layout = layout)

    fig.update_layout(
    autosize=False,
    width=1200,
    height=800,
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)'
    )
    #fig = sns.scatterplot(data = MER_sub, x= 'x', y = 'y', hue = "class_id_label")
    fig['layout']['yaxis']['autorange'] = "reversed"
    return fig

@callback(
    Output('graph_mer2', 'figure'),
    Input('graph3', 'clickData'),
    Input('resolution', 'value'),
    Input('mer_slider', 'value'),
    Input('dropdown_taxonomy_filter', 'value'),
    Input('tbl', 'active_cell'))


def update_merfish_2(clickData,resolution,mer_slider, taxonomy_filter,active_cell,
    cldf = cldf, df = df, MER = MER):
    if clickData:
        print (clickData["points"][0]["label"])
        selected_group = clickData["points"][0]["label"]
        MER['point_color'] = 'lightgray'
        indexer = MER[MER[resolution] == selected_group].index
        MER.loc[indexer, 'point_color'] = 'red'
        MER  = MER.sort_values('point_color', ascending = True)
    else:

        selected = active_cell if active_cell else "All"

        #adjust colors for selected group
        if(selected != "All"):
            groupBy = selected['column_id']
            active_row_id = active_cell['row_id'] if active_cell else None
            selected_value = pd.DataFrame(df).loc[selected['row_id'], groupBy]
            MER['point_color'] = 'lightgray'
            indexer = MER[MER[groupBy] == selected_value].index
            MER.loc[indexer, 'point_color'] = 'red'
            MER  = MER.sort_values('point_color', ascending = True)
    if ((taxonomy_filter != "All") ):
        cldf = cldf[(cldf[taxonomy_filter]==True)    ]
        clusters =list(set(cldf.cluster_id_label.unique().tolist()) - set(["0_WB"]))
    else:
        clusters = MER.cluster_id_label.unique().tolist()
    #subset for clusters and section
    MER_sub = MER[(MER['section'] == mer_slider) & MER['cluster_id_label'].isin(clusters)]
    layout = go.Layout(
    title = 'MERFISH - selected group',
    xaxis = go.XAxis(
        title = '',
        visible=False,
        showticklabels=False),
    yaxis = go.YAxis(
        title = '',
        visible=False,
        showticklabels=False),
    )

    fig = go.Figure(data=go.Scatter(x=MER_sub['x'], y=MER_sub['y'], mode='markers', marker_color = MER_sub['point_color'],
        text=MER_sub['class_id_label']), layout = layout)
    fig.update_layout(
    autosize=False,
    width=1200,
    height=800,
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)'
    )

    #fig = sns.scatterplot(data = MER_sub, x= 'x', y = 'y', hue = "class_id_label")
    fig['layout']['yaxis']['autorange'] = "reversed"
    return fig


@callback(
    Output('img', 'figure'),
    Input('img_slider', 'value')
    )


def update_merfish_2(img_slider,active_cell):
    img_sub = img_data.iloc[img_slider]
    img = get_lims_image(storage_directory = img_sub.storage_directory, barcode = img_sub.barcode)
    #response = requests.get(url)
    #img = Image.open(BytesIO(response.content))
    fig = px.imshow(np.array(img))
    # Hide the axes and the tooltips
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(t=2, b=2, l=2, r=2),
        xaxis=dict(
            showgrid=False,
            showticklabels=False,
            linewidth=0
        ),
        yaxis=dict(
            showgrid=False,
            showticklabels=False,
            linewidth=0
        ),
        hovermode=False
    )
    return fig



if __name__ == '__main__':
    app.run_server(host='127.0.0.1')

