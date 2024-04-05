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

#initial setup ---------------------------------------------------
<<<<<<< HEAD
#os.chdir('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/WB_hierarchy/')
#os.chdir('/home/mh/app/for_GRC/WB_hierarchy/') #system76
os.chdir('/home/mh/app/WB_hierarchy/')#local
=======
os.chdir('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/WB_hierarchy/')
#os.chdir('/home/mh/app/for_GRC/WB_hierarchy/') #system76
#os.chdir('/home/mh/app/WB_hierarchy/')#local
>>>>>>> update python version

exec(open("./scripts/set_up_graph.py").read())#local



#read data
<<<<<<< HEAD
cldf = pd.read_csv('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/Allen_GT_WB_data/AIT21_updated_cldf_for_BG_with_parent.csv')
cldf.index = cldf['cl']
data_file = '/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/Allen_GT_WB_data/test_dataset.feather'
data_file = '/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/Allen_GT_WB_data/all_mapped_metadata_combined_1M_cell_subsample.feather'
data_file = '/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/Allen_GT_WB_data/SmartSeq_all_annotation_hmapped.csv'
=======
data_file = './data/SmartSeq_all_annotation_hmapped.csv'
>>>>>>> update python version


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

<<<<<<< HEAD
=======

cldf = pd.read_csv('./data/AIT21_updated_cldf_for_BG_with_parent.csv')
#### read data and set up colors
clus = pd.read_table('./data/WB_colorpal - clusters 230815.tsv')
clus.columns = clus.columns.str.replace('cluster_label','cluster_id_label')
clus = clus[['cluster_id_label','cluster_color']]
sub = pd.read_table('./data/WB_colorpal - subclasses 230815.tsv')
sub = sub[['subclass_id_label','subclass_color']]
clas = pd.read_table('./data/WB_colorpal - classes 230815.tsv')
clas = clas[['class_id_label','class_color']]
sup = pd.read_table('./data/WB_colorpal - supertypes 230815.tsv')
sup = sup[['supertype_id_label','supertype_color']]
root = pd.DataFrame({'id' : ['0_WB'],'node_var' : [root_value]})

#add colors
cldf = cldf.merge(right = clus[['cluster_id_label','cluster_color']], on = 'cluster_id_label', how = 'left')
cldf = cldf.merge(right = sub[['subclass_id_label','subclass_color']], on = 'subclass_id_label', how = 'left')
cldf = cldf.merge(right = sup[['supertype_id_label','supertype_color']], on = 'supertype_id_label', how = 'left')
cldf = cldf.merge(right = clas[['class_id_label','class_color']], on = 'class_id_label', how = 'left')

cldf.index = cldf['cl']


>>>>>>> update python version
counts_data['LabTracksID'] = counts_data['external_donor_name'].astype('str')
counts_data['LabTracksID'] = counts_data['LabTracksID'].astype('str')
counts_data = counts_data.fillna('')
counts_data["unique_condition"] = counts_data[["LabTracksID",'full_genotype','injection_materials',"facs_population_plan",'Region','roi']].agg(" ".join, axis=1)
<<<<<<< HEAD
=======
counts_data = counts_data.merge(right = cldf[['cluster_id_label','cluster_color','subclass_color','supertype_color','class_color']], on = 'cluster_id_label')
>>>>>>> update python version
#temporary
#counts_data['unique_condition'] = counts_data['library_prep']


<<<<<<< HEAD
=======
#may not need this later- for click data display
styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

>>>>>>> update python version


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
<<<<<<< HEAD
                html.Label(['data filter'], style={'font-weight': 'bold', "text-align": "center"}),
=======
                html.Label(['EC ID'], style={'font-weight': 'bold', "text-align": "center"}),
>>>>>>> update python version
                dcc.Dropdown(
                    id='dropdown_dataset_filter',
                    options=counts_data.unique_condition.unique().tolist()+['All'],
                    value='All',
                    searchable=True,
                    clearable=False
                ),

                
            ], style=dict(width='33.33%')),

<<<<<<< HEAD
        
        # html.Div(children=[
        #         html.Label(['data filter'], style={'font-weight': 'bold', "text-align": "center"}),
        #         dcc.Dropdown(
        #             id='dropdown_dataset_filter',
        #             options=counts_data.unique_condition.unique().tolist()+['All'],
        #             value='All',
        #             searchable=True,
        #             clearable=False
        #         ),
                
        #     ],style=dict(width='33.33%')),

        html.Div(children=[
                html.Label(['Placeholder'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='dropdown_subclass_filter',
                    options=cldf.subclass_id_label.unique().tolist()+['All'],
                    value='All',
=======

        html.Div(children=[
                html.Label(['Resolution'], style={'font-weight': 'bold', "text-align": "center"}),
                dcc.Dropdown(
                    id='resolution',
                    options=['class_id_label','subclass_id_label','cluster_id_label']+['All'],
                    value='subclass_id_label',
>>>>>>> update python version
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

<<<<<<< HEAD
        ], style={'width': '49%', 'display': 'inline-block'}),
=======
        ], style={'width': '33%', 'display': 'inline-block'}),
>>>>>>> update python version

        html.Div([

            dcc.Graph(id = 'graph2'),

<<<<<<< HEAD
        ], style={'width': '49%', 'display': 'inline-block'}),
=======
        ], style={'width': '33%', 'display': 'inline-block'}),

        html.Div([

            dcc.Graph(id = 'graph3'),

        ], style={'width': '33%', 'display': 'inline-block'}),
>>>>>>> update python version

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
<<<<<<< HEAD
        html.Div([dcc.Graph(id = 'graph_mer2'),], style={'width': '40%','height': '40%', 'display': 'inline-block'})
=======
        html.Div([dcc.Graph(id = 'graph_mer2'),], style={'width': '40%','height': '40%', 'display': 'inline-block'}),

        #display click data

        html.Div([
            dcc.Markdown("""
                **Click Data**

                Click on points in the graph.
            """),
            html.Pre(id='click-data', style=styles['pre']),
        ], className='three columns'),
>>>>>>> update python version




    ],)

<<<<<<< HEAD
=======
#display click data
@callback(
    Output('click-data', 'children'),
    Input('graph3', 'clickData'))
def display_click_data(clickData):
    return json.dumps(clickData, indent=2)



>>>>>>> update python version

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
<<<<<<< HEAD
    Input('tbl', 'active_cell'))
#@callback(Output('tbl_sub', 'children'), Input('tbl', 'active_cell'))
def update_graphs(taxonomy_filter,dataset_filter, row_ids, selected_row_ids,active_cell, cldf = cldf, df = df, counts_data = counts_data):
=======
    Input('tbl', 'active_cell'),
    Input('resolution', 'value'))
#@callback(Output('tbl_sub', 'children'), Input('tbl', 'active_cell'))
def update_graphs(taxonomy_filter,dataset_filter, row_ids, selected_row_ids,active_cell,resolution, cldf = cldf, df = df, counts_data = counts_data ):
>>>>>>> update python version
    #subset to dataset of interest
    #get counts and percents by unique condition, currently library_prep
    
    # if dataset_filter on, get counts by unique condition, else, get counts for every expt
    #if dataset_filter == "All":
    #	counts_data['unique_condition'] = "All"
<<<<<<< HEAD

    counts_data['cluster_counts'] = counts_data.groupby(['unique_condition','cluster_id_label'])['cluster_id_label'].transform('count')
    counts_data['total_counts'] = counts_data.groupby(['unique_condition'])['unique_condition'].transform('count')
    counts_data['cluster_percent'] = (100*counts_data['cluster_counts']) / (counts_data['total_counts'])
    counts_data = counts_data[['cluster_id_label','cluster_counts','unique_condition','cluster_percent']]
    counts_data = counts_data[counts_data['unique_condition'] == dataset_filter]
    counts_data['include'] = counts_data['cluster_id_label']

    #merge with taxonomy, get counts and percents ready
    cldf = cldf.merge(right = counts_data, how = 'left', on = 'cluster_id_label')
    cldf['cluster_counts'] = cldf['cluster_counts'].replace(np.nan, 0)
    cldf['cluster_percent'] = cldf['cluster_percent'].replace(np.nan, 0)
    cldf['cluster_counts'] = cldf['cluster_counts'].astype('int')
    cldf = cldf.drop_duplicates(subset='cl')
    #subset graph
    selected = active_cell if active_cell else "All"
    if selected == "All":
        sub_df_graph = cldf
=======
    counts_data_cluster = counts_data
    counts_data_subclass = counts_data
    cldf_cluster = compute_counts_return_cldf(groupBy = "cluster_id_label", df = counts_data, cldf = cldf, dataset_filter = dataset_filter)
    cldf_select = compute_counts_return_cldf(groupBy = resolution, df = counts_data, cldf = cldf, dataset_filter = dataset_filter)
    #subset graph
    selected = active_cell if active_cell else "All"
    if selected == "All":
        sub_df_graph = cldf_cluster
>>>>>>> update python version
    else:
        groupBy = selected['column_id']
        active_row_id = active_cell['row_id'] if active_cell else None
        selected_value = pd.DataFrame(df).loc[selected['row_id'], groupBy]
<<<<<<< HEAD
        sub_df_graph = cldf[(cldf[groupBy] == selected_value) | (cldf['source'] == "WB")]
=======
        sub_df_graph = cldf_cluster[(cldf_cluster[groupBy] == selected_value) | (cldf_cluster['source'] == "WB")]
>>>>>>> update python version

    if ((taxonomy_filter != "All") ):
    	classes = sub_df_graph[sub_df_graph[taxonomy_filter]==True].class_id_label.unique().tolist()+["0_WB"]
    	sub_df_graph = sub_df_graph[(sub_df_graph[taxonomy_filter]==True) | ( (sub_df_graph['class_id_label'].isin(classes)) & (sub_df_graph['source'] == "WB")  ) ]
    if ( (selected == "All") & (taxonomy_filter == "All") & (dataset_filter == "All")):
    	graph1 = orig_full_graph
    	graph2 = orig_sub_graph
    else:
<<<<<<< HEAD
    	clusters =sub_df_graph.cluster_id_label.unique().tolist()
    	graph1 = build_plotly_taxonomy_full_graph(graph=G,data = cldf,clusters = clusters)
    	graph2 = build_plotly_taxonomy_sub_graph(graph=G,data = cldf,clusters = clusters)
=======
        clusters =sub_df_graph.cluster_id_label.unique().tolist()
        graph1 = build_plotly_taxonomy_full_graph(graph=G,data = cldf_cluster,clusters = clusters)
        graph2 = build_plotly_taxonomy_sub_graph(graph=G,data = cldf_cluster,clusters = clusters)
    #graph3 = build_plotly_bar(data = cldf_select,groupBy = resolution)
>>>>>>> update python version
    return graph1, graph2


@callback(
<<<<<<< HEAD
=======
    Output('graph3','figure'),
    Input('dropdown_taxonomy_filter', 'value'),
    Input('dropdown_dataset_filter', 'value'),
    Input('tbl', 'derived_virtual_row_ids'),
    Input('tbl', 'selected_row_ids'),
    Input('tbl', 'active_cell'),
    Input('resolution', 'value'))
#@callback(Output('tbl_sub', 'children'), Input('tbl', 'active_cell'))
def make_bar_chart(taxonomy_filter,dataset_filter, row_ids, selected_row_ids,active_cell,resolution, cldf = cldf, df = df, counts_data = counts_data ):
    #subset to dataset of interest
    cldf_select = compute_counts_return_cldf(groupBy = resolution, df = counts_data, cldf = cldf, dataset_filter = dataset_filter)
    graph3 = build_plotly_bar(data = cldf_select,groupBy = resolution)
    return graph3

# @callback(
#     Output('graph3','figure'),
#     Input('graph3', 'figure'),
#     Input('graph3', 'clickData'))
# def update_bar_chart(barchart,clickData ):
    
#     graph3 = graph3.update_traces
#     return graph3




@callback(
>>>>>>> update python version
    Output('graph_mer', 'figure'),
    Input('mer_slider', 'value'),
    Input('dropdown_taxonomy_filter', 'value'))


<<<<<<< HEAD
def update_merfish_1(mer_slider, taxonomy_filter,
    cldf = cldf, MER = MER):
    print(mer_slider)
=======


def update_merfish_1(mer_slider, taxonomy_filter,
    cldf = cldf, MER = MER):
>>>>>>> update python version
    if ((taxonomy_filter != "All") ):
        cldf = cldf[(cldf[taxonomy_filter]==True)    ]
        clusters =list(set(cldf.cluster_id_label.unique().tolist()) - set(["0_WB"]))
    else:
        clusters = MER.cluster_id_label.unique().tolist()
    #subset for clusters and section
    MER_sub = MER[(MER['section'] == mer_slider) & MER['cluster_id_label'].isin(clusters)]
<<<<<<< HEAD
    print(MER_sub.head())
=======
>>>>>>> update python version
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
<<<<<<< HEAD
=======
    Input('graph3', 'clickData'),
    Input('resolution', 'value'),
>>>>>>> update python version
    Input('mer_slider', 'value'),
    Input('dropdown_taxonomy_filter', 'value'),
    Input('tbl', 'active_cell'))


<<<<<<< HEAD
def update_merfish_2(mer_slider, taxonomy_filter,active_cell,
    cldf = cldf, df = df, MER = MER):
    print(mer_slider)
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
=======
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
>>>>>>> update python version
    if ((taxonomy_filter != "All") ):
        cldf = cldf[(cldf[taxonomy_filter]==True)    ]
        clusters =list(set(cldf.cluster_id_label.unique().tolist()) - set(["0_WB"]))
    else:
        clusters = MER.cluster_id_label.unique().tolist()
    #subset for clusters and section
    MER_sub = MER[(MER['section'] == mer_slider) & MER['cluster_id_label'].isin(clusters)]
<<<<<<< HEAD
    print(MER_sub.head())
=======
>>>>>>> update python version
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
<<<<<<< HEAD
    print(fig)
=======
>>>>>>> update python version
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

<<<<<<< HEAD

#fig1 = build_plotly_taxonomy_full_graph(graph=G)
#fig1.write_html("AIT21_BG_subset.html")


#build_plotly_taxonomy_sub_graph(graph=G,clusters = clusters)


    #return sub_df_graph.to_json(date_format='iso', orient='split')
    #return data

# @callback(
#     Output('graph1', 'figure'),
#     Input('tbl_sub', 'data'),
#     Input('cldf', 'data'))

# def update_graph(cldf_json, data_json ):
#     cldf = pd.read_json(cldf_json, orient='split')
#     updated_table = pd.read_json(data_json, orient='split')
#     clusters = updated_table.cluster_id_label.unique().tolist()
#     return build_plotly_taxonomy_full_graph(graph=G,data = cldf,clusters = clusters)


# @callback(
#     Output('graph2', 'figure'),
#     Input('tbl_sub', 'data'),
#     Input('cldf','data'))
# def update_graph2(cldf_json,data_json):
#     cldf = pd.read_json(cldf_json, orient='split')
#     updated_table = pd.read_json(data_json, orient='split')
#     clusters = updated_table.cluster_id_label.unique().tolist()
#     return build_plotly_taxonomy_sub_graph(graph=G,data = updated_table,clusters = clusters)






# def apply_taxonomy_filter(tbl_sub  ,taxonomy_filter ): #,
#     if taxonomy_filter == "All":
#         updated_table = tbl_sub
#     else:
#         classes = updated_table[updated_table[taxonomy_filter]==True].class_id_label.unique().tolist()+["0_WB"]
#         updated_table = updated_table[(updated_table[taxonomy_filter]==True) | ( (updated_table['class_id_label'].isin(classes)) & (updated_table['source'] == "WB")  ) ]
    
#     return updated_table

#old dropdowns
# def update_cldf(class_filter  ,taxonomy_filter, subclass_filter ): #,
#     if class_filter == "All":
#         updated_table = cldf
#     else:
#         updated_table = cldf[cldf.class_id_label==class_filter]
    
#     if subclass_filter == "All":
#          updated_table = updated_table
#     else:
#         updated_table = updated_table[updated_table.subclass_id_label.isin([subclass_filter]+["0_WB"])]
#          #subclasses = [subclass_filter]+["0_WB"]
#          #updated_table = updated_table[ ( (updated_table['subclass_id_label'].isin(subclasses)) | (updated_table['source'] == "WB")  ) ]
#     if taxonomy_filter == "All":
#         updated_table = updated_table
#     else:
#         classes = updated_table[updated_table[taxonomy_filter]==True].class_id_label.unique().tolist()+["0_WB"]
#         updated_table = updated_table[(updated_table[taxonomy_filter]==True) | ( (updated_table['class_id_label'].isin(classes)) & (updated_table['source'] == "WB")  ) ]
    
#     return updated_table

#test
#update_cldf(class_filter = "All", taxonomy_filter = "All",subclass_filter = cldf.subclass_id_label[500])




#old with dropdown filters
# @callback(
#     Output('graph1', 'figure'),
#     [Input('dropdown_class_filter', 'value'),
#      Input('dropdown_taxonomy_filter', 'value'),
#      Input('dropdown_subclass_filter', 'value')
#      ] 
# )

# def update_graph(dropdown_class_filter,  dropdown_taxonomy_filter, dropdown_subclass_filter):
#     updated_table = update_cldf(class_filter = dropdown_class_filter,taxonomy_filter = dropdown_taxonomy_filter, subclass_filter = dropdown_subclass_filter) # 
#     return build_plotly_taxonomy_full_graph(graph=G,data = updated_table,clusters = updated_table.cluster_id_label.unique().tolist())


# def update_graph(dropdown_class_filter, dropdown_subclass_filter, dropdown_taxonomy_filter):
#     if dropdown_class_filter == "All":
#         updated_table = cldf
#     else:
#         updated_table = cldf[cldf.class_id_label==dropdown_class_filter]

#     if dropdown_subclass_filter == "All":
#         updated_table = updated_table
#     else:
#         updated_table = updated_table[updated_table.subclass_id_label==dropdown_subclass_filter]
#         #subclasses = [dropdown_subclass_filter]+["0_WB"]
#         #updated_table = updated_table[ ( (updated_table['subclass_id_label'].isin(subclasses)) | (updated_table['source'] == "WB")  ) ]

#     if dropdown_taxonomy_filter == "All":
#         updated_table = updated_table
#     else:
#         classes = updated_table[updated_table[dropdown_taxonomy_filter]==True].class_id_label.unique().tolist()+["0_WB"]
#         updated_table = updated_table[(updated_table[dropdown_taxonomy_filter]==True) | ( (updated_table['class_id_label'].isin(classes)) & (updated_table['source'] == "WB")  ) ]

    
#     return build_plotly_taxonomy_full_graph(graph=G,clusters = updated_table.cluster_id_label.unique().tolist())



# @callback(
#     Output('graph2', 'figure'),
#     [Input('dropdown_class_filter', 'value'),
#      Input('dropdown_taxonomy_filter', 'value'),
#      Input('dropdown_subclass_filter', 'value')] #,Input('dropdown_subclass_filter', 'value')
# )

# def update_graph2(dropdown_class_filter, dropdown_taxonomy_filter, dropdown_subclass_filter):
#     updated_table = update_cldf(class_filter = dropdown_class_filter,taxonomy_filter = dropdown_taxonomy_filter, subclass_filter = dropdown_subclass_filter) 
#     return build_plotly_taxonomy_sub_graph(graph=G,data = updated_table,clusters = updated_table.cluster_id_label.unique().tolist())

# def update_graph2(dropdown_class_filter, dropdown_taxonomy_filter):
#     if dropdown_class_filter == "All":
#         updated_table = cldf
#     else:
#         updated_table = cldf[cldf.class_id_label==dropdown_class_filter]

#     if dropdown_taxonomy_filter == "All":
#         updated_table = updated_table
#     else:
#         classes = updated_table[updated_table[dropdown_taxonomy_filter]==True].class_id_label.unique().tolist()+["0_WB"]
#         updated_table = updated_table[(updated_table[dropdown_taxonomy_filter]==True) | ( (updated_table['class_id_label'].isin(classes)) & (updated_table['source'] == "WB")  ) ]
    
#     return build_plotly_taxonomy_sub_graph(graph=G,clusters = updated_table.cluster_id_label.unique().tolist())
=======
>>>>>>> update python version
