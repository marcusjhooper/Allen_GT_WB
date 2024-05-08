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
import pandas as pd
import seaborn as sns
import networkx
import numpy as np
import pandas as pd
from networkx.drawing.nx_agraph import graphviz_layout
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import networkx as nx
import igraph as ig
import plotly
from plotly.graph_objs import *
from itertools import compress
import s3fs
#pandas.read_feather('/allen/programs/celltypes/workgroups/rnaseqana')


#_--------------------------------General data setup
#variables
highlight_color = "springgreen"
regular_edge_color = "lightgray"
regular_edge_width = 1
root_value = "gray" # black
pd.options.mode.chained_assignment = None
log_transform = False

default_click_data = {
  "points": [
    {
      "curveNumber": 0,
      "pointNumber": 0,
      "pointIndex": 0,
      "x": "046 Vip Gaba",
      "y": 16,
      "label": "046 Vip Gaba",
      "value": 16,
      "bbox": {
        "x0": 2012.57,
        "x1": 2169.1,
        "y0": 776.0699999999999,
        "y1": 776.0699999999999
      }
    }
  ]
}

#os.chdir('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/WB_hierarchy/')

#os.chdir('/home/mh/app/WB_hierarchy/')#local
#os.chdir('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/WB_hierarchy/')


test = False

if test == True:
    os.chdir('/home/mh/app/WB_hierarchy/Allen_GT_WB/Allen_GT_WB') #local
    cldf = pd.read_csv('/home/mh/app/WB_hierarchy_data/AIT21_updated_cldf_for_BG_with_parent.csv')
    clus = pd.read_table('/home/mh/app/WB_hierarchy_data/WB_colorpal - clusters 230815.tsv')
    sub = pd.read_table('/home/mh/app/WB_hierarchy_data/WB_colorpal - subclasses 230815.tsv')
    clas = pd.read_table('/home/mh/app/WB_hierarchy_data/WB_colorpal - classes 230815.tsv')
    sup = pd.read_table('/home/mh/app/WB_hierarchy_data/WB_colorpal - supertypes 230815.tsv')
    MER = pd.read_feather('/home/mh/app/WB_hierarchy_data/app_MERFISH_data.feather')

elif test == False:
    #os.chdir('app/Allen_GT_WB')
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

#### read data and set up colors
root = pd.DataFrame({'id' : ['0_WB'],'node_var' : [root_value]})

#add colors
cldf = cldf.merge(right = clus[['cluster_id_label','cluster_color']], on = 'cluster_id_label')
cldf = cldf.merge(right = sub[['subclass_id_label','subclass_color']], on = 'subclass_id_label')
cldf = cldf.merge(right = sup[['supertype_id_label','supertype_color']], on = 'supertype_id_label')
cldf = cldf.merge(right = clas[['class_id_label','class_color']], on = 'class_id_label')

cldf.index = cldf['cl']

clusters = cldf[cldf.BG == "Keep"].cluster_id_label.unique().tolist()

#load pre-downsampled MERFISH data



MER['point_color'] = 'lightgrey'

# 
#os.chdir('/home/mh/app/WB_hierarchy/') #local



for attr in [clus,sub,clas,sup,root]:
    attr.columns = ['id','node_var']
node_cols = pd.concat([clus,sub,clas,sup,root], ignore_index=True)
node_cols = node_cols.drop_duplicates(subset = 'id')
node_cols.index = node_cols['id']

# set up graph
#wb data separate into connections between different cell type resolutions
sub_cldf = cldf[['subclass_id_label','class_id_label']].drop_duplicates(subset = 'subclass_id_label')
sup_cldf = cldf[['supertype_id_label','subclass_id_label']].drop_duplicates(subset = 'supertype_id_label')
clus_cldf = cldf[['supertype_id_label','cluster_id_label']].drop_duplicates(subset = 'cluster_id_label')
edges = list(zip(sub_cldf['subclass_id_label'],sub_cldf['class_id_label']))
edges_sup = list(zip(sup_cldf['supertype_id_label'],sup_cldf['subclass_id_label']))
edges_clus = list(zip(clus_cldf['supertype_id_label'],clus_cldf['cluster_id_label']))
root = [('0_WB',class_lab) for class_lab in  list(set(cldf['class_id_label'].unique()) - set(["0_WB"]))   ]

#make network graph
G = nx.Graph()
G.add_edges_from(edges)
G.add_edges_from(edges_sup)
G.add_edges_from(root)
G.add_edges_from(edges_clus) #test
nodes = list(G.nodes.keys())
degrees = dict(G.degree) # of lines connected to node - can be used for size instead
pos = graphviz_layout(G, prog="twopi")#plotly graph positions
#_--------------------------------General data setup





#---------------------------------functions
def get_graph_components(graph = G,data = cldf, clusters = clusters,
                        highlight_color = highlight_color,
                        regular_edge_width = regular_edge_width,
                        regular_edge_color = regular_edge_color,
                        attr = 'cluster_percent', pos = pos):
    #keep all clusters, supertypes, subclasses and classes from list of clusters
    nodes_to_keep = ["0_WB"] +clusters + data[data['cluster_id_label'].isin(clusters)].supertype_id_label.tolist() + data[data['cluster_id_label'].isin(clusters)].class_id_label.tolist() + data[data['cluster_id_label'].isin(clusters)].subclass_id_label.tolist()
    nodes_to_keep = pd.Series(nodes_to_keep).unique().tolist()
    
    #get all edges of interest from list of nodes to keep
    edges_to_keep = []
    for node in nodes_to_keep:
        edges_sub = list(G.edges(node))
        edge_indices = list(range(0,len(edges_sub))) 
        edge_filter = [(edges_sub[edge][0] in nodes_to_keep)  &  (edges_sub[edge][1] in nodes_to_keep) for edge in edge_indices  ]
        edge_filter
        edges_to_keep_sub = list(compress(edges_sub,edge_filter))
        edges_to_keep = edges_to_keep+edges_to_keep_sub

    #construct node attr column for size and color
    df = data[['class_id_label',attr]].drop_duplicates(subset = 'class_id_label')
    df2 = data[['subclass_id_label',attr]].drop_duplicates(subset = 'subclass_id_label')
    df3 = data[['supertype_id_label',attr]].drop_duplicates(subset = 'supertype_id_label')
    df4 = data[['cluster_id_label',attr]].drop_duplicates(subset = 'cluster_id_label')
    df5 = pd.DataFrame({'id' : ['root'],'node_var' : [1]})
    for attr in [df,df2,df3,df4,df5]:
        attr.columns = ['id','node_var']
    node_attrs = pd.concat([df,df2,df3,df4,df5], ignore_index=True)
    node_attrs = node_attrs.drop_duplicates(subset = 'id')
    node_attrs.index = node_attrs['id']
    node_attrs.loc[~ node_attrs.index.isin(clusters),'node_var' ] = 0# set non cluster values to 0 for now
    

    #set up aesthetics
    node_size = node_attrs.loc[nodes,:]['node_var']
    node_color = node_cols.loc[nodes,:]['node_var']
    node_labels = node_cols.loc[nodes,:]['id']
    node_labels.index = node_labels
    node_text = node_labels+" percent: "+node_size.astype('str')
    node_text.index = node_labels
    node_text.loc[~ node_text.index.isin(clusters) ] = node_labels.loc[~ node_labels.index.isin(clusters) ] # set non cluster values to 0 for now

    #node_text = [x+" percent: "+ str(y) for x,y in list(zip(node_labels, node_size))]

    log_transform = False
    if(log_transform == True):
        node_attrs['node_var'] = np.log10(node_attrs['node_var']+1)
    scalar = 100
    node_attrs['node_var'] = scalar * (node_attrs['node_var'])
    


    #set up edge colors
    non_highlighted_edges = list(set(G.edges) - set(nodes_to_keep) )
    non_highlighted_edge_col = {edge:regular_edge_color for edge in non_highlighted_edges}
    non_highlighted_edge_col
    edge_col = {edge:highlight_color for edge in edges_to_keep}
    edge_col = edge_col | non_highlighted_edge_col
    plot_components = {"nodes_to_keep": nodes_to_keep,"edges_to_keep":edges_to_keep,"non_highlighted_edges": non_highlighted_edges, "edge_col": edge_col,
                      "node_size": node_size, "node_color": node_color,"node_text": node_text}
    return plot_components


def update_image(clickData):
   color = ['blue'] * 4

   fig = {
       'data': [{
           'type': 'bar',
           'x': [1,2,3,4],
           'y': [10,8,11,7],
           'marker': {
               'color': color
           }
       }],
       'layout': {
           'title': 'click a bar'
       }
   }
   if clickData is not None:
       fig['data'][0]['marker']['color'][clickData['points'][0]['pointNumber']] = 'red'
   return fig


def compute_counts_return_cldf(groupBy, df, cldf, dataset_filter):
    prefix = groupBy.replace("_id_label","")
    if dataset_filter == "All":
        df = df
    else:
        df = df[df['unique_condition'] == dataset_filter]
        
    df[prefix+'_counts'] = df.groupby(['unique_condition',groupBy])[groupBy].transform('count')
    df['total_counts'] = df.groupby(['unique_condition'])['unique_condition'].transform('count')
    df[prefix+'_percent'] = (100*df[prefix+'_counts']) / (df['total_counts'])
    df = df[[groupBy,prefix+'_counts','unique_condition',prefix+'_percent']]
        
    df['include'] = df[groupBy]
    #merge with taxonomy, get counts and percents ready
    cldf_counts = cldf.merge(right = df, how = 'left', on = groupBy)
    cldf_counts[prefix+'_counts'] = cldf_counts[prefix+'_counts'].replace(np.nan, 0)
    cldf_counts[prefix+'_percent'] = cldf_counts[prefix+'_percent'].replace(np.nan, 0)
    cldf_counts[prefix+'_counts'] = cldf_counts[prefix+'_counts'].astype('int')
    cldf_counts = cldf_counts.drop_duplicates(subset=groupBy)
    return cldf_counts






def build_plotly_bar(data,groupBy,clickData = default_click_data, width=800, height=800):

    prefix = groupBy.replace("_id_label","")
    axis=dict(showline=False,zeroline=False,
              showgrid=False,showticklabels=False,
              title='')

    layout=Layout(title= dict(text = "Observed populations",font= dict(size=20)),
        showlegend=False,autosize=False,
        width=width,height=height,yaxis_title=dict(text = "N",font=dict(size=14)),
        hovermode='closest',
        annotations=[
               dict(
               showarrow=False,text='',
                xref='paper',yref='paper',
                x=0,y=-0.1,
                xanchor='left',yanchor='bottom',
                font=dict(size=14)
               )])

    color_col = groupBy.replace('_id_label','_color')
    print(clickData)
    data = data[data[prefix+'_counts'] !=0]
    data.index = data[groupBy]
    try:
        data.loc[clickData["points"][0]["label"], color_col] = '#FF0000' #change clicked bar to red
    except Exception:
        pass
    #data = [go.Bar(x = data[groupBy],y = data[prefix+"_counts"],color= data[color_col])]
    #fig1 = px.bar(data, x=groupBy, y=prefix+"_counts",
    #         hover_data=[groupBy,prefix+"_counts"], marker_color=color_col, height=height, width = width)
    fig1 = go.Figure(data=[go.Bar(
        x = data[groupBy],
        y = data[prefix+"_counts"],
        marker_color=data[color_col]
        )])
    fig1.update_layout(
    autosize=False,
    width=800,
    height=800,
    yaxis_title=dict(text = "N",font=dict(size=20))
    )
    
    #fig1=go.Figure(data=data, layout=layout)
    #fig1.update_layout(yaxis={'visible': False, 'showticklabels': False})
    #fig1.update_layout(xaxis={'visible': False, 'showticklabels': False})
    return fig1

def build_plotly_bar_init(data,groupBy, width=800, height=800,cldf = cldf):

    prefix = groupBy.replace("_id_label","")
    axis=dict(showline=False,zeroline=False,
              showgrid=False,showticklabels=False,
              title='')

    layout=Layout(title= dict(text = "Observed populations",font= dict(size=20)),
        showlegend=False,autosize=False,
        width=width,height=height,yaxis_title=dict(text = "N",font=dict(size=14)),
        hovermode='closest',
        annotations=[
               dict(
               showarrow=False,text='',
                xref='paper',yref='paper',
                x=0,y=-0.1,
                xanchor='left',yanchor='bottom',
                font=dict(size=14)
               )]
                 )


    color_col = groupBy.replace('_id_label','_color')
    data = data[data[prefix+'_counts'] !=0]
    print(data.head())
    data[groupBy]
    data = [go.Bar(x = data[groupBy],y = data[prefix+"_counts"])]
    fig1=go.Figure(data=data, layout=layout)
    fig1.update_layout(
    xaxis = dict(
        tickfont = dict(size=16)),
    yaxis = dict(
        tickfont = dict(size=16)))
    #fig1.update_layout(yaxis={'visible': False, 'showticklabels': False})
    #fig1.update_layout(xaxis={'visible': False, 'showticklabels': False})
    return fig1


def build_plotly_taxonomy_full_graph(graph,data,clusters, width=800, height=800,
                                     regular_edge_color = 'lightgrey',highlight_color = 'springgreen'):
    graph_components = get_graph_components(graph = graph,data = data, clusters = clusters)
    pos = graphviz_layout(graph, prog="twopi")
    nodes = list(graph.nodes)

    Xv=[pos[k][0] for k in nodes] #all edges
    Yv=[pos[k][1] for k in nodes] #all edges
    Xed=[]
    Yed=[]
    
    XHl=[pos[k][0] for k in graph_components['nodes_to_keep']] #highlighted_edges
    YHl=[pos[k][1] for k in graph_components['nodes_to_keep']] #highlighted_edges
    XHled=[]
    YHled=[]
    
    axis=dict(showline=False,zeroline=False,
              showgrid=False,showticklabels=False,
              title='')

    for edge in graph_components['non_highlighted_edges']:
        Xed+=[pos[edge[0]][0],pos[edge[1]][0], None]
        Yed+=[pos[edge[0]][1],pos[edge[1]][1], None]

    for edge in graph_components['edges_to_keep']:
        XHled+=[pos[edge[0]][0],pos[edge[1]][0], None]
        YHled+=[pos[edge[0]][1],pos[edge[1]][1], None]
    layout=Layout(title= "AIT21",font= dict(size=12),
        showlegend=False,autosize=False,
        width=width,height=height,
        hovermode='closest',
        annotations=[
               dict(
               showarrow=False,text='',
                xref='paper',yref='paper',
                x=0,y=-0.1,
                xanchor='left',yanchor='bottom',
                font=dict(size=14)
               )]
                 )
    trace3=Scatter(x=Xed,y=Yed,
                   mode='lines',
                   line=dict(color= regular_edge_color, width=1),
                   hoverinfo='none')
    trace4=Scatter(x=Xv,y=Yv,
                   mode='markers',name='net',
                   #marker_size = graph_components['node_size'],
                   marker=dict(symbol='circle-dot',
                                 size= 10,
                                 color = graph_components['node_color'],
                                 line=dict(color='rgb(50,50,50)', width=0.5)),
                   text=graph_components['node_text'],
                   hoverinfo='text')
    trace5=Scatter(x=XHled,y=YHled,
                   mode='lines',line=dict(color= highlight_color, width=2),
                   hoverinfo='none')
    
    annot="AIT21"
    
    data1=[trace3, trace4,trace5]
    fig1=Figure(data=data1, layout=layout)
    #fig1['layout']['annotations'][0]['text']=annot
    fig1.update_layout(yaxis={'visible': False, 'showticklabels': False})
    fig1.update_layout(xaxis={'visible': False, 'showticklabels': False})
    return fig1




def build_plotly_taxonomy_sub_graph(graph,data,clusters, width=800, height=800,
                                     regular_edge_color = 'darkgray',highlight_color = 'springgreen'):
    graph_components = get_graph_components(graph = graph,data = data, clusters = clusters)
    edges_to_keep = graph_components['edges_to_keep']
    graph = graph.edge_subgraph(edges_to_keep)
    graph_components = get_graph_components(graph = graph,data = data, clusters = clusters)
    sub_colors = graph_components['node_color'].loc[list(graph.nodes)]
    sub_sizes = graph_components['node_size'].loc[list(graph.nodes)]
    node_text = graph_components['node_text'].loc[list(graph.nodes)]
    #sub_colors = graph_components['node_color']
    #sub_sizes = graph_components['node_size']
    
    pos = graphviz_layout(graph, prog="twopi")
    nodes = list(graph.nodes)
    
    Xv=[pos[k][0] for k in nodes] #all edges
    Yv=[pos[k][1] for k in nodes] #all edges
    Xed=[]
    Yed=[]
    
    axis=dict(showline=False, # hide axis line, grid, ticklabels and  title
              zeroline=False,
              showgrid=False,
              showticklabels=False,
              title=''
              )

    #for highlighting - if doesn't work, delete and replace with graph.edges
    #list(edge_col.keys())
    for edge in graph.edges:
        Xed+=[pos[edge[0]][0],pos[edge[1]][0], None]
        Yed+=[pos[edge[0]][1],pos[edge[1]][1], None]

    layout=Layout(title= "AIT21 Subset",
        font= dict(size=12),
        showlegend=False,
        autosize=False,
        width=width,
        height=height,
    
    
        hovermode='closest',
        annotations=[
               dict(
               showarrow=False,
                text='',
                xref='paper',
                yref='paper',
                x=0,
                y=-0.1,
                xanchor='left',
                yanchor='bottom',
                font=dict(
                size=14
                )
                )
            ]
        )
    trace3=Scatter(x=Xed,
                   y=Yed,
                   mode='lines',
                   line=dict(color= regular_edge_color, width=1),
                   hoverinfo='none'
                   )
    trace4=Scatter(x=Xv,
                   y=Yv,
                   mode='markers',
                   name='net',
                   marker=dict(symbol='circle-dot',


                                 size= sub_sizes,
                                 sizeref=2.*max(sub_sizes)/(10.**2),
                                 color = sub_colors,

                                 line=dict(color='rgb(50,50,50)', width=0.5)
                                 ),
                   text=node_text,
                   hoverinfo='text'

                   )

    
    annot="subgraph"
    
    data1=[trace3, trace4]
    fig1=Figure(data=data1, layout=layout)
    #fig1['layout']['annotations'][0]['text']=annot
    fig1.update_layout(yaxis={'visible': False, 'showticklabels': False})
    fig1.update_layout(xaxis={'visible': False, 'showticklabels': False})
    return fig1


   


##orig_full_graph - for resetting
if test == True:
    cldf = pd.read_csv('/home/mh/app/WB_hierarchy_data/AIT21_updated_cldf_for_BG_with_parent.csv')
    data_file = '/home/mh/app/WB_hierarchy_data/test_dataset.csv'
elif test == False:
    cldf = pd.read_csv('s3://mh-allen-gt-wb/app/Allen_GT_WB/data/AIT21_updated_cldf_for_BG_with_parent.csv')
    data_file = 's3://mh-allen-gt-wb/app/Allen_GT_WB/data/test_dataset.csv'

cldf.index = cldf['cl']





counts_data = pd.read_csv(data_file,low_memory = False)

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

counts_data['LabTracksID'] = counts_data['external_donor_name'].astype('str')
#counts_data["unique_condition"] = counts_data[["LabTracksID",'full_genotype','injection_materials',"facs_population_plan",'Region','roi']].agg(" ".join, axis=1)
#temporary
counts_data['unique_condition'] = counts_data['library_prep']



counts_data['cluster_counts'] = counts_data.groupby(['unique_condition','cluster_id_label'])['cluster_id_label'].transform('count')
counts_data['total_counts'] = counts_data.groupby(['unique_condition'])['unique_condition'].transform('count')
counts_data['cluster_percent'] = (100*counts_data['cluster_counts']) / (counts_data['total_counts'])
counts_data = counts_data[['cluster_id_label','cluster_counts','unique_condition','cluster_percent']]
dataset_filter= counts_data['unique_condition'][0]
counts_data = counts_data[counts_data['unique_condition'] == dataset_filter]
counts_data['include'] = counts_data['cluster_id_label']
#merge with taxonomy, get counts and percents ready
cldf = cldf.merge(right = counts_data, how = 'left', on = 'cluster_id_label')
cldf['cluster_counts'] = cldf['cluster_counts'].replace(np.nan, 0)
cldf['cluster_percent'] = cldf['cluster_percent'].replace(np.nan, 0)
cldf['cluster_counts'] = cldf['cluster_counts'].astype('int')
cldf = cldf.drop_duplicates(subset='cl')



clusters =cldf.cluster_id_label.unique().tolist()
orig_full_graph = build_plotly_taxonomy_full_graph(graph=G,data = cldf,clusters = clusters)
orig_sub_graph = build_plotly_taxonomy_sub_graph(graph=G,data = cldf,clusters = clusters)



#to generate data

#MERFISH data - downsample to make plotting faster
# MER = pd.read_csv('/home/mh/app/WB_hierarchy/cell_metadata.csv')
# cell_metadata = pd.read_csv('data/cell_metadata_with_parcellation_annotation.csv')
# MER = MER[MER['average_correlation_score'] > 0.5]
# n_cells = 1000000
# MER = MER.sample(n = n_cells)
# MER['cl'] = MER['cluster_alias']
# MER = MER.merge(right= cldf[['cl','cluster_id_label','subclass_id_label', 'class_id_label','anatomical_annotation']],on= 'cl',how = 'left')
# sections = MER.brain_section_label.unique().tolist()
# MER['section'] = MER['brain_section_label']
# MER = MER.sort_values('z')

# indices = [7,8,9,10,11,20,21,28,30,34,52,58]
# exclude_sections = [sections[i] for i in indices]

# #add colors
# MER = MER.merge(right = clas, on = 'class_id_label')
# MER = MER.merge(right = cell_metadata, how = 'left')

# #exclude bad sections
# MER = MER[~MER['brain_section_label'].isin(exclude_sections)]
# #MER.merge(right = sub, on = 'subclass_id_label')
# #MER.merge(right = clus, on = 'cluster_id_label')
# rep = 0
# for section in MER.brain_section_label.unique().tolist():
#     MER['section'] = MER['section'].replace(section, rep)
#     rep = rep+1

# MER.to_feather('data/app_MERFISH_data.feather')
