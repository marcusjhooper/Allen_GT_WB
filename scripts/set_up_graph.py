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
import polars
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

#pandas.read_feather('/allen/programs/celltypes/workgroups/rnaseqana')


#_--------------------------------General data setup
#variables
highlight_color = "springgreen"
regular_edge_color = "lightgray"
regular_edge_width = 1
root_value = "gray" # black
pd.options.mode.chained_assignment = None
log_transform = False


#os.chdir('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/WB_hierarchy/')
os.chdir('/home/mh/app/WB_hierarchy/') #local





cldf = pd.read_csv('./allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/Allen_GT_WB_data/AIT21_updated_cldf_for_BG_with_parent.csv')
#### read data and set up colors
clus = pd.read_table('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/Allen_GT_WB_data/WB_colorpal - clusters 230815.tsv')
clus = clus[['cluster_label','cluster_color']]
sub = pd.read_table('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/Allen_GT_WB_data/WB_colorpal - subclasses 230815.tsv')
sub = sub[['subclass_id_label','subclass_color']]
clas = pd.read_table('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/Allen_GT_WB_data/WB_colorpal - classes 230815.tsv')
clas = clas[['class_id_label','class_color']]
sup = pd.read_table('/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/Allen_GT_WB_data/WB_colorpal - supertypes 230815.tsv')
sup = sup[['supertype_id_label','supertype_color']]
root = pd.DataFrame({'id' : ['0_WB'],'node_var' : [root_value]})


clusters = cldf[cldf.BG == "Keep"].cluster_id_label.unique().tolist()

#load pre-downsampled MERFISH data
MER = pd.read_feather('data/app_MERFISH_data.feather')

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
    nodes_to_keep = ["0_WB"] +clusters + cldf[cldf['cluster_id_label'].isin(clusters)].supertype_id_label.tolist() + cldf[cldf['cluster_id_label'].isin(clusters)].class_id_label.tolist() + cldf[cldf['cluster_id_label'].isin(clusters)].subclass_id_label.tolist()
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
    node_attrs.head()

    #set up aesthetics
    node_size = node_attrs.loc[nodes,:]['node_var']
    node_color = node_cols.loc[nodes,:]['node_var']
    node_text = [x+" percent: "+ str(y) for x,y in list(zip(nodes, node_size))]

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
    sub_colors = graph_components['node_color'].loc[list(graph.nodes)]
    sub_sizes = graph_components['node_size'].loc[list(graph.nodes)]

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
                                 size= graph_components['node_size'],
                                 sizeref=2.*max(graph_components['node_size'])/(10.**2),
                                 color = graph_components['node_color'],
                                 line=dict(color='rgb(50,50,50)', width=0.5)
                                 ),
                   text=graph_components['node_text'],
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
cldf = pd.read_csv('./data/AIT21_updated_cldf_for_BG_with_parent.csv')
cldf.index = cldf['cl']
data_file = '/allen/programs/celltypes/workgroups/mct-t200/marcus/VGT/app/Allen_GT_WB_data/test_dataset.csv'
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
