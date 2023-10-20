import typing as t 

import numpy as np
import visual_graph_datasets.typing as tv
from visual_graph_datasets.graph import copy_graph_dict
from visual_graph_datasets.graph import graph_edge_set, graph_node_adjacency
from visual_graph_datasets.graph import graph_add_edge, graph_remove_edge
from visual_graph_datasets.generation.colors import DEFAULT_COLORS
from visual_graph_datasets.processing.colors import ColorProcessing



def get_neighborhood(value: str,
                     colors: t.List[str] = DEFAULT_COLORS,
                     processing: ColorProcessing = ColorProcessing(),
                     ) -> t.List[dict]:
    graph = processing.process(value)
    
    results = []
    results += get_valid_node_replace(
        graph=graph,
        colors=colors,
        processing=processing
    )
    results += get_valid_remove_edge(
        graph=graph,
        processing=processing,
    )
    results += get_valid_add_edge(
        graph=graph,
        processing=processing,
    )
    
    return results



def get_valid_node_replace(graph: tv.GraphDict,
                           colors: t.List[t.Any],
                           processing: ColorProcessing,
                           ) -> t.List[dict]:
    
    results = []
    for i in graph['node_indices']:
        new_graph = copy_graph_dict(graph)
        attributes = graph['node_attributes'][i]
        replacements = [color for color in colors if not np.isclose(color, attributes).all()]
        for color in replacements:
            new_graph['node_attributes'][i] = color
            results.append({
                'type': 'replace_node',
                'value': processing.unprocess(new_graph),
                'org': (i, i),
                'mod': (i, i),
            })
            
    return results


def get_valid_remove_edge(graph: tv.GraphDict,
                          processing: ColorProcessing
                          ) -> t.List[dict]:
    results = []
    for i, j in graph_edge_set(graph):
        new_graph = copy_graph_dict(graph)
        new_graph = graph_remove_edge(new_graph, i, j, directed=False)
        value = processing.unprocess(new_graph)
        
        # Here we have two conditions in which case we accept the modified graph as a 
        # valid neighbor:
        # - if the graph is still fully connected, which can be easily checked form the COGILES 
        #   string as the colon will only be used for disconnected partial graphs
        # - if the graph is disconnected and one of the disconnected parts consists of only a single 
        #   node. In that case this is essentially the same as a single node removal which will 
        #   be done as well.
        if '.' not in value or 1 in [len(part) for part in value.split('.')]:
            results.append({
                'type': 'remove_edge',
                'value': value,
                'org': (i, j),
                'mod': (i, j),
            })
        
    return results


def get_valid_add_edge(graph: tv.GraphDict,
                       processing: ColorProcessing
                       ) -> t.List[dict]:
    results = []
    node_adjancency = graph_node_adjacency(graph)
    for i in graph['node_indices']:
        for j in graph['node_indices']:
            if i > j and not node_adjancency[i, j]:
                new_graph = copy_graph_dict(graph)
                graph_add_edge(
                    graph=new_graph, 
                    node_index_1=i, 
                    node_index_2=j, 
                    directed=False,
                    attributes={
                        'edge_attributes': [1]
                    }
                )
                
                results.append({
                    'type': 'add_edge',
                    'value': processing.unprocess(new_graph),
                    'org': (i, j),
                    'mod': (i, j)
                })

    return results