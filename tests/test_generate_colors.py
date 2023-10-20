import os
import tempfile
import time
import typing as t

from visual_graph_datasets.processing.colors import ColorProcessing
from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.testing import MockModel
from vgd_counterfactuals.generate.colors import get_neighborhood
from vgd_counterfactuals.generate.colors import get_valid_add_edge


def test_bug_get_valid_add_edge_not_actually_adding():
    """
    12.06.23 - There was a bug where the edges were not actually added to the modified graph candidates 
    in the "get_valid_add_edge" function
    """
    cogiles = 'RRRRR'
    
    processing = ColorProcessing()
    graph = processing.process(cogiles)
    num_edges = len(graph['edge_indices'])
    assert num_edges == 8
    
    results: t.List[dict] = get_valid_add_edge(
        graph=graph,
        processing=processing
    )
    for data in results:
        _graph = processing.process(data['value'])
        assert 'edge_indices' in _graph
        # This fails with the bug
        assert len(_graph['edge_indices']) == num_edges + 2



def test_counterfactual_generation_works():
    """
    12.06.23 - Tests if the whole counterfactual setup using the counterfactual generator works 
    properly for the color graphs.
    """
    cogiles = 'MY-1(B-2(H-3(Y)(HY-4(B-5-1(MY-5))(R-6(H)(MC-7-1(GM-2))(M)(G-2)(G)))))(B)-5-7'
    
    model = MockModel()
    processing = ColorProcessing()
    generator = CounterfactualGenerator(
        model=model,
        processing=processing,
        neighborhood_func=get_neighborhood,
        distance_func=lambda a, b: abs(a - b)
    )
    with tempfile.TemporaryDirectory() as path:
        index_data_map = generator.generate(
            original=cogiles,
            path=path,
            k_results=5
        )
        # Surface level testing if the resulting index data map contains the correct number of 
        # elements and for example if the images were at least created properly.
        assert isinstance(index_data_map, dict)
        assert len(index_data_map) == 5
        for key, data in index_data_map.items():
            assert 'image_path' in data
            assert os.path.exists(data['image_path'])
            assert 'metadata' in data
            assert 'graph' in data['metadata']


def test_get_neighborhood_basically_works():
    """
    12.06.23 - get_neighborhood should return all the valid 1-edit neighbors given an initial 
    domain representation (cogiles string) of a color graph.
    """
    cogiles = 'MY-1(B-2(H-3(Y)(HY-4(B-5-1(MY-5))(R-6(H)(MC-7-1(GM-2))(M)(G-2)(G)))))(B)-5-7'
    
    start_time = time.time()
    neighbors = get_neighborhood(cogiles)
    print(f'{len(neighbors)} neighbors in {time.time() - start_time:.2f} seconds')
    
    assert isinstance(neighbors, list)
    assert len(neighbors) != 0
    