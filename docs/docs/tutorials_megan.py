import pathlib
import tempfile
from typing import List, Dict

import matplotlib.pyplot as plt
from rich.pretty import pprint
from graph_attention_student.torch.megan import Megan
from visual_graph_datasets.util import dynamic_import
from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.generate.molecules import get_neighborhood

PATH = pathlib.Path(__file__).parent.parent

# ~ loading the pre-trained model & processing
CHECKPOINT_PATH = PATH / 'assets' / 'megan.ckpt'  # Replace with your path
PROCESSING_PATH = PATH / 'assets' / 'process.py'  # Replace with your path

model: Megan = Megan.load(CHECKPOINT_PATH)
processing = dynamic_import(PROCESSING_PATH, 'processing').processing

# ~ defining the necessary functions

def predict_func(model: Megan, graphs: List[dict]):
    results: List[dict] = model.forward_graphs(graphs)
    return results

def distance_func(org: dict, mod: dict):
    distance: float = abs(org['graph_output'] - mod['graph_output'])
    return distance

# ~ setting up the counterfactual generator

generator = CounterfactualGenerator(
    processing=processing,
    model=model,
    neighborhood_func=get_neighborhood,  
    distance_func=distance_func,
    predict_func=predict_func
)

# ~ generating counterfactuals

with tempfile.TemporaryDirectory() as temp_dir:
    
    # Generate counterfactuals for a specific molecule
    index_data_dict: Dict[int, dict] = generator.generate(
        original='CCO',                 # Ethanol
        path=tempfile.gettempdir(),     # Temporary directory for results
        k_results=5                     # Generate 5 counterfactuals
    )

    for index, data in index_data_dict.items():

        image_path = data['image_path']

        graph = data['metadata']['graph']       # Full graph dictionary
        pprint(graph)

        plt.imshow(plt.imread(image_path))
        plt.title(f"Counterfactual {index}")
        plt.axis('off')
        plt.show()
        break