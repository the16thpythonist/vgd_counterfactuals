# Using MEGAN

This section will showcase a simple application example of how to use the ``vgd_counterfactuals`` library with the [👩‍🏫 Multi-Explanation Graph Attention Network (MEGAN)](https://github.com/aimat-lab/graph_attention_student) model to generate counterfactuals on top of the attributional explanations provided by MEGAN.

## Pre-requisites

**Additional Packages.** This example assumes that you have the following libraries installed, in addition to the ``vgd_counterfactuals`` library:

```bash
pip install graph_attention_student
```

**Pre-trained Model.** This example also assumes that you have access to a pre-trained MEGAN model in the form of a checkpoint file. It's assumed that the model was trained on graph structured obtained through ``MoleculeProcessing`` class which is loaded from a separate ``process.py`` module. If that is not the case&mdash;e.g. different domain&mdash;the corresponding components of the example will have to be adapted accordingly.

## Example

```python title="MEGAN Counterfactuals"

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
    results: List[dict] = model.forward_graphs(graphs)  # (1)
    return results

def distance_func(org: dict, mod: dict):
    distance: float = abs(org['graph_output'] - mod['graph_output'])  # (2)
    return distance

# ~ setting up the counterfactual generator

generator = CounterfactualGenerator(
    processing=processing,
    model=model,
    neighborhood_func=get_neighborhood,  # (3)
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

```

1.  ``Megan`` model class offers the ``forward_graphs`` method which accepts a list of graph 
    dictionaries as an input and returns the corresponding list of model predictions. However, note that these model predictions are again dictionary structures which contain multiple fields&mdash;not just the actual prediction, but also the attributional explanation masks and the graph embedding vector, for example. Returning this list of of dictionaries for the ``predict_func`` means that these dicts will be the input types for the distance function later on.

2.  Since the model prediction results are dictionaries, the input type of the distance function 
    is also dictionaries. For ``Megan`` models the actual model prediction (in this case regression value) is stored in the ``graph_output`` field of the dictionary.

    Note that using a dictionary as the model prediction and consequently as the input type for the distance function allows us more flexibility in the definition of the counterfactual distance. For instance, one could define a distance function which also takes the explanations into account and not just the model prediction.

3.  In this case we can use the *molecule* based neighborhood function which is pre-defined in 
    the library. This function will attempt to do all possible single-edit atom additions, removals and bond additions and removals based on the chemical valence rules using the RDKit library.