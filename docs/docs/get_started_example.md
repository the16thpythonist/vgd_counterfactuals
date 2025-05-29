# Introduction by Example

This section will provide a simple mock example of how to use the VGD Counterfactuals library to generate counterfactuals for a given graph representation.

```python title="Basic Example"
import tempfile

from visual_graph_datasets.processing.molecules import MoleculeProcessing

from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.testing import MockModel
from vgd_counterfactuals.generate.molecules import get_neighborhood

processing = MoleculeProcessing() 
model = MockModel()

generator = CounterfactualGenerator(
    processing=processing,  # (1)
    model=model,  # (2)
    neighborhood_func=get_neighborhood,  # (3)
    distance_func=lambda orig, mod: abs(orig - mod),  # (4)
    predict_func=lambda model, graphs: model.predict_graphs(graphs),  # (5) 
)

with tempfile.TemporaryDirectory() as path: # (6)

    # The "generate" function will create all the possible neighbors of the
    # given "original" element, then query the model for to predict the
    # output for each of them, and sort them by their distance to the original.
    # The top k elements will be turned into a temporary visual graph dataset
    # within the given folder "path".
    index_data_map = generator.generate(
        original='CCCCCC',
        # Path to the folder into which to save the vgd element files
        path=path,
        # The number of counterfactuals to be returned.
        # Elements will be sorted by their distance.
        k_results=10,
    )

    # The keys of the resulting dict are the integer indices and the values
    # are dicts themselves which describe the corresponding vgd elements.
    # These dicts contain for example the absolute path to the PNG file,
    # the full graph representation and additional metadata.
    print(f'generated {len(index_data_map)} counterfactuals:')
    for index, data in index_data_map.items():
        print(f' * {data["metadata"]["name"]} '
              f' - distance: {data["metadata"]["distance"]:.2f}')

```

1.  The domain specific ``Processing`` instance from the VGD libarary which implements the conversion
    of the string domain representation into the graph dictionary format. This graph dictionary format 
    is then what is usually passed to the actual model for the prediction.

    In this example, we use the ``MoleculeProcessing`` class which is able to take a SMILES string 
    representation of a molecule and convert it into a graph dictionary format.

2.  The model which will be used to obtain the predictions for the original and the modified graphs.
    This model may in principle by any kind of model.

3.  The neighborhood function which defines how to find *similar* and *valid* graphs in the domain. 
    This function will be called with the original graph and should return a list of similar graphs which are exactly one modification away from the original graph&mdash;what exactly one modification means, however, relies on the specific domain.

    A valid neighborhood function ``get_neighborhood(value: str) -> List[str]`` accepts a single string domain representation and returns a list of string which are again valid domain representations.

4.  The distance function which defines how to measure the difference between the predictions 
    of the original and modified graphs. 

    A valid distance function ``distance_func(original: float | np.ndarray, modified: float | np.ndarray) -> float`` accepts the predictions of the original and modified graphs and returns a numeric value which quantifies the difference between them. Note that the predictions can either be single values (regression) or vectors (classification) depending on the model used.

5.  The predict function which defines how to obtain the predictions for a given *list* of graphs. 
    This function will be called with the model and a list of graph dictionaries and should return a list of predictions for each graph.

    A valid predict function ``predict_func(model: Any, graphs: List[Dict]) -> List[float | np.ndarray]`` accepts the model and a list of graph dictionaries and returns a list of predictions for each graph. Note that the graphs which are passed here are the graph dictionaries *not* the string representations. Also note, that whatever kind of data structure (float, np array, etc) the model returns here as the result of the prediction will be the arguments passed to the distance function.

6.  For the sake of the example we will create a temporary directory to store the generated files 
    in. In a real application, it might make sense to store the results to a more permanent location.

    In essence the counterfactual generation will store the counterfactuals as a visual graph dataset
    which consists of 2 files per graph: A visualization PNG and a JSON file containing the graph 
    representation and metadata.


## The ``CounterfactualGenerator`` Class

The ``CounterfactualGenerator`` class is the central component of the library to orchestrate the counterfactual generation process. An instance of this class will have to be initialized once for the generation of counterfactuals for a specific domain and based on a specific model.

The counterfactual generator requires a series of arguments to be passed to its constructor which define the various steps of the generation process, such as how the domain representations are processed into graph structures, which model to use for the predictions, how to find similar graphs in the domain and how to measure the difference between the predictions of the original and modified graphs.

## The ``generate`` Method

The generate method of the counterfactual generator will then actually generate the counterfactuals for an original graph&mdash;given in its string domain representation. The method will return the ``k_results`` top ranked counterfactuals based on the distance function defined in the constructor.

Internally, the method will do the following:

1. Convert the original string representatio into a graph and query the model to obtain the original prediction.

2. Use the neighborhood function to find all valid modifications of the original graph.

3. Query the model for the predictions of all modified graphs.

4. Calculate the distance between the original prediction and each modified prediction using the distance function.

5. Sort the modified graphs by their distance to the original graph.

6. Return the top `k_results` modified graphs as a dictionary mapping indices to graph data.

Specifically, the return value of the `generate` method is a dictionary in the visual graph dataset format. Moreover, the visual graph dataset will also directly be stored persistently to the given ``path`` directory, including the PNG visualizations and the JSON graph representations.