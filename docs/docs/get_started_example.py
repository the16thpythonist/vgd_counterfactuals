import tempfile

from visual_graph_datasets.processing.molecules import MoleculeProcessing

from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.testing import MockModel
from vgd_counterfactuals.generate.molecules import get_neighborhood

processing = MoleculeProcessing() 
model = MockModel()

generator = CounterfactualGenerator(
    processing=processing,  # 
    model=model,  # 
    neighborhood_func=get_neighborhood,  # 
    distance_func=lambda orig, mod: abs(orig - mod),  # 
    predict_func=lambda model, graphs: model.predict_graphs(graphs),  #  
)

with tempfile.TemporaryDirectory() as path: # 

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