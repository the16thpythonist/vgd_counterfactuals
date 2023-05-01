"""
This example shows how the generated counterfactuals can easily be visualized in a single PDF file,
which contains the graph visualization of the original input element as well as all the
counterfactuals.
"""
import os
import pathlib

from pycomex.functional.experiment import Experiment
from pycomex.utils import folder_path, file_namespace
from visual_graph_datasets.processing.molecules import MoleculeProcessing

from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.testing import MockModel
from vgd_counterfactuals.generate.molecules import get_neighborhood
from vgd_counterfactuals.visualization import create_counterfactual_pdf

# The smiles representation of the molecule for which to generate the
# counterfactuals
SMILES = 'CCCCC'
# The number of counterfactuals with the highest distance values to return from the
# generation process
NUM_COUNTERFACTUALS = 10

__DEBUG__ = True


@Experiment(base_path=folder_path(__file__),
            namespace=file_namespace(__file__),
            glob=globals())
def experiment(e: Experiment):

    e.log('setting up the generator...')
    generator = CounterfactualGenerator(
        processing=MoleculeProcessing(),
        model=MockModel(),
        neighborhood_func=get_neighborhood,
        distance_func=lambda a, b: abs(a - b)
    )

    e.log('generating counterfactuals...')
    path = os.path.join(e.path, 'counterfactuals')
    os.mkdir(path)
    index_data_map = generator.generate(
        original=e.SMILES,
        path=path,
        k_results=e.NUM_COUNTERFACTUALS,
    )

    # ~ visualizing the results
    e.log('visualizing counterfactuals...')
    counterfactual_elements = list(index_data_map.values())
    original_element = generator.create(e.SMILES, path, index='-1')

    output_path = os.path.join(e.path, 'counterfactuals.pdf')
    create_counterfactual_pdf(
        counterfactual_elements=counterfactual_elements,
        counterfactual_labels=[f'Prediction: {element["metadata"]["prediction"]}'
                               for element in counterfactual_elements],
        original_element=original_element,
        original_label=f'Prediction: {original_element["metadata"]["prediction"]}',
        output_path=output_path,
    )


experiment.run_if_main()
