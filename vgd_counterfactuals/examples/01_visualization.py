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

    # The generator object also offers a "create" method, which will use the
    # given domain-specific graph representation to create new VGD element
    # files in the given "path".
    # The method will return the VGD element dictionary representation for
    # that element.
    # We need to do this separately here for the original element because
    # the VGD element dict of it is required for the visualization and is NOT
    # included when generating the counterfactuals
    original_element = generator.create(e.SMILES, path, index='-1')

    output_path = os.path.join(e.path, 'counterfactuals.pdf')
    create_counterfactual_pdf(
        # This needs to be a list of counterfactual VGD element dictionaries, which
        # is easy since that is exactly what the
        counterfactual_elements=counterfactual_elements,
        # Optionally it is possible to add a custom label to each of the PDF pages
        # belonging to the various counterfactuals by providing a list of strings here.
        # We use that feature to also display the raw prediction for each element.
        counterfactual_labels=[f'Prediction: {element["metadata"]["prediction"]}'
                               for element in counterfactual_elements],
        # Providing the original element is optional. It is also valid to just
        # visualize the counterfactuals, although it is highly recommended to
        # supply the original as well!
        original_element=original_element,
        # it is also possible to add a custom label for the original element
        original_label=f'Prediction: {original_element["metadata"]["prediction"]}',
        output_path=output_path,
    )


experiment.run_if_main()
