import os
import tempfile

from visual_graph_datasets.processing.molecules import MoleculeProcessing

from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.testing import MockModel
from vgd_counterfactuals.generate.molecules import get_neighborhood
from vgd_counterfactuals.visualization import create_counterfactual_pdf

from .util import ARTIFACTS_PATH


def test_create_counterfactual_pdf():
    # First of all we actually need to create some counterfactuals
    original = 'CCCCCC'

    model = MockModel()
    processing = MoleculeProcessing()
    generator = CounterfactualGenerator(
        model=model,
        processing=processing,
        neighborhood_func=get_neighborhood,
        distance_func=lambda a, b: abs(a - b),
    )

    with tempfile.TemporaryDirectory() as path:
        index_data_map = generator.generate(original, path)
        counterfactual_elements = list(index_data_map.values())
        original_element = generator.create(original, path, index='-1')

        output_path = os.path.join(ARTIFACTS_PATH, 'counterfactual_visualization.pdf')
        create_counterfactual_pdf(
            counterfactual_elements=counterfactual_elements,
            counterfactual_labels=[f'Prediction {element["metadata"]["prediction"]:.3f}'
                                   for element in counterfactual_elements],
            original_element=original_element,
            output_path=output_path
        )

