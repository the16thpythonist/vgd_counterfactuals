import os
import tempfile

from visual_graph_datasets.processing.molecules import MoleculeProcessing

from vgd_counterfactuals.testing import MockModel
from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.generate.molecules import get_neighborhood


def test_counter_factual_generator_basically_works():
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
        print(len(index_data_map))
        print(index_data_map)
