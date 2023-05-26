import os
import time
import tempfile

from visual_graph_datasets.processing.molecules import MoleculeProcessing

from vgd_counterfactuals.testing import MockModel
from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.generate.molecules import get_neighborhood


def test_counterfactual_generator_custom_predict():
    """
    26.05.23 - The generator class should now support a custom prediction implementation as parameter,
    where a function can be provided that determines how the prediction results are obtained from the
    model
    """
    original = 'CCC'
    value = 10

    model = MockModel()
    processing = MoleculeProcessing()
    generator = CounterfactualGenerator(
        model=model,
        processing=processing,
        neighborhood_func=get_neighborhood,
        distance_func=lambda a, b: abs(a - b),
        predict_func=lambda modl, graphs: [value for _ in graphs]
    )

    with tempfile.TemporaryDirectory() as path:
        index_data_map = generator.generate(
            original,
            k_results=5,
            k_neighborhood=2,
            path=path,
        )
        for index, data in index_data_map.items():
            # Since we modified the predict function to always return the constant "value" all the
            # predictions correspondingly should now be that value!
            assert data['metadata']['prediction'] == value


def test_counterfactual_generator_multiprocessing_speed():
    """
    16.05.23 - Due to the exponential growth of the number of counterfactuals when increasing the
    neighborhood size, multiprocessing support was added to the generator. The generation of the
    counterfactuals is now done in a parallel manner. Testing if this is actually able to provide
    better speed.
    """
    original = 'NCCNCC(=O)C1=CC=CC=C1C'

    model = MockModel()
    processing = MoleculeProcessing()
    generator = CounterfactualGenerator(
        model=model,
        processing=processing,
        neighborhood_func=get_neighborhood,
        distance_func=lambda a, b: abs(a - b),
    )

    with tempfile.TemporaryDirectory() as path:
        start_time = time.time()
        generator.num_processes = None
        generator.generate(
            original,
            k_results=5,
            k_neighborhood=2,
            path=path,
        )
        duration_1 = time.time() - start_time
        print(f'Duration for 1 processes: {duration_1}')

        start_time = time.time()
        generator.num_processes = 8
        generator.generate(
            original,
            k_results=5,
            k_neighborhood=2,
            path=path,
        )
        duration_8 = time.time() - start_time
        print(f'Duration for 8 processes: {duration_8}')

    # We would at least assume that more processes are actually faster for this case
    assert duration_1 > duration_8


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
        index_data_map = generator.generate(original, path, k_results=5)
        assert len(index_data_map) == 5

