import os
import time
import tempfile

from visual_graph_datasets.processing.molecules import MoleculeProcessing

from vgd_counterfactuals.testing import MockModel
from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.generate.molecules import get_neighborhood


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
        generator.num_processes = 2
        generator.generate(
            original,
            k_results=5,
            k_neighborhood=2,
            path=path,
        )
        duration_2 = time.time() - start_time
        print(f'Duration for 2 processes: {duration_2}')

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
    assert duration_2 > duration_8


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
