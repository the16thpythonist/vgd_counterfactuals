import logging
import os
import typing as t
from copy import deepcopy

from visual_graph_datasets.processing.base import ProcessingBase
from visual_graph_datasets.data import load_visual_graph_dataset

from vgd_counterfactuals.utils import NULL_LOGGER


class CounterfactualGenerator:
    """

    """

    def __init__(self,
                 model,
                 processing: ProcessingBase,
                 neighborhood_func: t.Callable,
                 distance_func: t.Callable[[t.Any, t.Any], float],
                 logger: logging.Logger = NULL_LOGGER,
                 ):
        self.model = model
        self.processing = processing
        self.neighborhood_func = neighborhood_func
        self.distance_func = distance_func
        self.logger = logger

    def generate(self,
                 original: t.Union[str],
                 path: str,
                 k_results: int = 5,
                 k_neighborhood: int = 1,
                 image_width: int = 1000,
                 image_height: int = 1000,
                 ):
        """

        """
        graph = self.processing.process(original)

        original_prediction = self.model.predict_graph(graph)

        neighbors_set = set([original])
        for i in range(k_neighborhood):
            for value in deepcopy(neighbors_set):
                neighbors_set.update(set(self.neighborhood_func(value)))

        values = list(neighbors_set)
        graphs = [self.processing.process(value) for value in values]
        predictions = self.model.predict_graphs(graphs)
        distances = [self.distance_func(original_prediction, pred) for pred in predictions]

        sorted_results = sorted(
            zip(values, graphs, distances),
            key=lambda tupl: tupl[2],
            reverse=True
        )
        num = min(len(sorted_results), k_results)
        top_results = sorted_results[:num]

        # For these top results we now want to create a visual graph dataset folder so that they can be
        # visualized and processed further
        for index, (value, graph, distance) in enumerate(top_results):
            self.processing.create(
                value,
                index=str(index),
                name=value,
                output_path=path,
                width=image_width,
                height=image_height,
            )

        metadata_map, index_data_map = load_visual_graph_dataset(
            path=path,
            logger=self.logger,
            log_step=10,
        )
        return index_data_map

