import logging
import os
import typing as t
from copy import deepcopy

import visual_graph_datasets.typing as tv
from visual_graph_datasets.processing.base import ProcessingBase
from visual_graph_datasets.data import load_visual_graph_dataset
from visual_graph_datasets.data import load_visual_graph_element

from vgd_counterfactuals.utils import NULL_LOGGER


class CounterfactualGenerator:
    """

    """

    DEFAULT_IMAGE_WIDTH = 1000
    DEFAULT_IMAGE_HEIGHT = 1000

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
            zip(distances, values, predictions, graphs),
            key=lambda tupl: tupl[0],
            reverse=True
        )
        num = min(len(sorted_results), k_results)
        top_results = sorted_results[:num]

        # For these top results we now want to create a visual graph dataset folder so that they can be
        # visualized and processed further
        for index, (distance, value, prediction, graph) in enumerate(top_results):
            self.processing.create(
                value,
                index=str(index),
                name=value,
                output_path=path,
                width=image_width,
                height=image_height,
                additional_metadata={
                    'distance': distance,
                    'prediction': prediction
                }
            )

        metadata_map, index_data_map = load_visual_graph_dataset(
            path=path,
            logger=self.logger,
            log_step=10,
        )
        return index_data_map

    def create(self,
               value: t.Any,
               path: str,
               index: str,
               image_width=DEFAULT_IMAGE_WIDTH,
               image_height=DEFAULT_IMAGE_HEIGHT,
               ) -> dict:

        self.processing.create(
            value,
            index=str(index),
            name=value,
            output_path=path,
            width=image_width,
            height=image_height,
        )

        element = load_visual_graph_element(
            path=path,
            name=index
        )
        element["metadata"]["prediction"] = self.model.predict_graph(element["metadata"]["graph"])

        return element





