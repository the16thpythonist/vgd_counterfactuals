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
                 original: tv.DomainRepr,
                 path: str,
                 k_results: int = 5,
                 k_neighborhood: int = 1,
                 image_width: int = DEFAULT_IMAGE_WIDTH,
                 image_height: int = DEFAULT_IMAGE_HEIGHT,
                 ) -> tv.VisGraphIndexDict:
        """
        Creates the ``k_results`` top counterfactuals for the prediction of the ``original`` element. The
        resulting top counterfactuals will be created as VGD elements into the folder given by the
        absolute ``path``.

        :param original: The domain-specific representation of the graph, for whose prediction results the
            counterfactuals should be generated.
        :param path: The absolute path to the existing folder, into which the VGD element files for all
            the counterfactuals are to be saved into.
        :param k_results: The integer number of the top counterfactuals to be returned. The VGD element
            files will only be created for these top results!
        :param k_neighborhood: The integer number of how large a graph edit neighborhood to explore. Default
            value is 1, in which case only the immediate neighbors (single edits) will be explored for
            possible counterfactuals. If this value would be 2, for example, the single-edit neighborhood
            of each of those neighbors would be explored as well.
            NOTE: Computational load rises exponentially, increasing this value is not encouraged.
        :param image_height: The integer width in pixels for the visualization images.
        :param image_width: The integer height in pixels for the visualization images.

        :returns: The visual graph dataset index_data_map for the created VGD elements of all the
            generated counterfactuals.
        """
        # To be able to generate counterfactuals at all we first need to start from
        # the prediction of the original element. To get this prediction we first need
        # to convert the domain-specific representation into a GraphDict
        graph = self.processing.process(original)
        original_prediction = self.model.predict_graph(graph)

        # Based on the original we need ALL possible neighbors of this in a k-neighborhood
        neighbors_set = set([original])
        for i in range(k_neighborhood):
            for value in deepcopy(neighbors_set):
                neighbors_set.update(set(self.neighborhood_func(value)))

        # This will be a list of the domain-spec. representations of all the generated neighbors of
        # the original graph.
        neighbors: t.List[tv.DomainRepr] = list(neighbors_set)
        graphs = [self.processing.process(value) for value in neighbors]
        predictions = self.model.predict_graphs(graphs)
        distances = [self.distance_func(original_prediction, pred) for pred in predictions]

        sorted_results = sorted(
            zip(distances, neighbors, predictions, graphs),
            key=lambda tupl: tupl[0],  # sort by distance
            reverse=True,  # sort in descending order
        )
        # If there are less neighbors in total than the desired number of results, we can only
        # provide as many as there are.
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
               value: tv.DomainRepr,
               path: str,
               index: str,
               image_width=DEFAULT_IMAGE_WIDTH,
               image_height=DEFAULT_IMAGE_HEIGHT,
               ) -> tv.VgdElementDict:
        """
        Will create and return a VGD element dict for the given graph with the domain-specific
        representation ``value``. The VGD element files (JSON, PNG) will be saved into the existing folder
        identified by the absolute ``path``. The files will be called according to the given ``index``.

        NOTE: The created visual graph element will be annotated with the subject model's prediction, saved
        as the field "prediction" in the "metadata" dictionary

        :param value: The domain specific representation of the graph which to generate as a VGD element
        :param path: The absolute path to an existing folder, into which the VGD element files
            will be saved into
        :param index: A string of an integer that defines the index of the element.
        :param image_height: The integer width in pixels for the visualization images.
        :param image_width: The integer height in pixels for the visualization images.

        :returns: A VGD element dict representing the given graph
        """
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





