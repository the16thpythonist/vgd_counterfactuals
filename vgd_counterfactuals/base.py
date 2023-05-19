import logging
import os
import itertools
import multiprocessing
import typing as t
from copy import deepcopy


import visual_graph_datasets.typing as tv
from visual_graph_datasets.processing.base import ProcessingBase
from visual_graph_datasets.data import load_visual_graph_dataset
from visual_graph_datasets.data import load_visual_graph_element
from visual_graph_datasets.util import Batched

from vgd_counterfactuals.utils import NULL_LOGGER


class CounterfactualGenerator:
    """
    Generic implementation for the generation of counterfactuals.

    An instance of this generator class has to be constructed with the following parameters: ``model`` must
    be a (machine learning) model which implements the ``PredictGraphMixin`` and which predicts the
    property of the graphs in question. ``processing`` must be an instance of ``ProcessingBase``, which
    can process the domain specific representations of a graph (most of the time a string representation)
    into a GraphDict and a visualization image. In that regard it is important that the correct processing
    instance is used which is compatible with the used predictor model. ``distance_func`` should be a
    callable function that determines the distance between two predictions of the model. The counterfactual
    generation will try to maximize that measure of distance. ``neighborhood_func`` has to be a function
    which receives the domain specific representation of an original graph and then returns the entire
    1-edit neighborhood of that graph, which is legible to be a counterfactual.

    :param model: A machine learning model which implements the PredictGraphMixin
    :param processing: An instance of BaseProcessing, which can transform a domain specific
        graph representation into a valid GraphDict for the model to make its prediction from.
    :param neighborhood_func: A function which must accept a single positional parameter which is the
        domain specific representation of a single graph. It should return a list of domain specific
        representations of valid graphs which represent the entire 1-edit neighborhood of the original.
    :param distance_func: A function which must accept two arguments: The original prediction and the
        prediction of a counterfactual. It should return a single numeric value indicating the prediction
        distance between the two elements. This distance measure will be used as the basis for the
        counterfactual choice.
    """

    DEFAULT_IMAGE_WIDTH = 1000
    DEFAULT_IMAGE_HEIGHT = 1000

    def __init__(self,
                 model,
                 processing: ProcessingBase,
                 neighborhood_func: t.Callable,
                 distance_func: t.Callable[[t.Any, t.Any], float],
                 logger: logging.Logger = NULL_LOGGER,
                 num_processes: int = 2,
                 batch_size: int = 5_000,
                 ):
        self.model = model
        self.processing = processing
        self.neighborhood_func = neighborhood_func
        self.distance_func = distance_func
        self.logger = logger
        self.num_processes = num_processes
        self.batch_size = batch_size

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
        # 16.05.23 - For a larger neighborhood size the number of counterfactuals increases exponentially
        # and generating them takes forever. That is why we use multiprocessing to do the generation in
        # parallel and hopefully be a bit faster.
        with multiprocessing.Pool(processes=self.num_processes) as pool:
            neighbors_set = set([original])
            for i in range(k_neighborhood):
                neighbors_set = set(itertools.chain(*pool.map(self.neighborhood_func, neighbors_set)))

        # This will be a list of the domain-spec. representations of all the generated neighbors of
        # the original graph.
        neighbors: t.List[tv.DomainRepr] = list(neighbors_set)
        graphs = [self.processing.process(value) for value in neighbors]

        # 16.05.23 - I have noticed that for larger neighborhood sizes there can be A LOT of counterfactuals
        # so many that a model cannot possibly predict them all at once without causing a memory overflow
        # for the resulting vector. This is why we need to process this in a batched manner now!
        # NOTE: One might think that this could be parallelized as well, like the counterfactual generation
        # itself, but that would actually not work because we cant serialize the model to send it to another
        # process and loading the model from memory takes so long that it would not make sense to do that
        # in each process either.
        predictions = []
        for graph_batch in Batched(graphs, batch_size=self.batch_size):
            predictions += self.model.predict_graphs(graph_batch)

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





