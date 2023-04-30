import typing as t

import numpy as np
import visual_graph_datasets.typing as vt


class MockModel:

    def predict_graph(self, graph: vt.GraphDict) -> float:
        node_averages = np.mean(graph['node_attributes'], axis=-1)
        return float(np.sum(node_averages))

    def predict_graphs(self, graphs: t.List[vt.GraphDict]) -> t.List[float]:
        return [self.predict_graph(graph) for graph in graphs]
