# Core Concepts

This section will introduce the core concepts of the ``vgd_counterfactuals`` library including the idea of **⚒️ counterfactual mining**, **🗺️ graph domains** and **🧩 modularization**.

## ⚒️ Counterfactual Mining

### What are Counterfactuals?

The core idea of counterfactuals&mdash;independent of the data structure for now&mdash;is to find elements which represent only small modifications regarding the original element, but cause a large deviation in the model's property prediction.

Mathematically, this can be expressed as an optimization problem in which we want to simultaneously *minimize the distance* to the original element and *maximize the difference* to the original model prediction.

### Counterfactuals in other Domains

In other application domains, the previously described optimization problem is non-trivial to solve, as the space of possible modifications is extremely large and only a small subset of modifications leads to a valid element.

In image processing, for example, pixels are represented as continuous values, and the space of possible modifications is infinite. Furthermore, a counterfactual modification should itself represent a valid image in the context of the underlying task to be meaningful. The modification of a handwritten digit should result in something that still looks like a digit, and&mdash;crucially&mdash;not a random noise pattern. Due to these constraints, these domains have to use more sophisticated methods to solve the non-trivial optimization problem of finding counterfactual samples. In image generation, for instance one possible approach would be to use gradient descent to iteratively modify the original image using a joint objective function that combines the distance to the original image and the difference in model prediction.

In natural language processing, the space of possible modifications is at least discrete, but still very large since there is an extremely large number of words which could be inserted, removed or replaced at any position in a given text.

### Graph Counterfactuals

In contrast to these previously mentioned domains, graph data structures are much easier to handle. For one thing, graphs are discrete structures for which nodes and edges can either exist or not exist. Additionally, the *vocabulary* of different node and edge types is usually much smaller than with natural language, for instance.

Therefore, one possible approach to generate counterfactuals&mdash;infeasible in other domains&mdash;is to generate all possible modifications of the original element, obtain the model's prediction for each of these modifications and then filter the modifications that lead to the largest deviation in the model's prediction. This is the basic approach taken by the ``vgd_counterfactuals`` library.

## 🗺️ Graph Domains

One important consideration for the generation of counterfactuals as described in the previous paragraph is that *not all graphs are the same*. The rules of what constitutes a valid graph structure&mdash;and what modifications are allowed&mdash;depend on the *domain* of where the graph originates from.

For example, a molecular graph is a specific type of graph where nodes represent atoms and edges represent chemical bonds. In this case, there are very specific rules about which kinds of connections between atoms are allowed to constitute a *valid* molecule. Such a molecular graph needs completely different rules than an electrical circuit graph, where nodes represents connections and edges represent electrical components.

Ultimately, the 🗺️ graph domain heavily influences how counterfactuals need to be generated, which is why they are a core concept of the ``vgd_counterfactuals`` library. Currently, the following domains are supported:

- **Simple Color Graphs**: For testing and demonstration purposes, where nodes are colored and edges connect nodes of different colors.
- **Molecular Graphs**: Using RDKit and SMILES notation, where nodes represent atoms and edges represent chemical bonds.

### 📝 Domain Representation

The ``vgd_counterfactuals`` library builds on top of the **Visual Graph Datasets (VGD)** library for it's graph representation. Consequently, it borrows the VGD concept of the *📝 domain representation*. In the VGD library, each graph can be represented either generically as a graph structure with numeric node and edge features (ideal for machine learning applications) or in a domain-specific *string* format.

An example of such a domain-specific string representation is the SMILES (Simplified Molecular Input Line Entry System) notation for molecular graphs, which is a widely used format to represent chemical structures as strings. The caffeine molecule, for example, can be represented as the SMILES string `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`. This format is very compact and allows for easy sharing and storage of molecular graphs. More importantly, this exchange format is the foundation for many domain-specific tools such as the cheminformatics toolkit RDKit.

## 🧩 Modularization

The goal of the ``vgd_counterfactuals`` library is to be *customizable* and to provide a **modular framework** for counterfactual generation. This means that individual components such as the model, the neighborhood function etc. can be swapped out and/or be replaced with custom implementations.

The library is designed around the following core components which can be assembled to interact with each other:

- **Neighborhood Function.** Given the *📝 domain representation* of a graph structure, this function should return the domain representations of it's immediate neighbors&mdash;meaning all valid graphs with an edit distance of 1. The neighborhood function is domain-specific and will have to be implemented or customized for each new graph domain.
- **Counterfactual Distance Function.** Given an original graph and a modified graph, as well as the corresponding model predictions, this function is meant to return a numeric value that represents a measure of how different the two *predictions* are. This function will have to be supplied by the user as it depends on the specific circumstances. For instance, different distance functions will be necessary for regression tasks than for classification tasks. Another important consideration is that sometimes only one direction of the distance matters, e.g. increasing/decreasing the value in a regression task.
- **Predictor Model.** In principle, the library is meant to be agnostic to the specific model implementation, meaning that it should be possible to generate counterfactuals based on *any* model. However, to make a given model work with the library, one needs to supply an adapter function which accepts the domain representation of a graph, queries the model and returns the model's prediction.