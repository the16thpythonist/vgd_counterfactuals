Changelog
=========

0.1.0 - 01.05.2023
------------------

Initial version


0.1.1 - 10.05.2023
------------------

- Fixed a crucial bug in the generation of chemical counterfactuals where single atom "molecules" were
  allowed and

0.1.2 - 15.05.2023
------------------

- Added the possibility to add additional filters to the generation of a molecular neighborhood which
  exclude certain patterns for the resulting neighbors.
- Added a default molecular filter which excludes bridge head carbon structures from appearing

0.1.3 - 16.05.2023
------------------

- The counterfactual generation now uses multiprocessing for the generation of the graph neighborhood
  itself.
- The predictions for the counterfactuals are now batched to avoid memory overflows when too many
  elements have to be assessed at once.

0.1.4 - 19.05.2023
------------------

- Fixed the functionality for filtering the molecule neighborhood with SMARTS patterns.

0.1.5 - 26.05.2023
------------------

- Added the possibility to not use multiprocessing in the generator because I have the suspicion that
  this causes major problems on the web server

0.1.6 - 26.05.2023
------------------

- Added the optional ``predict_func`` parameter to the generator class which makes it possible to provide
  a custom implementation of how the predictions are obtained from the model.
- Other minor changes

0.2.0 - 01.06.2023
------------------

- Changed the basic interface of how the ``neighborhood_func`` is supposed to work. The input is still
  supposed to be a single original domain representation of the graph, but the output is now a list of
  dictionaries where each dictionary not only encodes the modified graph itself but also additional
  information about the modification such as what the indices of origin for the modification are and what
  type of modification has resulted in that specific graph.
- Also added an additional feature for the molecule ``get_neighborhood`` function where it is possible to
  apply ph range specific protonation to the resulting SMILES to make them more "realistic"
- Added ``dimorphite-dl`` to the dependencies

TODO:

- Fork the dimorphite package and add options to modify the way the SMILES are saved to better work to
  not destroy the node indicing too much.