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


0.3.0 - 12.06.2023
------------------

- changed the ``visual_graph_datasets`` version to 0.13.1 which implements the "unprocessing" feature
- Added the module ``generate.colors`` which now implements counterfactual candidate generation for the 
  the near future
  - Added test cases for the color graph counterfactuals

0.3.1 - 12.06.2023
------------------

- fixed a bug in the generation of color graph counterfactuals where the insertion of an edge added the 
  wrong shape of edge_attributes

0.3.3 - 20.10.2023
------------------

- Added another filter to the generation of the molecular neighborhood which prevents the existance of single 
  atom neighbors as they would cause problems in downstream AI applications.
- Added some more documentation
- updated the readme file.

0.3.4 - 23.10.2023
------------------

- Fixed a rare bug where get_neighborhood for molecules would produce invalid SMILES codes that then led to 
  a processing error down the line.

0.3.5 - 14.11.2023
------------------

- For the generation of molecule neighbors, added parameters which can control what operations to be applied 
  to the original molecule to generate the neighbors
  - Set the default value for the "bond addition" operation to False as that operation often leads to molecules 
    which do not make a lot of realisitic sense.