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
