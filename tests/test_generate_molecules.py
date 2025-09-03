import os
import time
import tempfile
import typing as t

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from rdkit import Chem
from visual_graph_datasets.visualization.base import create_frameless_figure
from visual_graph_datasets.processing.molecules import MoleculeProcessing
from visual_graph_datasets.visualization.molecules import mol_from_smiles
from visual_graph_datasets.visualization.molecules import visualize_molecular_graph_from_mol

from vgd_counterfactuals.visualization import plot_modification
from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.generate.molecules import DEFAULT_ATOM_VALENCE_MAP
from vgd_counterfactuals.generate.molecules import get_neighborhood
from vgd_counterfactuals.generate.molecules import get_free_valence_map
from vgd_counterfactuals.generate.molecules import get_valid_atom_additions
from vgd_counterfactuals.generate.molecules import get_valid_bond_additions
from vgd_counterfactuals.generate.molecules import get_valid_bond_removals
from vgd_counterfactuals.generate.molecules import is_bridge_head_carbon
from vgd_counterfactuals.generate.molecules import is_nitrogen_nitrogen_sulfur
from .util import ARTIFACTS_PATH


def test_bug_none_molecule_neighbor():
    """
    23.10.23 - There is a bug with the creation of molecular neigbors. It turns out that in some very 
    rare cases the neighborhood creation function can actually produce SMILES codes that are completely 
    invalid such that RDKIT cant even parse them which then leads to an error in processing.process() 
    for downstream tasks.
    """
    # The bug specifically appeared in this SMILES code for one of the generated neighbors.
    smiles = 'Cc1ccc(cc1N\\N=C2\\C(=O)C(=Cc3ccccc23)C(=O)Nc4cc(C)c(NC(=O)C5=Cc6ccccc6C(=N/Nc7c(C)cccc7C(=O)OCCCl)/C5=O)cc4C)C(=O)OCCCl'
    neighbors: t.List[dict] = get_neighborhood(smiles)
    processing = MoleculeProcessing()
    for index, data in enumerate(neighbors):
        value = data['value']
        print(index, value)
        assert value is not None
        assert len(value) > 1
        # The bug would have happened at this line - the processing from a smiles into a molecular graph - for the 
        # ~400th neighbor.
        graph = processing.process(value)
        assert isinstance(graph, dict)


def test_bug_single_atom_counterfactuals():
    """
    20.10.23 - There is a problem with the counterfactuals of small molecules (with 2 atoms) which 
    might generate a single atom neighbor throug edge deletion modifcations. Such single atom neighbors 
    however cause problems in downstream AI applications.
    Therefore we want to make sure that there are no single atom neighbors anymore throug a filter.
    """
    # This operation would cause an error
    neighbors: t.List[dict] = get_neighborhood('CF')
    assert len(neighbors) != 0
    # We want to make sure that none of the generated neighbors is just a single atom!
    for data in neighbors:
        assert len(data['value']) > 1


def test_get_neighborhood():

    smiles = 'CCCC(S)=O'
    mol = Chem.MolFromSmiles(smiles)

    neighbors = get_neighborhood(smiles, fix_protonation=False)
    assert isinstance(neighbors, list)
    assert len(neighbors) != 0

    pdf_path = os.path.join(ARTIFACTS_PATH, 'molecule_get_neighborhood.pdf')
    with PdfPages(pdf_path) as pdf:
        for data in neighbors:
            assert isinstance(data, dict)
            assert 'value' in data
            assert 'mol' in data

            fig, (ax_org, ax_mod) = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))
            fig.suptitle(f'SMILES: {data["value"]}\n'
                         f'modification: {data["type"]}')

            node_positions_org, _ = visualize_molecular_graph_from_mol(ax_org, mol, 1000, 1000)
            plot_modification(ax_org, node_positions_org, *data['org'])

            new_mol = Chem.MolFromSmiles(data['value'])
            node_positions_mod, _ = visualize_molecular_graph_from_mol(ax_mod, new_mol, 1000, 1000)
            plot_modification(ax_mod, node_positions_mod, *data['mod'])

            pdf.savefig(fig)
            plt.close(fig)


def test_get_valid_bond_additions():
    """
    ``get_valid_bond_additions`` is supposed to return a list of dictionaries which describe the valid
    bond addition neighbors and the specific modification that was done to arrive at each new molecule.
    """
    smiles = 'CCCCC'
    mol = Chem.MolFromSmiles(smiles)

    free_valence_map = get_free_valence_map(mol)
    results = get_valid_bond_additions(
        mol,
        free_valence_indices_map=free_valence_map
    )

    # This section will create a PDF file that visualizes the original and each new molecule on each page
    # These visualizations will also point out the exact locations at which the modifications have been
    # done to the molecule.
    pdf_path = os.path.join(ARTIFACTS_PATH, 'get_valid_bond_additions.pdf')
    with PdfPages(pdf_path) as pdf:
        for data in results:
            assert isinstance(data, dict)
            assert 'value' in data
            assert 'mol' in data

            fig, (ax_org, ax_mod) = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))
            node_positions_org, _ = visualize_molecular_graph_from_mol(ax_org, mol, 1000, 1000)
            plot_modification(ax_org, node_positions_org, *data['org'])

            node_positions_mod, _ = visualize_molecular_graph_from_mol(ax_mod, data['mol'], 1000, 1000)
            plot_modification(ax_mod, node_positions_mod, *data['mod'])

            pdf.savefig(fig)
            plt.close(fig)


def test_get_valid_atom_additions():
    """
    ``get_valid_atom_additions`` is supposed to return a list of dictionaries which describe the valid
    bond addition neighbors and the specific modification that was done to arrive at each new molecule.
    """
    smiles = 'CCC'
    mol = Chem.MolFromSmiles(smiles)

    free_valence_map = get_free_valence_map(mol)
    results = get_valid_atom_additions(
        mol,
        atom_valence_map=DEFAULT_ATOM_VALENCE_MAP,
        free_valence_indices_map=free_valence_map
    )

    # This section will create a PDF file that visualizes the original and each new molecule on each page
    # These visualizations will also point out the exact locations at which the modifications have been
    # done to the molecule.
    pdf_path = os.path.join(ARTIFACTS_PATH, 'get_valid_atom_additions.pdf')
    with PdfPages(pdf_path) as pdf:
        for data in results:
            assert isinstance(data, dict)
            assert 'value' in data
            assert 'mol' in data

            fig, (ax_org, ax_mod) = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))
            node_positions_org, _ = visualize_molecular_graph_from_mol(ax_org, mol, 1000, 1000)
            plot_modification(ax_org, node_positions_org, *data['org'])

            node_positions_mod, _ = visualize_molecular_graph_from_mol(ax_mod, data['mol'], 1000, 1000)
            plot_modification(ax_mod, node_positions_mod, *data['mod'])

            pdf.savefig(fig)
            plt.close(fig)


def test_get_valid_bond_removals():
    """
    ``get_valid_bond_removals`` is supposed to return a list of dictionaries which describe the valid
    bond removal neighbors and the specific modification that was done to arrive at each new molecule.
    """
    smiles = 'CCC=CC(CC=O)CCC(CC=C)CCC'
    mol = Chem.MolFromSmiles(smiles)

    results = get_valid_bond_removals(mol)
    assert isinstance(results, list)
    assert len(results) != 0

    # This section will create a PDF file that visualizes the original and each new molecule on each page
    # These visualizations will also point out the exact locations at which the modifications have been
    # done to the molecule.
    pdf_path = os.path.join(ARTIFACTS_PATH, 'get_valid_bond_removals.pdf')
    with PdfPages(pdf_path) as pdf:
        for data in results:
            assert isinstance(data, dict)
            assert 'value' in data
            assert 'mol' in data

            fig, (ax_org, ax_mod) = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))
            fig.suptitle(f'NEW SMILES: {data["value"]}')

            node_positions_org, _ = visualize_molecular_graph_from_mol(ax_org, mol, 1000, 1000)
            plot_modification(ax_org, node_positions_org, *data['org'])

            node_positions_mod, _ = visualize_molecular_graph_from_mol(ax_mod, data['mol'], 1000, 1000)
            plot_modification(ax_mod, node_positions_mod, *data['mod'])

            pdf.savefig(fig)
            plt.close(fig)


def test_bridgehead_carbons_exclusion():
    """
    15.05.23 - Tests a problem which Hunter has noticed about the generation of molecular counterfactuals.
    He pointed out a specific molecular configuration which RDkit technically says is valid solely based
    on the valence rules etc., but which never practically happens in chemistry. Tests a filter which
    detects this kind of pattern and excludes it from the neighborhood of an atom.
    """
    # This molecule contains the structure which is not supposed to happen and should therefore be
    # detected by the filter.
    smiles = 'NCCCC1=CC2=CC=C1C2'
    mol = Chem.MolFromSmiles(smiles)

    processing = MoleculeProcessing()
    fig, node_positions = processing.visualize_as_figure(smiles, width=1000, height=1000)
    fig.savefig(os.path.join(ARTIFACTS_PATH, 'bridgehead_carbons.pdf'))

    match = is_bridge_head_carbon(mol)
    assert match is True

    # Now this is molecule contains only a normal carbon ring and should therefore not be detected by
    # the filter!
    smiles = 'NCCCCC1=CC=CC=C1'
    mol = Chem.MolFromSmiles(smiles)

    match = is_bridge_head_carbon(mol)
    assert match is False


def test_nitrogen_nitrogen_sulfur_exclusion():
    smiles = 'SNNCCCCc1ccccc1'
    mol = Chem.MolFromSmiles(smiles)
    match = is_nitrogen_nitrogen_sulfur(mol)
    assert match is True



