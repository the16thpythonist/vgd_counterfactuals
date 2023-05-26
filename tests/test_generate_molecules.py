import os
import time
import tempfile

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from rdkit import Chem
from visual_graph_datasets.visualization.base import create_frameless_figure
from visual_graph_datasets.processing.molecules import MoleculeProcessing
from visual_graph_datasets.visualization.molecules import mol_from_smiles
from visual_graph_datasets.visualization.molecules import visualize_molecular_graph_from_mol

from vgd_counterfactuals.base import CounterfactualGenerator
from vgd_counterfactuals.generate.molecules import get_neighborhood
from vgd_counterfactuals.generate.molecules import get_free_valence_map
from vgd_counterfactuals.generate.molecules import get_valid_bond_additions
from vgd_counterfactuals.generate.molecules import get_valid_bond_removals
from vgd_counterfactuals.generate.molecules import is_bridge_head_carbon
from vgd_counterfactuals.generate.molecules import is_nitrogen_nitrogen_sulfur
from vgd_counterfactuals.generate.molecules import molecule_differences
from vgd_counterfactuals.generate.molecules import molecule_differences_ecfp
from .util import ARTIFACTS_PATH


def test_molecule_delta():

    # mol_0 = Chem.MolFromSmiles('COc1ccc(C2NC(=S)N3C(c4ccc(OC)cc4)NC(=S)N23)cc1')
    # mol_1 = Chem.MolFromSmiles('COc1ccc(C2NC(S)N3C(c4ccc(OC)cc4)NC(=S)N23)cc1')
    mol_0 = Chem.MolFromSmiles('COc1ccc2nc(N=C(N)N=C(S)Nc3cccc(C)c3)nc(C)c2c1')
    mol_1 = Chem.MolFromSmiles('COc1ccc2nc(N=C(N)N=CNc3cccc(C)c3)nc(C)c2c1')

    fig, rows = plt.subplots(ncols=2, nrows=2, figsize=(20, 20), squeeze=False)
    ax_0, ax_1 = rows[0]
    visualize_molecular_graph_from_mol(ax_0, mol_0, 1000, 1000)
    visualize_molecular_graph_from_mol(ax_1, mol_1, 1000, 1000)

    delta_0, delta_1 = molecule_differences(mol_0, mol_1)

    ax_0, ax_1 = rows[1]
    visualize_molecular_graph_from_mol(ax_0, mol_from_smiles(delta_0), 1000, 1000)
    visualize_molecular_graph_from_mol(ax_1, mol_from_smiles(delta_1), 1000, 1000)

    fig_path = os.path.join(ARTIFACTS_PATH, 'molecule_delta.pdf')
    fig.savefig(fig_path)


def test_get_neighborhood():

    smiles = 'CCCCCC'

    neighbors = get_neighborhood(smiles)
    assert len(neighbors) != 0

    processing = MoleculeProcessing()
    pdf_path = os.path.join(ARTIFACTS_PATH, 'molecule_neighborhood.pdf')
    with PdfPages(pdf_path) as pdf:
        for smiles in neighbors:
            fig, ax = create_frameless_figure(1000, 1000)
            image = processing.visualize(smiles, 1000, 1000)
            ax.imshow(image)
            pdf.savefig(fig)
            plt.close(fig)

    smiles = 'O=O'
    neighbors = get_neighborhood(smiles)
    assert len(neighbors) != 0
    print(neighbors)


def test_get_valid_bond_additions():

    smiles = 'CCCCCC'
    mol = Chem.MolFromSmiles(smiles)

    free_valence_map = get_free_valence_map(mol)
    result = get_valid_bond_additions(
        mol,
        free_valence_indices_map=free_valence_map
    )
    print(result)


def test_get_valid_bond_removals():
    """
    get_valid_bond_removals is supposed to return a set of smiles which each represent a valid bond removal
    action. These actions are either the removal of a single atom or the downgrade of a bond to a different
    type.
    """
    smiles = 'CCC=CC'
    mol = Chem.MolFromSmiles(smiles)

    # For the given smiles we know that there are only three valid node removal actions, which are either
    # the removal of an atom at either end of the chain or the downgrade of the double bond.
    # ['CCCCC', 'C=CCC', 'CC=CC']
    result = get_valid_bond_removals(mol)
    assert len(result) != 0
    assert len(result) == 3


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



