import os
import time
import tempfile

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from rdkit import Chem
from visual_graph_datasets.visualization.base import create_frameless_figure
from visual_graph_datasets.processing.molecules import MoleculeProcessing

from vgd_counterfactuals.generate.molecules import get_neighborhood
from vgd_counterfactuals.generate.molecules import get_free_valence_map
from vgd_counterfactuals.generate.molecules import get_valid_bond_additions
from vgd_counterfactuals.generate.molecules import get_valid_bond_removals
from .util import ARTIFACTS_PATH


def test_get_neighborhood():

    smiles = 'CCCCCC'

    start_time = time.time()
    neighbors = get_neighborhood(smiles)
    end_time = time.time()
    print(f'generated {len(neighbors)} molecular mutations in {end_time - start_time:.3f} seconds')
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


