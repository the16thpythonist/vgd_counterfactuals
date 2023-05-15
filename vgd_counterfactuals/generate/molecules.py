import itertools
import typing as t
from rdkit import Chem

from vgd_counterfactuals.utils import invert_dict


DEFAULT_ATOM_VALENCE_MAP = {
    'C': 4,
    'N': 5,
    'P': 5,
    'O': 6,
    'F': 1,
    'Cl': 7,
    'S': 6,

}


def is_bridge_head_carbon(mol: Chem.Mol, pattern: str = '*1=**2=**=*1*2'):
    """
    Checks if the given molecule ``mol`` is a bridge head carbon structure which apparently does not
    appear in chemistry.
    """
    smarts = Chem.MolFromSmarts(pattern)
    is_match = mol.HasSubstructMatch(smarts)
    return is_match


def get_neighborhood(smiles: str,
                     atom_valence_map=DEFAULT_ATOM_VALENCE_MAP,
                     mol_filters: t.Sequence[t.Callable[[Chem.Mol], bool]] = (
                        is_bridge_head_carbon,
                     ),
                     ) -> t.List[str]:
    """
    Given a ``smiles`` representation of a molecule, this function will return a list of the SMILES
    representations of all the valid molecules within a 1-edit neighborhood. Optionally, a list of boolean
    functions can be provided for ``mol_filters`` to further limit the kinds of molecules included in the
    result.

    :param smiles: The SMILES string representation of
    :param atom_valence_map: A dictionary, whose keys are strings that identify atom types of the SMILES
        notation (O, N, Cl ...) and the values are the integer valence of the corresponding atom. Only
        atoms that are listed in this dict will be considered for edit operations such as adding or
        replacing and atom!
    :param mol_filters: A list of functions which each take a Mol object as input and return a boolean value
        to determine whether that atom should be excluded (True) or not (False).

    :returns: A list of strings
    """
    neighbors_set = set()

    mol = Chem.MolFromSmiles(smiles)
    free_valence_indices_map = get_free_valence_map(mol)

    neighbors_set.update(get_valid_atom_additions(
        mol=mol,
        atom_valence_map=atom_valence_map,
        free_valence_indices_map=free_valence_indices_map,
    ))

    neighbors_set.update(get_valid_bond_additions(
        mol=mol,
        free_valence_indices_map=free_valence_indices_map,
    ))

    neighbors_set.update(get_valid_bond_removals(
        mol=mol
    ))

    # 15.05.23 - All of the above functions will create "valid" molecular SMILES in the sense that RDKit
    # does not tell us that they are completely wrong, but the molecules that are created might still not
    # realistically exist in chemistry, which is why we additionally apply a set of filters that decide
    # for each molecule if it should be included
    neighbors_filtered = []
    for smiles in neighbors_set:
        mol = Chem.MolFromSmiles(smiles)
        for func in mol_filters:
            if not func(mol):
                continue

        neighbors_filtered.append(smiles)

    return neighbors_filtered


def get_free_valence_map(mol: Chem.Mol, max_valence: int = 8) -> dict:
    result = {}
    for i in range(1, max_valence):
        result[i] = [
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetNumImplicitHs() >= i
        ]

    return result


def get_valid_atom_additions(mol: Chem.Mol,
                             atom_valence_map: t.Dict[str, int],
                             free_valence_indices_map: t.Dict[int, t.List[int]]
                             ) -> t.Set[str]:
    """
    Given a molecule ``mol``, this function will return a list of SMILES strings that represent valid
    atom additions to the base molecule.

    :returns: A list of strings.
    """
    results = set()

    value_bond_type_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    for i, bond_type in value_bond_type_map.items():
        # We iterate through all the atom indices
        for atom in free_valence_indices_map[i]:
            for element, valence in atom_valence_map.items():
                new_mol = Chem.RWMol(mol)
                idx = new_mol.AddAtom(Chem.Atom(element))
                new_mol.AddBond(atom, idx, bond_type)
                sanitization_result = Chem.SanitizeMol(new_mol, catchErrors=True)
                if sanitization_result:
                    continue

                smiles = Chem.MolToSmiles(new_mol)
                results.add(smiles)

    return results


def get_valid_bond_additions(mol: Chem.Mol,
                             free_valence_indices_map: t.Dict[int, t.List[int]]
                             ) -> t.Set[str]:
    """
    Given a molecule ``mol``, this function will return a list of SMILES strings which results from valid
    bond additions to that base molecule. Valid bond additions in this case are defined as allowed in
    terms of the atom valences. bond additions may connect two atoms which are not yet directly connected
    or upgrade an existing bond (single to double).

    Also disallowed are bonds between two atoms which are themselves already part of a ring.

    :param mol: The base molecule
    :param free_valence_indices_map: A dict whose keys are integers starting from 0. THe key values
        represent the *free number of valences*. The values are lists of atom indices where each atom in
        the list has to have the corresponding number of free valences given by the dict key.

    :returns: A list of strings
    """
    results = set()

    value_bond_type_map = {
        0: None,
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    bond_type_value_map = invert_dict(value_bond_type_map)

    for valence, atoms in free_valence_indices_map.items():
        for atom1, atom2 in itertools.combinations(atoms, 2):

            # We don't even want to entertain bonds between two atoms which are already part of rings
            if mol.GetAtomWithIdx(atom1).IsInRing() and mol.GetAtomWithIdx(atom2).IsInRing():
                continue

            bond = Chem.Mol(mol).GetBondBetweenAtoms(atom1, atom2)
            new_mol = Chem.RWMol(mol)

            Chem.Kekulize(new_mol, clearAromaticFlags=True)
            if bond is not None:
                bond_type = bond.GetBondType()
                if bond_type not in value_bond_type_map.values():
                    continue

                bond_value = bond_type_value_map[bond_type]
                bond_value += valence
                if bond_value > 3:
                    continue

                index = bond.GetIdx()
                bond.SetBondType(value_bond_type_map[bond_value])
                new_mol.ReplaceBond(index, bond)

            else:
                bond_type = value_bond_type_map[valence]
                new_mol.AddBond(atom1, atom2, bond_type)

            sanitization_result = Chem.SanitizeMol(new_mol, catchErrors=True)
            if sanitization_result:
                continue

            smiles = Chem.MolToSmiles(new_mol)
            results.add(smiles)

    return results


def get_valid_bond_removals(mol: Chem.Mol
                            ) -> t.Set[str]:
    """
    Given a molecule ``mol``, this function will return a list of SMILES strings which represent valid
    bond removals. A bond removal is either a downgrade of an existing bond (e.g. double to single) or the
    removal of a single bond which would mean to disconnect at most a single atom from the rest of the
    molecule!

    :param mol: The base molecule for the removals.

    :returns: A list of strings
    """
    results = set()

    value_bond_type_map = {
        0: None,
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    bond_type_value_map = invert_dict(value_bond_type_map)

    for valence in [1, 2, 3]:
        for bond in mol.GetBonds():
            atom1, atom2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bond = Chem.Mol(mol).GetBondBetweenAtoms(atom1, atom2)

            bond_type = bond.GetBondType()
            if bond_type not in value_bond_type_map.values():
                continue

            new_mol = Chem.RWMol(mol)
            Chem.Kekulize(new_mol, clearAromaticFlags=True)

            bond_value = bond_type_value_map[bond_type]
            bond_value -= valence

            if bond_value == 0:
                new_mol.RemoveBond(atom1, atom2)

            elif bond_value > 0:
                bond.SetBondType(value_bond_type_map[bond_value])
                index = bond.GetIdx()
                new_mol.ReplaceBond(index, bond)

            else:
                continue

            sanitization_result = Chem.SanitizeMol(new_mol, catchErrors=True)
            if sanitization_result:
                continue

            smiles = Chem.MolToSmiles(new_mol)
            if '.' in smiles:
                parts = sorted(smiles.split('.'), key=len)
                if len(parts[0]) > 1:
                    continue

                # 10.05.23 - Not filtering this has caused a bug in the down stream machine learning because
                # we would have essentially allowed "molecules" only consisting of a single atom and those
                # "graphs" would have no edges at all, which is a case that graph neural networks are not
                # prepared to handle.
                if len(Chem.MolFromSmiles(parts[1]).GetAtoms()) < 2:
                    continue

                smiles = parts[1]

            results.add(smiles)

    return results

