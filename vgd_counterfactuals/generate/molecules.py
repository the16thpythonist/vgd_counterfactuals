import itertools
import typing as t
from rdkit import Chem

from vgd_counterfactuals.utils import invert_dict


DEFAULT_ATOM_VALENCE_MAP = {
    'C': 4,
    'N': 5,
    'P': 5,
    'O': 6,
    'F': 7,
    'Cl': 7,
}


def get_neighborhood(smiles: str,
                     atom_valence_map=DEFAULT_ATOM_VALENCE_MAP
                     ) -> t.List[str]:
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

    return list(neighbors_set)


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

                smiles = parts[1]

            results.add(smiles)

    return results

