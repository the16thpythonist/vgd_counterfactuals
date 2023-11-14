import itertools
import typing as t
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem

from dimorphite_dl import DimorphiteDL
from vgd_counterfactuals.utils import invert_dict


DEFAULT_ATOM_VALENCE_MAP = {
    'C': 4,
    'N': 5,
    # 'P': 5,
    'O': 6,
    'F': 1,
    'Cl': 7,
    'S': 6,
}


# Molecule Filters
# ----------------
# There are some situations where the counterfactual generation might generate certain kinds of molecules 
# that might be in the strictest sense chemically valid, but which still do not make realistic sense. These 
# cases will be filtered from the results by applying a set of filter rules to the final generated counterfactual 
# list. These filters are implemented as the possible callback functions below.

def is_bridge_head_carbon(mol: Chem.Mol, pattern: str = '*1=**2=**=*1*2'):
    """
    Checks if the given molecule ``mol`` is a bridge head carbon structure which apparently does not
    appear in chemistry.
    """
    smarts = Chem.MolFromSmarts(pattern)
    is_match = mol.HasSubstructMatch(smarts)
    return is_match


def is_nitrogen_nitrogen_sulfur(mol: Chem.Mol, pattern: str = 'SNN'):
    """
    This checks if a molecules contains two nitrogen atoms connected to a sulfur atom in a chain. This is a 
    configuration which does not really appear in nature.
    """
    smarts = Chem.MolFromSmarts(pattern)
    is_match = mol.HasSubstructMatch(smarts)
    return is_match


def is_single_atom(mol: Chem.Mol):
    """
    Checks if the given molecule is essentially just a single atom.
    
    We want to filter this trivial case, because it will likely cause problems for the downstream AI applications, 
    since it is not possible to extract a valid graph representation that consists of just a single node.
    """
    num_atoms = len(mol.GetAtoms())
    return num_atoms < 2


# Protonation
# -----------
# Oftentimes, molecules exist in various different protonation states. That means that for a given molecule 
# not all of the are always properly bonded to a hydrogen atom, but that they might exist in their charged 
# form. For moderately large molecules, there exist many different permutations of which atoms are charged 
# and which are properly protonated. Which of those states are possible at all depend on the pH of the 
# environment.
# In certain applications is might make sense to consider these different protonation states when generating 
# the local neighborhood, which is implemented here.


def fix_protonation_dimorphite(smiles_list: t.List[str],
                               min_ph: float,
                               max_ph: float,
                               max_variants: float = 10,
                               ) -> t.List[str]:
    """
    This method generates the different protonation states in a given ph range ``min_ph`` and ``max_ph`` for 
    each molecule represented as a SMILES string in the given ``smiles_list``.
    
    :param smiles_list: A list of smiles strings where each element represents a molecule for which to generate 
        the protonation states.
    :param min_ph: The min pH
    :param max_ph: The max pH
    :param max_variants: The max number of different protonation states to generate for each molecule.
    
    :returns: A list of SMILES which has either the same number of elements as the originally given list 
        but more likely has more elements.
    """
    dmph = DimorphiteDL(
        min_ph=min_ph,
        max_ph=max_ph,
        max_variants=max_variants,
        label_states=False,
        pka_precision=1.0,
    )
    result = list(itertools.chain.from_iterable(
        [dmph.protonate(smiles) for smiles in smiles_list]
    ))

    return result


# Counterfactual generation
# --------------------------
# The following functions implement the neighborhood generation process for the molecular graphs 


def get_neighborhood(smiles: str,
                     atom_valence_map=DEFAULT_ATOM_VALENCE_MAP,
                     mol_filters: t.Sequence[t.Callable[[Chem.Mol], t.List[str]]] = (
                        is_bridge_head_carbon,
                        is_nitrogen_nitrogen_sulfur,
                        is_single_atom,
                     ),
                     use_atom_additions: bool = True,
                     use_bond_additions: bool = False,
                     use_bond_removals: bool = True,
                     fix_protonation: bool = False,
                     max_ph: float = 6.4,
                     min_ph: float = 6.4,
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
    :param use_atom_additions: Whether to generate neighbors through the addition of new atoms.
    :param use_bond_additions: Whether to generate neighbors through the formation of new bonds between existing 
        elements.
    :param use_bond_removals: Whether to generate neighbors through the modification of existing bonds
    :param fix_protonation: Whether to apply the dimorphite_dl tool to generate the accurate protonation states 
        for the molecules within a certain pH scale.
        
    :returns: A list of strings
    """
    neighbors = []

    mol = Chem.MolFromSmiles(smiles)
    free_valence_indices_map = get_free_valence_map(mol)

    if use_atom_additions:
        neighbors += get_valid_atom_additions(
            mol=mol,
            atom_valence_map=atom_valence_map,
            free_valence_indices_map=free_valence_indices_map,
        )

    if use_bond_additions:
        neighbors += get_valid_bond_additions(
            mol=mol,
            free_valence_indices_map=free_valence_indices_map,
        )

    if use_bond_removals:
        neighbors += get_valid_bond_removals(
            mol=mol
        )

    # 15.05.23 - All of the above functions will create "valid" molecular SMILES in the sense that RDKit
    # does not tell us that they are completely wrong, but the molecules that are created might still not
    # realistically exist in chemistry, which is why we additionally apply a set of filters that decide
    # for each molecule if it should be included
    neighbors_filtered = []
    for data in neighbors:
        mol = data['mol']
        
        # 23.10.23 - This additional filter is added to address a bug. It turns out that in some very rare 
        # cases, the procedure actually produces completely invalid molecules that RDKIT cant even parse.
        # and with this addtional line we filter these from the output.
        if not Chem.MolFromSmiles(data['value']):
            continue
        
        # If any of the filters returns "True" then we will omitt that corresponding molecules
        if any([func(mol) for func in mol_filters]):
            continue

        neighbors_filtered.append(data)

    neighbors = neighbors_filtered

    # 01.06.2023 - One problem we have encountered when dealing with the counterfactuals for the aggregation
    # task is that the generated molecules are often times not "realistic" in the sense that they are not
    # correctly protonated for the target pH range of the task.
    # This is why we introduce the option to fix the protonation state of these with an external tool here.
    if fix_protonation:
        neighbors_protonated = []
        dmph = DimorphiteDL(
            min_ph=min_ph,
            max_ph=max_ph,
            max_variants=10,
            label_states=False,
            pka_precision=1.0,
        )
        for data in neighbors:
            smiles_protonated = dmph.protonate(data['value'])
            for smiles in smiles_protonated:
                neighbors_protonated.append({**data, 'value': smiles})

        neighbors = neighbors_protonated

    return neighbors


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
                             ) -> t.List[dict]:
    """
    Given a molecule ``mol``, this function will return a list of SMILES strings that represent valid
    atom additions to the base molecule.

    :param mol: The original molecule for which to calculate the modifications
    :param atom_valence_map: A dict whose keys are string identifiers of atom types (O, N, S...) and the
        values are the integer valences associated with these atom types.
    :param free_valence_indices_map: A dict whose keys are integer values for possible atom valences in the
        molecule and the values are lists containing integer node indices of all the atoms which have that
        corresponding valence.

    :returns: A list of dictionaries which define the valid modifications. Each dictionary has the keys
        - mol: The mol object of the NEW molecule
        - smiles: The SMILES string of the NEW molecule
        - org: The integer atom index of the ORIGINAL molecule where the modification originates
        - mod: The integer atom index of the NEW molecule where the modification occurred
    """
    results = []

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

                smiles = Chem.MolToSmiles(new_mol, canonical=False, rootedAtAtom=0)
                results.append({
                    'type':     'add_node',
                    'mol':      new_mol,
                    'value':    smiles,
                    'org':      (atom, atom),
                    'mod':      (atom, idx),
                })

    return results


def get_valid_bond_additions(mol: Chem.Mol,
                             free_valence_indices_map: t.Dict[int, t.List[int]]
                             ) -> t.List[dict]:
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
    results = []

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

            smiles = Chem.MolToSmiles(new_mol, canonical=False, rootedAtAtom=0)
            results.append({
                'type':     'modify_edge',
                'mol':      new_mol,
                'value':    smiles,
                'org':      (atom1, atom2),
                'mod':      (atom1, atom2),
            })

    return results


def get_valid_bond_removals(mol: Chem.Mol
                            ) -> t.List[dict]:
    """
    Given a molecule ``mol``, this function will return a list of SMILES strings which represent valid
    bond removals. A bond removal is either a downgrade of an existing bond (e.g. double to single) or the
    removal of a single bond which would mean to disconnect at most a single atom from the rest of the
    molecule!

    :param mol: The base molecule for the removals.

    :returns: A list of strings
    """
    results = []

    value_bond_type_map = {
        0: None,
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    bond_type_value_map = invert_dict(value_bond_type_map)

    for valence in [1, 2, 3]:
        for bond in mol.GetBonds():
            atom_removed = None
            atom1, atom2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bond = Chem.Mol(mol).GetBondBetweenAtoms(atom1, atom2)

            bond_type = bond.GetBondType()
            if bond_type not in value_bond_type_map.values():
                continue

            new_mol = Chem.RWMol(mol)
            Chem.Kekulize(new_mol, clearAromaticFlags=True)

            bond_value = bond_type_value_map[bond_type]
            bond_value -= valence

            # In this case, if the bond value has been evaluated to zero then that means that we want to
            # completely remove that bond.
            if bond_value == 0:
                new_mol.RemoveBond(atom1, atom2)

                if new_mol.GetAtomWithIdx(atom1).GetDegree() == 0:
                    new_mol.RemoveAtom(atom1)
                    atom_removed = atom1
                elif new_mol.GetAtomWithIdx(atom2).GetDegree() == 0:
                    new_mol.RemoveAtom(atom2)
                    atom_removed = atom2

            elif bond_value > 0:
                bond.SetBondType(value_bond_type_map[bond_value])
                index = bond.GetIdx()
                new_mol.ReplaceBond(index, bond)

            else:
                continue

            sanitization_result = Chem.SanitizeMol(new_mol, catchErrors=True)
            if sanitization_result:
                continue

            smiles = Chem.MolToSmiles(new_mol, canonical=False, rootedAtAtom=0)
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

            mod1 = atom1 if atom1 != atom_removed else atom2
            mod2 = atom2 if atom2 != atom_removed else atom1
            results.append({
                'type':     'remove_edge' if atom_removed else 'modify_edge',
                'mol':      new_mol,
                'value':    smiles,
                'org':      (atom1, atom2),
                'mod':      (mod1, mod2),
            })

    return results

