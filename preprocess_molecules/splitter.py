from rdkit.Chem import BRICS
from rdkit import Chem
import numpy as np

class Splitter:
    def __init__(self):
        self.i = 0
        self.matrix = []
        self.broken_bonds = []
        self.count = 0

    def _get_broken_bonds(self, mol: Chem.rdchem.Mol):
        brics_bonds = list(BRICS.FindBRICSBonds(mol))

        broken_bonds = [(bond[0][0], bond[0][1]) for bond in brics_bonds]
        self.broken_bonds = broken_bonds

    def _mol_to_adjacency_matrix(self, mol: Chem.rdchem.Mol):
        num_atoms = mol.GetNumAtoms()

        adj_matrix = np.zeros((num_atoms, num_atoms), dtype=int)
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            bond_type = int(bond.GetBondTypeAsDouble())

            adj_matrix[i][j] = bond_type
            adj_matrix[j][i] = bond_type
        return adj_matrix

    def _get_subgraph_adjacency_matrix(self, full_adj_matrix, subgraph_vertices):
        full_adj = np.array(full_adj_matrix)

        sub_adj = full_adj[np.ix_(subgraph_vertices, subgraph_vertices)]

        return sub_adj

    def _get_masked_adjacency_matrix(self, mol: Chem.rdchem.Mol):
        matrix = self._mol_to_adjacency_matrix(mol)
        for a,b in self.broken_bonds:
            matrix[a][b] = -1
            matrix[b][a] = -1
        self.matrix = matrix
        return matrix

    def _find_connected_components(self, adj_matrix: np.array):
        n = len(adj_matrix)
        visited = [False] * n

        def dfs(v, component):
            visited[v] = True
            component.append(v)
            for u in range(n):
                if adj_matrix[v][u] != 0 and adj_matrix[v][u] != -1 and not visited[u]:
                    dfs(u, component)

        components = []
        for v in range(n):
            if not visited[v]:
                component = []
                dfs(v, component)
                components.append(component)

        return components

    def _get_merged_substructures(self, list_of_atoms_substructures):
        merged_substr = []
        cnt = 0
        for pair in self.broken_bonds:
            a,b = pair
            res_substr = []
            for substr in list_of_atoms_substructures:
                if a in substr or b in substr:
                    res_substr += substr
                    cnt += 1
                if cnt == 2:
                    cnt = 0
                    break
            cnt = 0
            merged_substr.append(res_substr)
        return merged_substr

    def _get_substructure_smiles_fragments_not_merged(self, mol, components): 
        emol = Chem.RWMol(mol)

        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            in_same = any(a1 in comp and a2 in comp for comp in components)
            if not in_same:
                emol.RemoveBond(a1, a2)

                for i in range(int(mol.GetBondBetweenAtoms(a1, a2).GetBondTypeAsDouble())):
                    new_atom1 = emol.AddAtom(Chem.Atom(1))
                    emol.AddBond(a1, new_atom1, Chem.rdchem.BondType.SINGLE)

                    new_atom2 = emol.AddAtom(Chem.Atom(1))
                    emol.AddBond(a2, new_atom2, Chem.rdchem.BondType.SINGLE)

        frags = Chem.GetMolFrags(emol.GetMol(), asMols=True)
        return [Chem.MolToSmiles(Chem.RemoveHs(frag)) for frag in frags]

    def _removing_bonds(self, current_component, adj_matr):
        removing_bonds = []
        for i in range(len(adj_matr)):
            for j in range(i, len(adj_matr)):

                if adj_matr[i][j] != 0:
                    if i in current_component and j not in current_component:
                        removing_bonds.append((i,j))

                    elif i not in current_component and j in current_component:
                        removing_bonds.append((i,j))
        return removing_bonds

    def _get_true_substr(self, adj_matr, frags, component):
        true_adj_matr = self._get_subgraph_adjacency_matrix(adj_matr, component)
        for frag in frags:
            if np.array_equal(self._mol_to_adjacency_matrix(Chem.RemoveHs(frag)), true_adj_matr):
                return frag

    def _get_substructure_smiles_fragments_merged(self, mol, components): 
        adj_matr = self._mol_to_adjacency_matrix(mol)
        bonds = []
        smiles_of_frags = []
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            bonds.append((a1,a2))

        for component in components:
            removing_bnds = []
            removing_bnds = self._removing_bonds(component, adj_matr)
            emol = Chem.RWMol(mol)
            for a1, a2 in removing_bnds:
                emol.RemoveBond(a1, a2)

                for i in range(int(mol.GetBondBetweenAtoms(a1, a2).GetBondTypeAsDouble())):
                    new_atom1 = emol.AddAtom(Chem.Atom(1))
                    emol.AddBond(a1, new_atom1, Chem.rdchem.BondType.SINGLE)

                    new_atom2 = emol.AddAtom(Chem.Atom(1))
                    emol.AddBond(a2, new_atom2, Chem.rdchem.BondType.SINGLE)

            frags = Chem.GetMolFrags(emol.GetMol(), asMols=True)
            true_frag = self._get_true_substr(adj_matr, frags, component)
            smi = Chem.MolToSmiles(Chem.RemoveHs(true_frag))
            smiles_of_frags.append(smi)
        return smiles_of_frags

    def _get_substructures_merged(self, smiles: str):
        mol = Chem.MolFromSmiles(smiles)
        mol.UpdatePropertyCache(strict=True)

        broken_bonds = self._get_broken_bonds(mol)
        adj_matrix = self._get_masked_adjacency_matrix(mol, broken_bonds)

        list_of_connected_components = self._find_connected_components(adj_matrix)
        merged_substructures = self._get_merged_substructures(broken_bonds, list_of_connected_components)

        return merged_substructures

    def get_substructures_smiles_and_merged(self, smiles: str):
        mol = Chem.MolFromSmiles(smiles)
        mol.UpdatePropertyCache(strict=True)

        self._get_broken_bonds(mol)
        adj_matrix = self._get_masked_adjacency_matrix(mol)


        list_of_connected_components = self._find_connected_components(adj_matrix)
        list_of_connected_components = [sorted(i) for i in list_of_connected_components]

        list_of_smiles = self._get_substructure_smiles_fragments_not_merged(mol, list_of_connected_components)
        merged_substructures = self._get_merged_substructures(list_of_connected_components)

        merged_substructures = [sorted(i) for i in merged_substructures]

        list_of_smiles.append('$')
        list_of_smiles += self._get_substructure_smiles_fragments_merged(mol, merged_substructures)
        return list_of_smiles