import argparse
from functools import lru_cache
import io

import networkx as nx
import pandas as pd
from rdkit import Chem, DataStructs, Geometry
from rdkit.Chem import (AllChem, rdDepictor, rdFMCS, rdqueries,
                        rdRGroupDecomposition)
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Fingerprints import FingerprintMols
rdDepictor.SetPreferCoordGen(True)


@lru_cache(maxsize=4)
def get_palette(name="OKABE"):
    """Build a palette from a set of selected ones."""
    # "Tol" colormap from https://davidmathlogic.com/colorblind
    TOL = [(51, 34, 136), (17, 119, 51), (68, 170, 153), (136, 204, 238),
           (221, 204, 119), (204, 102, 119), (170, 68, 153), (136, 34, 85)]
    # "IBM" colormap from https://davidmathlogic.com/colorblind
    IBM = [(100, 143, 255), (120, 94, 240), (220, 38, 127), (254, 97, 0),
           (255, 176, 0)]
    # Okabe_Ito colormap from https://jfly.uni-koeln.de/color/
    OKABE = [(230, 159, 0), (86, 180, 233), (0, 158, 115), (240, 228, 66),
             (0, 114, 178), (213, 94, 0), (204, 121, 167)]

    result = [(0, 0, 0)]
    for color in locals().get(name):
        result.append(tuple(channel / 255 for channel in color))

    return result


class MoleculeImage:
    """Creates and writes Images of molecules."""
    # Workflow mostly based on
    # https://rdkit.blogspot.com/2020/10/molecule-highlighting-and-r-group.html
    def __init__(self, molecule, core=None, palette="OKABE"):
        self.molecule = Chem.Mol(molecule)  # RDKit modifies the original
        # Create a new molecule to store the drawing changes
        self.draw_mol = self.molecule

        self.palette = get_palette(name=palette)
        try:
            self.name = molecule.GetProp("_Name")
        except KeyError:
            self.name = ""

        self.core = core  # This is the part of the molecule common to others
        self.size = (400, 400)
        self.fill_rings = False
        self.idx_property = "SourceAtomIdx"
        self.bond_width = 2
        self.atom_radius = 0.4

        self.rings = []
        self.old_new = {}

        # First step: flat the molecule and the core
        self._flatten_molecule()

    def _flatten_molecule(self):
        """Flatten the molecule into 2D space."""
        rdDepictor.Compute2DCoords(self.molecule)
        self.molecule.UpdatePropertyCache()

        if self.core and hasattr(self.core, "smartsString"):
            self.core = Chem.MolFromSmarts(self.core.smartsString)
            rdDepictor.Compute2DCoords(self.core)

    def _find_core(self):
        """Find the MCS (core) structure in our molecule and set up groups."""
        #  then tag matching atom indices in each ligand.
        ps = Chem.AdjustQueryParameters.NoAdjustments()
        ps.makeDummiesQueries = True

        qcore = Chem.AdjustQueryProperties(self.core, ps)
        if self.molecule.HasSubstructMatch(qcore):
            for atom in self.molecule.GetAtoms():
                atom.SetIntProp("SourceAtomIdx", atom.GetIdx())

        self.groups = rdRGroupDecomposition.RGroupDecompose(
            [qcore], [self.molecule], asSmiles=False, asRows=True)[0][0]

    def _fix_isotopes(self):
        """Include the atom map numbers in the substructure search in order to
        try to ensure a good alignment of the molecule to symmetric cores."""
        for atom in (_ for _ in self.core.GetAtoms() if _.GetAtomMapNum()):
            atom.ExpandQuery(
                rdqueries.IsotopeEqualsQueryAtom(200 + atom.GetAtomMapNum()))

        for label, group in self.groups.items():
            if label == 'Core':
                continue
            for atom in group.GetAtoms():
                if (not atom.GetAtomicNum() and atom.GetAtomMapNum() and
                    atom.HasProp('dummyLabel') and atom.GetProp('dummyLabel') == label):
                    # attachment point. the atoms connected to this
                    # should be from the molecule
                    for nbr in [_ for _ in atom.GetNeighbors()
                                if _.HasProp(self.idx_property)]:
                        mAt = self.molecule.GetAtomWithIdx(
                            nbr.GetIntProp(self.idx_property))
                        if mAt.GetIsotope():
                            mAt.SetIntProp('_OrigIsotope', mAt.GetIsotope())
                        mAt.SetIsotope(200 + atom.GetAtomMapNum())

    def _remove_hs(self):
        """Remove unmapped hs so that they don't mess up the depiction.

        Updates the old indexes mapped to new ones"""

        rhps = Chem.RemoveHsParameters()
        rhps.removeMapped = False
        self.draw_mol = Chem.RemoveHs(self.molecule, rhps)
        # Reset the original isotope values and account for the fact that
        #  removing the Hs changed atom indices
        for i, atom in enumerate(self.draw_mol.GetAtoms()):
            if atom.HasProp(self.idx_property):
                self.old_new[atom.GetIntProp(self.idx_property)] = i
                if atom.HasProp("_OrigIsotope"):
                    atom.SetIsotope(atom.GetIntProp("_OrigIsotope"))
                    atom.ClearProp("_OrigIsotope")
                else:
                    atom.SetIsotope(0)

        rdDepictor.GenerateDepictionMatching2DStructure(self.draw_mol, self.core)

    def _add_ring(self, group, color):
        """Add rings found in group to self.rings."""
        Chem.GetSSSR(group)  # This finds and sets the rings.
        ring_info = group.GetRingInfo()
        for aring in ring_info.AtomRings():
            tring = []
            allFound = True
            for aid in aring:
                atom = group.GetAtomWithIdx(aid)
                if not atom.HasProp(self.idx_property):
                    allFound = False
                    break
                tring.append(self.old_new[atom.GetIntProp(self.idx_property)])
            if allFound:
                self.rings.append((tring, color))

    def _highlight_residue_atoms(self, group, color):
        """Highlight the atoms of a group with a given color.

        Return the dict with the atoms to highlight."""
        highlights = {}
        atom_radius = {}
        for atom in group.GetAtoms():
            if atom.HasProp(self.idx_property):
                origIdx = self.old_new[atom.GetIntProp(self.idx_property)]
                highlights[origIdx] = [color]
                atom_radius[origIdx] = self.atom_radius

        return (highlights, atom_radius)

    def _highlight_residue_bonds(self, group, color):
        """Highlight the bonds of a group with a given color.

        Return the dict with the atoms to highlight."""
        highlights = {}
        width_multiplier = {}
        for bond in group.GetBonds():
            batom = bond.GetBeginAtom()
            eatom = bond.GetEndAtom()
            if batom.HasProp(self.idx_property) and eatom.HasProp(self.idx_property):
                origBnd = self.draw_mol.GetBondBetweenAtoms(
                    self.old_new[batom.GetIntProp(self.idx_property)],
                    self.old_new[eatom.GetIntProp(self.idx_property)])
                bndIdx = origBnd.GetIdx()
                highlights[bndIdx] = [color]
                width_multiplier[bndIdx] = self.bond_width

        return (highlights, width_multiplier)

    def _highlight_residue(self, group, color, highlight):
        atoms, atom_radius = self._highlight_residue_atoms(group, color)
        highlight[0].update(atoms)
        highlight[2].update(atom_radius)

        bonds, bond_widths = self._highlight_residue_bonds(group, color)
        highlight[1].update(bonds)
        highlight[3].update(bond_widths)

    def _colorfill_rings(self, highlights):
        """Fill the rings in self.rings with the color requested."""
        if self.fill_rings and self.rings:
            # Set the molecule scale
            self.canvas.DrawMoleculeWithHighlights(self.draw_mol, "", *highlights)
            self.canvas.ClearDrawing()
            canvas_options = self.canvas.drawOptions()

            conf = self.draw_mol.GetConformer()
            for (ring, color) in self.rings:
                ps = [Geometry.Point2D(conf.GetAtomPosition(_)) for _ in ring]
                self.canvas.SetFillPolys(True)
                self.canvas.SetColour(color)
                self.canvas.DrawPolygon(ps)
            canvas_options.clearBackground = False

    def png(self):
        """Create a PNG with the data."""
        # Store which atoms, bonds, and rings we'll be highlighted
        highlights = [{}, {}, {}, {}]

        # Initialize the drawing
        self.canvas = rdMolDraw2D.MolDraw2DCairo(*self.size)
        canvas_options = self.canvas.drawOptions()
        canvas_options.useBWAtomPalette()

        if self.core:
            self._find_core()
            self._fix_isotopes()
            self._remove_hs()

            # Loop over R groups.
            for color, (label, group) in zip(self.palette, self.groups.items()):
                if label == "Core":
                    continue

                self._highlight_residue(group, color, highlights)

                if self.fill_rings:
                    self._add_ring(group, color)

            if self.fill_rings and self.rings:
                self._colorfill_rings(highlights)

        # Draw the molecule, with highlights
        self.canvas.DrawMoleculeWithHighlights(self.draw_mol, "", *highlights)
        self.canvas.FinishDrawing()

        return self.canvas.GetDrawingText()  # This is a Png as a b"" string


class MapGen():
    def __init__(self, in_sdf, network_obj=None):
        if hasattr(in_sdf, "seek"):
            in_sdf.seek(0)
            self.suppl = Chem.ForwardSDMolSupplier(in_sdf)
        else:
            self.suppl = Chem.SDMolSupplier(str(in_sdf))

        self.network = network_obj
        self.metric = network_obj.metric
        self.ligands = {}
        self.simF = None
        self.molecules = None
        self.mcs = None

        self._set_ligands()
        self._set_similarity_function()
        self._set_mcs_core()
        self._set_similarity_matrix()

    def fingerprint(self, molecule):
        if self.metric == self.network.Tanimoto:
            return FingerprintMols.FingerprintMol(molecule)
        elif self.metric == self.network.MFP:
            return AllChem.GetMorganFingerprintAsBitVect(
                molecule, 2, nBits=2048)
        elif self.metric == self.network.SMILES:
            return Chem.MolToSmiles(molecule, isomericSmiles=True)

    def _set_mcs_core(self):
        """Set the core MCS (Maximum Common Substructure) for the network."""
        # Flatten all molecules in all charge groups
        self.molecules = [mol for charge in self.ligands.values()
                          for mol in charge["Mol"]]

        self.mcs = rdFMCS.FindMCS(self.molecules,
                                  matchValences=False,
                                  ringMatchesRingOnly=True,
                                  completeRingsOnly=True,
                                  matchChiralTag=False)

    def _set_similarity_function(self):
        """Set the similarity function to be used with selected metric."""
        if self.metric in [self.network.Tanimoto, self.network.MFP]:
            self.simF = DataStructs.FingerprintSimilarity
        elif self.metric == self.network.MCS:
            self.simF = Chem.rdFMCS.FindMCS
        elif self.metric == self.network.SMILES:
            from Bio import pairwise2
            self.simF = pairwise2.align.globalms

    def _ligands_score(self, data, lig_i, lig_j):
        """Return a similarity score between lig_i and lig_j using self.simF."""
        if lig_i == lig_j:
            return 100.0 if self.metric in [self.network.SMILES, self.network.MCS] else 1.0
        if self.metric in [self.network.Tanimoto, self.network.MFP]:
            return self.simF(data['FP'][lig_i], data['FP'][lig_j])
        if self.metric == self.network.MCS:
            score = self.simF([data['Mol'][lig_i],
                               data['Mol'][lig_j]],
                              atomCompare=rdFMCS.AtomCompare.CompareAny)
            return score.numAtoms + score.numBonds
        if self.metric == self.network.SMILES:
            return self.simF(
                data['FP'][lig_i], data['FP'][lig_j],
                1, -1, -0.5, -0.05)[0].score

    def _set_ligands(self):
        for mol in self.suppl:
            charge = Chem.rdmolops.GetFormalCharge(mol)
            v = self.ligands.setdefault(charge, {'Name': [], 'Mol': [], 'FP': []})
            v['Name'].append(mol.GetProp('_Name'))
            v['Mol'].append(mol)
            if self.metric != self.network.MCS:
                v['FP'].append(self.fingerprint(mol))

    def _set_similarity_matrix(self):
        # Ensure the self.ligands has been created
        if not self.ligands:
            self.set_ligands()
        for charge, ligand in self.ligands.items():
            df = pd.DataFrame()
            ligands_done = []
            for i, lig1 in enumerate(ligand['Name']):
                for j, lig2 in enumerate(ligand['Name']):
                    if lig2 not in ligands_done:
                        df.loc[lig2, lig1] = self._ligands_score(ligand, i, j)
                ligands_done.append(lig1)
            self.ligands[charge]['df'] = df

    def set_ligpairs(self):
        for charge, value in self.ligands.items():
            pairs_dict = {}
            for i in value['df'].index:
                for j in value['df'].columns:
                    if value['df'].loc[i, j] not in [1.0, 100] and not pd.isnull(value['df'].loc[i, j]):
                        pairs_dict['{} {}'.format(i, j)] = round(value['df'].loc[i, j], 3)

            if self.metric in [self.network.Tanimoto, self.network.MFP, self.network.MCS]:
                pairs_dict = {k: v for k, v in sorted(pairs_dict.items(), key=lambda item: item[1], reverse=True)}
            elif self.metric == self.network.SMILES:
                pairs_dict = {k: v for k, v in sorted(pairs_dict.items(), key=lambda item: item[1], reverse=False)}

            value['pairs_dict'] = pairs_dict

    def process_map(self):
        self._set_similarity_matrix()
        self.set_ligpairs()
        self.make_map()

    def intersection(self, edge_list, candidate_edge):
        r1, r2 = candidate_edge.split()[0], candidate_edge.split()[1]
        for edge in edge_list:
            if r1 == edge[0] or r1 == edge[1] or r2 == edge[0] or r2 == edge[1]:
                return True # Shortcut comparing: it's already True
        return False

    def not_ingraph(self, node_list, candidate_edge):
        r1, r2 = candidate_edge.split()[0], candidate_edge.split()[1]
        return r1 not in node_list or r2 not in node_list

    def outer_nodes(self, G):
        node_list = []
        for node in G.nodes:
            if len(G.edges(node)) == 1:
                node_list.append(node)
        return node_list

    def make_map(self):
        for charge, lig in self.ligands.items():
            if "pairs_dict" not in lig:
                # Call set_ligpairs the first iteration if it wasn't called before
                self.set_ligpairs()
            H = nx.Graph()
            lpd = self.ligands[charge]['pairs_dict']
            simdf = self.ligands[charge]['df']
            if len(lig['Name']) == 1:  # In case one ligand is found alone in a charge group
                ligcol = simdf[lig['Name']].sort_values(by=[lig['Name'][0]]) #complete similarity matrix
                H.add_edge(lig['Name'][0], ligcol.index[1])
                H.add_edge(lig['Name'][0], ligcol.index[2])
                lig['Graph'] = H
                break

            # 1. Make SPT
            incomplete = True
            while incomplete:
                for pert, score in lpd.items():
                    if len(H.nodes) == len(lig['Name']):
                        incomplete = False
                        break
                    l1, l2 = pert.split()[0], pert.split()[1]
                    if H.has_edge(l1, l2) or H.has_edge(l2, l1):
                        continue
                    if len(H.nodes) == 0 or self.intersection(H.edges, pert) and self.not_ingraph(H.nodes, pert):
                        H.add_edge(l1, l2, weight=score)
                        break

            # 2. Close Cycles
            while len(self.outer_nodes(H)) != 0:
                for pert, score in lpd.items():
                    l1, l2 = pert.split()[0], pert.split()[1]
                    if l1 in self.outer_nodes(H) or l2 in self.outer_nodes(H):
                        if (l1, l2) not in H.edges or (l2, l1) not in H.edges:
                            H.add_edge(l1, l2, weight=score)
                            break

            # 3. Add influence edges
            eig_cent = nx.eigenvector_centrality(H, max_iter=1000)
            eig_cent = {k: v for k, v in sorted(eig_cent.items(), key=lambda item: item[1], reverse=True)}
            per_nodes = [k for k,v in nx.eigenvector_centrality(H, max_iter=1000).items() if v < 0.01]
            per_len = len(per_nodes)
            while per_len > 1:
                for pert, score in lpd.items():
                    l1, l2 = pert.split()[0], pert.split()[1]
                    if l1 in per_nodes and l2 not in per_nodes or l1 not in per_nodes and l2 in per_nodes:
                        if (l1, l2) not in H.edges or (l2, l1) not in H.edges and intersection(H.edges, pert):
                            H.add_edge(l1, l2, weight=score)
                            nlen = len([v for k,v in nx.eigenvector_centrality(H, max_iter=1000).items() if v < 0.01])
                            if nlen > per_len:
                                H.remove_edge(l1, l2)
                                continue
                            else:
                                per_len = nlen
                                break
            lig['Graph'] = H


# The following code is the CLI
def getParser():
    parser = argparse.ArgumentParser(
        prog='MapGen',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='FEP map generator based on selected distance metrics.')
    parser.add_argument('-i', '--insdf',
                        dest="isdf",
                        required=True,
                        help=".sdf file name")
    parser.add_argument('-m', '--metric',
                        dest="metric",
                        default=g.MFP,
                        choices=[g.MFP, g.Tanimoto, g.MCS, g.SMILES],
                        required=False,
                        help="Distance metric for ligand pairwairse comparison")
    return parser

def main():
    args = getParser().parse_args()
    ## Put file in memory stream. This allows the server to read uploaded file
    ##  into memory and pass it as an io.BytesIO() to MapGen
    with open(args.isdf, "rb") as f:
        with io.BytesIO(f.read()) as fio:
            mg = MapGen(fio, args.metric)
            mg.process_map()

if __name__ == "__main__":
    main()
