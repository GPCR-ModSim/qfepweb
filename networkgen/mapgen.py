import argparse
from collections import namedtuple
from functools import cached_property, lru_cache
import io
import itertools
import operator

import networkx as nx
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

Ligand = namedtuple("Ligand", ["name", "pool_idx", "fingerprint"])

class MoleculePool:
    """Keeps track of a group of molecules, adding the properties needed.

    mcs: The Maximum Common Substructure (the core).
    core: The common core for all molecules as a Mol object.
    query_core: Same as core, but to use with decomposition.
    groups: A list of Core + Residues found for each molecule."""

    def __init__(self, molecules=None):
        self.molecules = []
        if molecules:
            self.molecules = molecules

    def __getitem__(self, idx):
        return self.molecules[idx]

    def __len__(self):
        return len(self.molecules)

    def append(self, item):
        """Adds a new molecule to the pool."""
        self.molecules.append(item)
        # Clear the cached properties so they get re-calculated
        for attr in ["mcs", "core", "query_core", "groups"]:
            try:
                delattr(self, attr)
            except AttributeError:
                pass

    @cached_property
    def mcs(self):
        if len(self) == 0:
            return
        if len(self) == 1:
            # If the "pool" is only one, create a pool of two to trick FindMCS.
            #  The MCS returned will be the whole molecule.
            # This is neeced because FindMCS returns a non-instantiable
            #  ResultMCS object.
            molecules = [self.molecules[0], self.molecules[0]]
        else:
            molecules = self.molecules

        return rdFMCS.FindMCS(molecules,
                              matchValences=False,
                              ringMatchesRingOnly=True,
                              completeRingsOnly=True,
                              matchChiralTag=False)

    @cached_property
    def core(self):
        if self.mcs:
            return Chem.MolFromSmarts(self.mcs.smartsString)

    @cached_property
    def query_core(self):
        if self.core:
            ps = Chem.AdjustQueryParameters.NoAdjustments()
            ps.makeDummiesQueries = True
            return Chem.AdjustQueryProperties(self.core, ps)

    @cached_property
    def groups(self):
        molecule_matches = []
        if self.query_core:
            for molecule in self.molecules:
                if molecule.HasSubstructMatch(self.query_core):
                    molecule_matches.append(molecule)
                    for atom in molecule.GetAtoms():
                        atom.SetIntProp("SourceAtomIdx", atom.GetIdx())

            return Chem.rdRGroupDecomposition.RGroupDecompose(
                [self.query_core], molecule_matches, asSmiles=False,
                asRows=True)[0]
        else:
            return []


class MoleculeImage:
    """Creates and writes Images of molecules."""
    # Workflow mostly based on
    # https://rdkit.blogspot.com/2020/10/molecule-highlighting-and-r-group.html
    def __init__(self, pool_idx=0, pool=None, palette="OKABE"):
        """Loads a molecule to create its image.

        As each molecule needs a core to refer to (orienting, coloring...),
        load a `pool` of them, pointing which one of the pool are you interested:

            pool = MoleculePool([Mol1, Mol2, ..., Moln])

            # Load Image for Mol1 above
            MoleculeImage(pool=pool, pool_idx=0)
            # For Mol2
            MoleculeImage(pool=pool, pool_idx=1)

        """
        self.pool = pool  # This is the info of all other molecules
        self.pool_idx = pool_idx
        self.pool.groups  # This initializes the pool attributes

        # Create a new molecule to store the drawing changes
        # RDKit modifies the original
        self.molecule = Chem.Mol(self.pool[self.pool_idx])
        self.draw_mol = self.molecule

        self.palette = get_palette(name=palette)
        try:
            self.name = self.molecule.GetProp("_Name")
        except KeyError:
            self.name = ""

        self.core = Chem.Mol(self.pool.core)  # Copy the core
        rdDepictor.Compute2DCoords(self.core) # This reorients the core
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

    def _fix_isotopes(self):
        """Include the atom map numbers in the substructure search in order to
        try to ensure a good alignment of the molecule to symmetric cores."""
        for atom in (_ for _ in self.core.GetAtoms() if _.GetAtomMapNum()):
            atom.ExpandQuery(
                rdqueries.IsotopeEqualsQueryAtom(200 + atom.GetAtomMapNum()))

        for label, group in self.pool.groups[self.pool_idx].items():
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

    def _reorient_molecule(self):
        """Reorients the drawing molecule following the core."""
        rdDepictor.GenerateDepictionMatching2DStructure(
            self.draw_mol, self.core)

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
        """Create a PNG with the data.

            imgr = mapgen.MoleculeImage(pool_idx=1, pool=pool)

            with open("MyMolecule.png", "wb") as r:
                h.write(imgr.png())

        """
        # Store which atoms, bonds, and rings will be highlighted
        highlights = [{}, {}, {}, {}]

        # Initialize the drawing
        self.canvas = rdMolDraw2D.MolDraw2DCairo(*self.size)
        canvas_options = self.canvas.drawOptions()
        canvas_options.useBWAtomPalette()

        if self.core:
            self._fix_isotopes()      #
            self._remove_hs()         # Does this three belongs here ??
            self._reorient_molecule() #
            # Loop over R groups.
            for color, (label, group) in zip(
                    self.palette, self.pool.groups[self.pool_idx].items()):
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


class MapGen:
    def __init__(self, network_obj=None, in_sdf=""):
        """Creates a Network from a Network Generator object.

        A Network Generator is a model, in which we need one field: metric

        If this model doesn't have a File-like at in_sdf, an alternative
        in_sdf parameter can be provided to work with."""
        if network_obj and hasattr(network_obj, "in_sdf") and \
                hasattr(network_obj.in_sdf, "seek"):
            network_obj.in_sdf.seek(0)
            self.suppl = Chem.ForwardSDMolSupplier(network_obj.in_sdf)
        else:
            self.suppl = Chem.SDMolSupplier(str(in_sdf))

        self.network = network_obj
        self.metric = network_obj.metric
        self.pool = MoleculePool()
        self.ligands = {}
        self.simF = None

        self._set_ligands()
        self._set_similarity_function()
        self._set_similarity_matrix()

    def fingerprint(self, molecule):
        if self.metric == self.network.Tanimoto:
            return FingerprintMols.FingerprintMol(molecule)
        elif self.metric == self.network.MFP:
            return AllChem.GetMorganFingerprintAsBitVect(
                molecule, 2, nBits=2048)
        elif self.metric == self.network.SMILES:
            return Chem.MolToSmiles(molecule, isomericSmiles=True)

    def _set_similarity_function(self):
        """Set the similarity function to be used with selected metric."""
        if self.metric in [self.network.Tanimoto, self.network.MFP]:
            self.simF = DataStructs.FingerprintSimilarity
        elif self.metric == self.network.MCS:
            self.simF = Chem.rdFMCS.FindMCS
        elif self.metric == self.network.SMILES:
            from Bio import pairwise2
            self.simF = pairwise2.align.globalms

    def _ligands_score(self, l1_idx, l2_idx):
        """Return a similarity score between l1 and l2 using self.simF."""
        if l1_idx == l2_idx:
            return 100.0 if self.metric in [self.network.SMILES, self.network.MCS] else 1.0

        if self.metric == self.network.MCS:
            molecule_1 = self.pool[l1_idx]
            molecule_2 = self.pool[l2_idx]
            score = self.simF([molecule_1, molecule_2],
                              atomCompare=rdFMCS.AtomCompare.CompareAny)
            return score.numAtoms + score.numBonds
        else:
            fingerprint_1 = self.fingerprint(self.pool[l1_idx])
            fingerprint_2 = self.fingerprint(self.pool[l2_idx])
            if self.metric in [self.network.Tanimoto, self.network.MFP]:
                return round(self.simF(fingerprint_1, fingerprint_2), 3)
            if self.metric == self.network.SMILES:
                return round(self.simF(fingerprint_1, fingerprint_2,
                                       1, -1, -0.5, -0.05)[0].score, 3)

    def _set_ligands(self):
        for idx, mol in enumerate(self.suppl):
            self.pool.append(mol)
            charge = Chem.rdmolops.GetFormalCharge(mol)
            v = self.ligands.setdefault(charge, {"Ligand": []})
            ligand = Ligand(name=mol.GetProp('_Name'),
                            pool_idx=idx,
                            fingerprint=self.fingerprint(mol))

            v["Ligand"].append(ligand)

    def _set_similarity_matrix(self):
        # TODO document this and change function name (not a matrix anymore).
        # Ensure the self.ligands has been created
        if not self.ligands:
            self.set_ligands()

        reverse = self.metric == self.network.SMILES
        for charge, ligands in self.ligands.items():
            ligands["Scores"] = {}
            for l1, l2 in itertools.combinations(ligands["Ligand"], 2):
                ligands["Scores"][(l1.pool_idx, l2.pool_idx)] = \
                    self._ligands_score(l1.pool_idx, l2.pool_idx)

            ligands["Scores"] = \
                {k: v for k, v in sorted(ligands["Scores"].items(),
                                         key=operator.itemgetter(1),
                                         reverse=reverse)}

    def build_charge_tree(self, ligands):
        """Build a graph with all ligands provided."""
        #if "pairs_dict" not in ligands:
        #    # Call set_ligpairs the first iteration if it wasn't called before
        #    self.set_ligpairs()
        H = nx.Graph()

        if len(ligands['Ligand']) == 1:
            # In case one ligand is found alone in a charge group
            # A "graph" of one node and no edges is created.
            H.add_node(ligands["Ligand"][0].pool_idx)
        elif len(ligands['Ligand']) == 2:
            # In case two ligands are found in a charge group
            # Complete similarity matrix. At this point, stop graph
            # generation because two nodes in a graph will always result
            # in the same graph.
            H.add_edge(ligands['Ligand'][0].pool_idx,
                       ligands['Ligand'][1].pool_idx,
                       weight=ligands["Scores"][(0, 1)])
        else:
            incomplete = True
            while incomplete:
                for (l1, l2), score in ligands["Scores"].items():
                    if len(H.nodes) == len(ligands['Ligand']):
                        # All nodes has been added to the graph
                        incomplete = False
                        break
                    if H.has_edge(l1, l2) or H.has_edge(l2, l1):
                        # Both Nodes already in graph and connected
                        continue
                    if len(H.nodes) == 0 or self.intersection(H.edges, l1, l2) \
                            and self.not_ingraph(H.nodes, l1, l2):
                        H.add_edge(l1, l2, weight=score)
                        break
        return H

    def close_cycles(self, ligands):
        """Close open cycles in graph.

        An outer node (open cycle) is defined as a Node with only one Edge.
        """
        outer_nodes = self.outer_nodes(ligands["Graph"])
        while True:
            added_edge = False
            for (l1, l2), score in ligands["Scores"].items():
                if l1 in outer_nodes or l2 in outer_nodes:
                    if (l1, l2) not in ligands["Graph"].edges or \
                            (l2, l1) not in ligands["Graph"].edges:
                        ligands["Graph"].add_edge(l1, l2, weight=score)
                        # signal that an edge has been added in this iteration
                        # of the while loop.
                        added_edge = True
                        break

            outer_nodes = self.outer_nodes(ligands["Graph"])
            # if no edge has been added, the while loop will hang. Exit.
            if not added_edge or len(outer_nodes) == 0:
                break

    def add_influence_edges(self, ligands):
        def under_cent(graph):
            return [v for k, v in
                    nx.eigenvector_centrality(graph, max_iter=1000).items()
                    if v < 0.01]

        eig_cent = {k: v for k, v in sorted(
            nx.eigenvector_centrality(ligands["Graph"], max_iter=1000).items(),
            key=lambda item: item[1], reverse=True)}
        per_nodes = [k for k, v in eig_cent.items() if v < 0.01]
        per_len = len(per_nodes)

        # FIXME: Potential hang loop here
        while per_len > 1:
            for (l1, l2), score in ligands["Scores"].items():
                if (l1 in per_nodes and l2 not in per_nodes) or \
                        (l1 not in per_nodes and l2 in per_nodes):
                    if (l1, l2) not in ligands["Graph"].edges or \
                            (l2, l1) not in ligands["Graph"].edges and \
                            intersection(ligands["Graph"].edges, l1, l2):
                        ligands["Graph"].add_edge(l1, l2, weight=score)
                        nlen = len(under_cent(ligands["Graph"]))
                        if nlen > per_len:
                            ligands["Graph"].remove_edge(l1, l2)
                            continue
                        else:
                            per_len = nlen
                            break

    def process_map(self):
        self._set_similarity_matrix()
        self.make_map()

    def intersection(self, edge_list, r1, r2):
        for edge in edge_list:
            if r1 == edge[0] or r1 == edge[1] or r2 == edge[0] or r2 == edge[1]:
                return True # Shortcut comparing: it's already True
        return False

    def not_ingraph(self, node_list, r1, r2):
        return r1 not in node_list or r2 not in node_list

    def outer_nodes(self, G):
        node_list = []
        for node in G.nodes:
            if len(G.edges(node)) == 1:
                node_list.append(node)
        return node_list

    def make_map(self):
        for charge, ligands in self.ligands.items():
            # 1. Make SPT (Shortest Path Tree)
            ligands["Graph"] = self.build_charge_tree(ligands)
            # 2. Close Cycles
            self.close_cycles(ligands)
            # 3. Add influence edges
            self.add_influence_edges(ligands)

# The following code is the CLI
def getParser():
    # FIXME This is broken cause models.py cannot be imported (circular)
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
    # FIXME This is broken cause models.py cannot be imported (circular)
    args = getParser().parse_args()
    ## Put file in memory stream. This allows the server to read uploaded file
    ##  into memory and pass it as an io.BytesIO() to MapGen
    with open(args.isdf, "rb") as f:
        with io.BytesIO(f.read()) as fio:
            mg = MapGen(args.metric)
            mg.process_map()

if __name__ == "__main__":
    main()
