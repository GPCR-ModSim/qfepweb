import argparse
import io
import json
import os
import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Fingerprints import FingerprintMols


class MapGen():
    def __init__(self, in_sdf, metric):
        self.suppl = Chem.ForwardSDMolSupplier(in_sdf)
        self.metric = metric.lower()
        self.lig_dict = {}
        self.simF = None  # The similarity function to calculate distance
                          #  between ligands
        self._set_similarity_function()

    def make_fp(self, mol):
        if self.metric == 'tanimoto':
            fp = FingerprintMols.FingerprintMol(mol)
        if self.metric == 'mfp':
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        if self.metric == 'smiles':
            fp = Chem.MolToSmiles(mol, isomericSmiles=True)
        return fp

    def _set_similarity_function(self):
        """Set the similarity function to be used with selected metric."""
        if self.metric in ['tanimoto', 'mfp']:
            from rdkit.DataStructs import FingerprintSimilarity as f
        elif self.metric == 'mcs':
            from rdkit.Chem.rdFMCS import FindMCS as f
        elif self.metric == 'smiles':
            from Bio.pairwise2 import align
            f = align.globalms

        self.simF = f

    def _ligands_score(self, data, lig_i, lig_j):
        """Return a similarity score between lig_i and lig_j using self.simF."""
        if lig_i == lig_j:
            return 100.0 if self.metric in ["smiles", "mcs"] else 1.0
        if self.metric in ['tanimoto', 'mfp']:
            return self.simF(data['FP'][lig_i], data['FP'][lig_j])
        if self.metric == 'mcs':
            score = self.simF([data['Mol'][lig_i],
                               data['Mol'][lig_j]],
                              atomCompare=rdFMCS.AtomCompare.CompareAny)
            return score.numAtoms + score.numBonds
        if self.metric == 'smiles':
            return self.simF(
                data['FP'][lig_i], data['FP'][lig_j], 1, -1, -0.5, -0.05)[0][2]

    def set_ligdict(self):
        self.lig_dict = {}
        for mol in self.suppl:
            charge = Chem.rdmolops.GetFormalCharge(mol)
            v = self.lig_dict.setdefault(charge, {'Name': [], 'Mol': [], 'FP': []})
            v['Name'].append(mol.GetProp('_Name'))
            v['Mol'].append(mol)
            if self.metric != 'mcs':
                v['FP'].append(self.make_fp(mol))

    def sim_mx(self):
        # Ensure the self.lig_dict has been created
        if not self.lig_dict:
            self.set_ligdict()

        for charge, ligand in self.lig_dict.items():
            df = pd.DataFrame()
            ligands_done = []
            for i, lig1 in enumerate(ligand['Name']):
                for j, lig2 in enumerate(ligand['Name']):
                    if lig2 not in ligands_done:
                        df.loc[lig2, lig1] = self._ligands_score(ligand, i, j)
                ligands_done.append(lig1)
            self.lig_dict[charge]['df'] = df

    def clean_mxs(self):
        for charge, value in self.lig_dict.items():
            for i, j in zip(value['df'].index, value['df'].idxmax()):
                vlist = value['df'].loc[i, :].tolist()
                index = vlist.index(1.0)
                vlist[index:] = [0.0] * len(vlist[index:])
                value['df'].loc[i, :] = vlist

            if self.metric in ['tanimoto', 'mfp']:
                df = 1 - value['df']  # get dissimilarity matrix

            if self.metric == 'mcs':
                df = 100 - value['df']  # get dissimilarity matrix
            value['df'] = df.replace(1.0, 0.0)  # set diagonal to 0
            value['df'] = df.replace(0.0, 1.0)  # set zeroes into 1 (in order to search shortest path)

    def set_ligpairs(self):
        for charge, value in self.lig_dict.items():
            pairs_dict = {}
            for i in value['df'].index:
                for j in value['df'].columns:
                    if value['df'].loc[i, j] != 1.0:
                        pairs_dict['{} {}'.format(i, j)] = round(value['df'].loc[i, j], 3)

            if self.metric in ['tanimoto', 'mfp', 'mcs']:
                pairs_dict = {k: v for k, v in sorted(pairs_dict.items(), key=lambda item: item[1], reverse=False)}
            elif self.metric == 'smiles':
                pairs_dict = {k: v for k, v in sorted(pairs_dict.items(), key=lambda item: item[1], reverse=True)}

            value['pairs_dict'] = pairs_dict

    def process_map(self):
        self.set_ligdict()
        self.sim_mx()
        self.clean_mxs()
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
        for charge, lig in self.lig_dict.items():
            H = nx.Graph()
            if len(lig['Name']) == 1:  # In case one ligand is found alone in a charge group
                ligcol = lig['Name'].sort_values(by=[lig['Name'][0]]) #complete similarity matrix
                H.add_edge(lig['Name'][0], ligcol.index[1])
                H.add_edge(lig['Name'][0], ligcol.index[2])
                lig['Graph'] = H
                break

            # 1. Make SPT
            incomplete = True
            while incomplete:
                for pert, score in lig['pairs_dict'].items():
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
                for pert, score in lig['pairs_dict'].items():
                    l1, l2 = pert.split()[0], pert.split()[1]
                    if l1 in self.outer_nodes(H) or l2 in self.outer_nodes(H):
                        if (l1, l2) not in H.edges or (l2, l1) not in H.edges:
                            H.add_edge(l1, l2, weight=score)
                            break

            # 3. Add influence edges
            eig_cent = nx.eigenvector_centrality(H, max_iter=1000)
            eig_cent = {k: v for k, v in sorted(eig_cent.items(), key=lambda item: item[1], reverse=True)}
            cent_nodes = [k for k,v in nx.eigenvector_centrality(H, max_iter=1000).items() if v > 0.15]
            per_nodes = [k for k,v in nx.eigenvector_centrality(H, max_iter=1000).items() if v < 0.01]
            per_len = len(per_nodes)
            while per_len > 1:
                for pert, score in lig['pairs_dict'].items():
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

    def as_json(self):
        ## Return the nodes and edges of the graph as a Json string
        ##  The "keys" are lists of compatible keys for edges and nodes as read
        ##  at /networkgen/static/js/networkgen.js
        edge_keys = ["label", "freenrg", "sem", "crashes", "from", "to"]
        node_keys = ["shape", "label", "image", "id"]
        nodes = set([])
        result = {"nodes": [], "edges": []}

        for charge, lig in self.lig_dict.items():
            # Add unique nodes for this charge to the final nodes
            nodes = nodes | set(
                [node for edge in lig['Graph'].edges for node in edge])
            for edge in lig['Graph'].edges:
                result["edges"].append(
                    {"from": edge[0],
                     "to": edge[1]})

        for node in nodes:
            result["nodes"].append(
                {"label": node,
                 "image": "",  # TODO: networkgen.js requires an image. Make it
                               # use some kind of generic image.
                 "id": node})

        return json.dumps(result, indent=2)


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
                        default='MFP',
                        choices=['MFP', 'Tanimoto', 'MCS', 'SMILES'],
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
    print(mg.as_json())  # TODO: This gets printed, but should go into some
                         # model field.

if __name__ == "__main__":
    main()
