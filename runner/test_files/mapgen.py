import argparse
import io
import json
import os
import networkx as nx
import numpy as np
import pandas as pd
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from rdkit import Chem


class MapGen():
    def __init__(self, in_sdf, metric):
        self.suppl = Chem.ForwardSDMolSupplier(in_sdf)
        self.metric = metric.lower()
        self.lig_dict = {}
        self.sim_dfs = {}

    def make_fp(self, mol):
        if self.metric == 'tanimoto':
            fp = FingerprintMols.FingerprintMol(mol)
        if self.metric == 'mfp':
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        if self.metric == 'smiles':
            fp = Chem.MolToSmiles(mol, isomericSmiles=True)
        return fp

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
        if self.metric in ['tanimoto', 'mfp']:
            from rdkit import DataStructs
        elif self.metric == 'mcs':
            from rdkit.Chem import rdFMCS
        elif self.metric == 'smiles':
            from Bio import pairwise2
        for charge, item in self.lig_dict.items():
            df = pd.DataFrame()
            for index, i in enumerate(item['Name']):
                for jndex, j in enumerate(item['Name']):
                    if i == j:
                        df.loc[i, j] = 1.0
                    else:
                        if self.metric in ['tanimoto', 'mfp']:
                            df.loc[i, j] = DataStructs.FingerprintSimilarity(
                                item['FP'][index], item['FP'][jndex])
                        if self.metric == 'mcs':
                            mcs = rdFMCS.FindMCS(
                                [item['Mol'][index],
                                 item['Mol'][jndex]],
                                atomCompare=rdFMCS.AtomCompare.CompareAny)
                            df.loc[i, j] = mcs.numAtoms + mcs.numBonds
                        if self.metric == 'smiles':
                            alignments = pairwise2.align.globalms(
                                item['FP'][index], item['FP'][jndex], 1, -1, -0.5, -0.05)
                            df.loc[i, j] = alignments[0][2]
            self.lig_dict[charge]['df'] = df

    def clean_mxs(self):
        for charge in self.lig_dict.keys():
            for i, j in zip(self.lig_dict[charge]['df'].index, self.lig_dict[charge]['df'].idxmax()):
                vlist = self.lig_dict[charge]['df'].loc[i, :].tolist()
                index = vlist.index(1.0)
                vlist[index:] = [0.0] * len(vlist[index:])
                self.lig_dict[charge]['df'].loc[i, :] = vlist

            if self.metric in ['tanimoto', 'mfp']:
                df = 1 - self.lig_dict[charge]['df']  # get dissimilarity matrix

            if self.metric == 'mcs':
                df = 100 - self.lig_dict[charge]['df']  # get dissimilarity matrix
            self.lig_dict[charge]['df'] = df.replace(1.0, 0.0)  # set diagonal to 0
            self.lig_dict[charge]['df'] = df.replace(0.0, 1.0)  # set zeroes into 1 (in order to search shortest path)

    def set_ligpairs(self):
        for charge in self.lig_dict.keys():
            pairs_dict = {}
            for i in self.lig_dict[charge]['df'].index:
                for j in self.lig_dict[charge]['df'].columns:
                    if self.lig_dict[charge]['df'].loc[i, j] != 1.0:
                        pairs_dict['{} {}'.format(i, j)] = round(self.lig_dict[charge]['df'].loc[i, j], 3)

            if self.metric in ['tanimoto', 'mfp', 'mcs']:
                pairs_dict = {k: v for k, v in sorted(pairs_dict.items(), key=lambda item: item[1], reverse=False)}
            elif self.metric == 'smiles':
                pairs_dict = {k: v for k, v in sorted(pairs_dict.items(), key=lambda item: item[1], reverse=True)}

            self.lig_dict[charge]['pairs_dict'] = pairs_dict

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
                ligcol = self.sim_dfs[lig['Name']].sort_values(by=[lig['Name'][0]]) #complete similarity matrix
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
