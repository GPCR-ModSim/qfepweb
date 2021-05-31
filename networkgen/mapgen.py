import argparse
import io
import json
import networkx as nx
import pandas as pd
from networkgen.models import Generator as g 
from networkgen.models import Ligand
from collections import defaultdict
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import rdqueries
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Geometry
rdDepictor.SetPreferCoordGen(True)


class MapGen():
    def __init__(self, in_sdf, metric, network_obj=None):
        in_sdf.seek(0)
        self.suppl = Chem.ForwardSDMolSupplier(in_sdf)
        self.metric = metric
        self.network = network_obj
        self.lig_dict = {}
        self.sim_dfs = {}

    def make_fp(self, mol):
        if self.metric == g.Tanimoto:
            fp = FingerprintMols.FingerprintMol(mol)
        elif self.metric == g.MFP:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        elif self.metric == g.SMILES:
            fp = Chem.MolToSmiles(mol, isomericSmiles=True)
        else: fp = None
        return fp
    
    def get_ligand(self, rdmol):
        n = rdmol.GetProp('_Name')
        return Ligand(charge=Chem.rdmolops.GetFormalCharge(rdmol), atom_number=len([atom for atom in rdmol.GetAtoms()]), name=n, SMILES=Chem.MolToSmiles(rdmol, isomericSmiles=True), image='networkgen/molimages/{}.png'.format(n), network=self.network)

    def generateImages(self, mol, row, core, width=350, height=200,
                        fillRings=False,legend="",
                        sourceIdxProperty="SourceAtomIdx",
                        lbls=None):

        # copy the molecule and core
        mol = Chem.Mol(mol)
        core = Chem.Mol(core)

        # -------------------------------------------
        # include the atom map numbers in the substructure search in order to 
        # try to ensure a good alignment of the molecule to symmetric cores
        for at in core.GetAtoms():
            if at.GetAtomMapNum():
                at.ExpandQuery(rdqueries.IsotopeEqualsQueryAtom(200+at.GetAtomMapNum()))
                
        for lbl in row:
            if lbl== 'Core':
                continue
            rg = row[lbl]
            for at in rg.GetAtoms():
                if not at.GetAtomicNum() and at.GetAtomMapNum() and \
                at.HasProp('dummyLabel') and at.GetProp('dummyLabel')==lbl:
                    # attachment point. the atoms connected to this
                    # should be from the molecule
                    for nbr in at.GetNeighbors():
                        if nbr.HasProp(sourceIdxProperty):
                            mAt = mol.GetAtomWithIdx(nbr.GetIntProp(sourceIdxProperty))
                            if mAt.GetIsotope():
                                mAt.SetIntProp('_OrigIsotope',mAt.GetIsotope())
                            mAt.SetIsotope(200+at.GetAtomMapNum())

        # remove unmapped hs so that they don't mess up the depiction
        rhps = Chem.RemoveHsParameters()
        rhps.removeMapped = False
        tmol = Chem.RemoveHs(mol,rhps)
        rdDepictor.GenerateDepictionMatching2DStructure(tmol,core)

        oldNewAtomMap={}
        # reset the original isotope values and account for the fact that
        # removing the Hs changed atom indices
        for i,at in enumerate(tmol.GetAtoms()):
            if at.HasProp(sourceIdxProperty):
                oldNewAtomMap[at.GetIntProp(sourceIdxProperty)] = i
                if at.HasProp("_OrigIsotope"):
                    at.SetIsotope(at.GetIntProp("_OrigIsotope"))
                    at.ClearProp("_OrigIsotope")
                else:
                    at.SetIsotope(0)
        
        # ------------------
        #  set up our colormap
        #   the three choices here are all "colorblind" colormaps
        
        # "Tol" colormap from https://davidmathlogic.com/colorblind
        colors = [(51,34,136),(17,119,51),(68,170,153),(136,204,238),(221,204,119),(204,102,119),(170,68,153),(136,34,85)]
        # "IBM" colormap from https://davidmathlogic.com/colorblind
        colors = [(100,143,255),(120,94,240),(220,38,127),(254,97,0),(255,176,0)]
        # Okabe_Ito colormap from https://jfly.uni-koeln.de/color/
        colors = [(230,159,0),(86,180,233),(0,158,115),(240,228,66),(0,114,178),(213,94,0),(204,121,167)]
        for i,x in enumerate(colors):
            colors[i] = tuple(y/255 for y in x)
    
        #----------------------
        # Identify and store which atoms, bonds, and rings we'll be highlighting
        highlightatoms = defaultdict(list)
        highlightbonds = defaultdict(list)
        atomrads = {}
        widthmults = {}

        rings = []
        # loop over R groups.
        for i,lbl in enumerate(lbls):    
            color = colors[i%len(colors)]
            try:
                rquery = row[lbl]
            # we don't know the number of R-groups, so just quit this loop if we've reached the end.
            except KeyError:
                continue

            Chem.GetSSSR(rquery)
            rinfo = rquery.GetRingInfo()
            for at in rquery.GetAtoms():
                if at.HasProp(sourceIdxProperty):
                    origIdx = oldNewAtomMap[at.GetIntProp(sourceIdxProperty)]
                    highlightatoms[origIdx].append(color)
                    atomrads[origIdx] = 0.4
            if fillRings:
                for aring in rinfo.AtomRings():
                    tring = []
                    allFound = True
                    for aid in aring:
                        at = rquery.GetAtomWithIdx(aid)
                        if not at.HasProp(sourceIdxProperty):
                            allFound = False
                            break
                        tring.append(oldNewAtomMap[at.GetIntProp(sourceIdxProperty)])
                    if allFound:
                        rings.append((tring,color))
            for qbnd in rquery.GetBonds():
                batom = qbnd.GetBeginAtom()
                eatom = qbnd.GetEndAtom()
                if batom.HasProp(sourceIdxProperty) and eatom.HasProp(sourceIdxProperty):
                    origBnd = tmol.GetBondBetweenAtoms(oldNewAtomMap[batom.GetIntProp(sourceIdxProperty)],
                                                    oldNewAtomMap[eatom.GetIntProp(sourceIdxProperty)])
                    bndIdx = origBnd.GetIdx()
                    highlightbonds[bndIdx].append(color)
                    widthmults[bndIdx] = 2

        d2d = rdMolDraw2D.MolDraw2DCairo(width,height)
        dos = d2d.drawOptions()
        dos.useBWAtomPalette()
                    
        #----------------------
        # if we are filling rings, go ahead and do that first so that we draw
        # the molecule on top of the filled rings
        if fillRings and rings:
            # a hack to set the molecule scale
            d2d.DrawMoleculeWithHighlights(tmol,legend,dict(highlightatoms),
                                        dict(highlightbonds),
                                        atomrads,widthmults)
            d2d.ClearDrawing()
            conf = tmol.GetConformer()
            for (aring,color) in rings:
                ps = []
                for aidx in aring:
                    pos = Geometry.Point2D(conf.GetAtomPosition(aidx))
                    ps.append(pos)
                d2d.SetFillPolys(True)
                d2d.SetColour(color)
                d2d.DrawPolygon(ps)
            dos.clearBackground = False

        #----------------------
        # now draw the molecule, with highlights:
        d2d.DrawMoleculeWithHighlights(tmol,"",dict(highlightatoms),dict(highlightbonds),
                                    atomrads,widthmults)
        d2d.FinishDrawing()
        png = d2d.GetDrawingText()

        # save to file:
        d2d.WriteDrawingText("static/molimages/{}.png".format(legend))
            
        return png

    def writeLigandImages(self):
        """
        Takes a list of ligand names and generates .png images in ./data/mol_images/ using RDKit.
        Workflow mostly based on https://rdkit.blogspot.com/2020/10/molecule-highlighting-and-r-group.html
        """

        lignames = [self.lig_dict[k]['Name'] for k in self.lig_dict.keys()]
        lignames = [lig for charge_group in lignames for lig in charge_group] #flatten ligand names

        mols = [self.lig_dict[k]['Mol'] for k in self.lig_dict.keys()]
        mols = [m for charge_group in mols for m in charge_group] #flatten mols

        if len(lignames) < 2:
            raise Exception("Number of ligands ({}) should be > 1".format(len(lignames)))

        ################## RDKIT IMAGE GENERATION ################
        # flatten all ligands into 2D space. Note that most rdkit mol operations are in-place.
        for m in mols:
            rdDepictor.Compute2DCoords(m)
            m.UpdatePropertyCache()


        # find the core using standard MCS algorithm. Flatten the core again just to be sure.
        mcs = rdFMCS.FindMCS(mols, matchValences=False,
                                        ringMatchesRingOnly=True,
                                        completeRingsOnly=True,
                                        matchChiralTag=False)

        core = Chem.MolFromSmarts(mcs.smartsString)
        rdDepictor.Compute2DCoords(core)

        # find subtructure matches with MCS per ligand, then tag matching atom indices in each ligand.
        ps = Chem.AdjustQueryParameters.NoAdjustments()
        ps.makeDummiesQueries=True
        qcore = Chem.AdjustQueryProperties(core,ps)
        #mhs = [Chem.AddHs(x,addCoords=True) for x in ms]
        mms = [x for x in mols if x.HasSubstructMatch(qcore)]
        for m in mms:
            for atom in m.GetAtoms():
                atom.SetIntProp("SourceAtomIdx",atom.GetIdx())

        # do an RDKit R-group decomposition.	
        groups,_ = rdRGroupDecomposition.RGroupDecompose([qcore],mms,asSmiles=False,asRows=True)

        # call the writer function with each molecule.
        for i, m in enumerate(mols):
            png = self.generateImages(m,groups[i],qcore,lbls=('R1','R2','R3','R4', 'R5', 'R6', 'R7'),
                                    legend=lignames[i],
                                    width=400,height=400)

    def set_ligdict(self):
        for mol in self.suppl:
            charge = Chem.rdmolops.GetFormalCharge(mol)
            ligand = self.get_ligand(mol)
            ligand.save()
            v = self.lig_dict.setdefault(charge, {'Name': [], 'Mol': [], 'FP': []})
            v['Name'].append(mol.GetProp('_Name'))
            v['Mol'].append(mol)
            if self.metric != g.MCS:
                v['FP'].append(self.make_fp(mol))
        self.writeLigandImages()

    def sim_mx(self):
        if self.metric in [g.Tanimoto, g.MFP]:
            from rdkit import DataStructs
        elif self.metric == g.MCS:
            from rdkit.Chem import rdFMCS
        elif self.metric == g.SMILES:
            from Bio import pairwise2
        for charge, item in self.lig_dict.items():
            df = pd.DataFrame()
            for index, i in enumerate(item['Name']):
                for jndex, j in enumerate(item['Name']):
                    if i == j:
                        df.loc[i, j] = 1.0
                    else:
                        if self.metric in [g.Tanimoto, g.MFP]:
                            df.loc[i, j] = DataStructs.FingerprintSimilarity(
                                item['FP'][index], item['FP'][jndex])
                        if self.metric == g.MCS:
                            mcs = rdFMCS.FindMCS(
                                [item['Mol'][index],
                                 item['Mol'][jndex]],
                                atomCompare=rdFMCS.AtomCompare.CompareAny)
                            df.loc[i, j] = mcs.numAtoms + mcs.numBonds
                        if self.metric == g.SMILES:
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

            if self.metric in [g.Tanimoto, g.MFP]:
                self.lig_dict[charge]['df'] = 1 - self.lig_dict[charge]['df']  # get dissimilarity matrix

            if self.metric == g.MCS:
                self.lig_dict[charge]['df'] = 100 - self.lig_dict[charge]['df']  # get dissimilarity matrix
            self.lig_dict[charge]['df'] = self.lig_dict[charge]['df'].replace(1.0, 0.0)  # set diagonal to 0
            self.lig_dict[charge]['df'] = self.lig_dict[charge]['df'].replace(0.0, 1.0)  # set zeroes into 1 (in order to search shortest path)

    def set_ligpairs(self):
        for charge in self.lig_dict.keys():
            pairs_dict = {}
            for i in self.lig_dict[charge]['df'].index:
                for j in self.lig_dict[charge]['df'].columns:
                    if self.lig_dict[charge]['df'].loc[i, j] != 1.0:
                        pairs_dict['{} {}'.format(i, j)] = round(self.lig_dict[charge]['df'].loc[i, j], 3)

            if self.metric in [g.Tanimoto, g.MFP, g.MCS]:
                pairs_dict = {k: v for k, v in sorted(pairs_dict.items(), key=lambda item: item[1], reverse=False)}
            elif self.metric == g.SMILES:
                pairs_dict = {k: v for k, v in sorted(pairs_dict.items(), key=lambda item: item[1], reverse=True)}

            self.lig_dict[charge]['pairs_dict'] = pairs_dict

    def process_map(self):
        self.set_ligdict()
        self.sim_mx()
        self.clean_mxs()
        self.set_ligpairs()
        self.make_map()
        self.as_json()

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
                 "image": "/static/molimages/{}.png".format(node),
                 "id": node})

        self.network.network = json.dumps(result)
        with open('networkgen/templates/networkgen/graph.json', 'w') as outfile:
            json.dump(result, outfile)
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
    print(mg.as_json())  # TODO: This gets printed, but should go into some
                         # model field.

if __name__ == "__main__":
    main()
