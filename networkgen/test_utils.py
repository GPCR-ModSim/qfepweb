import io
import json
from pathlib import Path
from tempfile import NamedTemporaryFile
from uuid import UUID

from django.test import TestCase
from model_bakery import baker
from diffimg import diff
from rdkit import Chem
from rdkit.DataStructs import FingerprintSimilarity

from networkgen import mapgen


class MapGenerator(TestCase):
    @staticmethod
    def loadSdf(path):
        """Load a SDF file and return a io.Bytes object."""
        with open(path, "rb") as f:
            return io.BytesIO(f.read())

    def setUp(self):
        self.sdfPath = Path(__file__).parent / "test_files" / "CDK2_ligands.sdf"
        self.sdfIo = self.loadSdf(self.sdfPath)
        self.genObj = baker.make_recipe("networkgen.network")

    def test_mapgen_can_be_inited_with_io_stream(self):
        m = mapgen.MapGen(network_obj=self.genObj)

        assert m.metric == self.genObj.MFP

    def test_mapgen_can_be_inited_with_string(self):
        """Passing a file path looks more human than passing the io around."""
        m = mapgen.MapGen(in_sdf=self.sdfPath, network_obj=self.genObj)

        assert m.metric == self.genObj.MFP

    def test_mapgen_setup_of_ligands(self):
        """After MapGen loads the SDF data, the ligand setup must be called."""
        #m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)
        m = mapgen.MapGen(network_obj=self.genObj)

        assert list(m.ligands.keys()) == [0]
        assert list(m.ligands[0].keys()) == ["Ligand", "Scores"]
        assert len(list(m.ligands[0]["Ligand"])) == 16  # Items in the file

    def test_set_the_similarity_function(self):
        self.genObj.metric = self.genObj.MCS
        m = mapgen.MapGen(network_obj=self.genObj)

        m._set_similarity_function()
        assert m.simF.__name__ == "FindMCS"

        for metric in [self.genObj.Tanimoto, self.genObj.MFP]:
            m.metric = metric
            m._set_similarity_function()
            assert m.simF.__name__ == "FingerprintSimilarity"

        m.metric = self.genObj.SMILES
        m._set_similarity_function()
        assert m.simF.__name__ == "globalms"

    def test_simmilarity_matrix(self):
        """Compute a similarity matrix for each ligand against all the others.

        This method is the meat of the class, so test thorougly. Any fails here
        and the whole web is doomed with bad input.
        """
        ## MFP
        m = mapgen.MapGen(network_obj=self.genObj)

        scores = m.ligands[0]["Scores"]
        assert len(scores) == sum(range(1, 16))
        # Lets check some precalculated values
        assert scores[(0, 1)] == 0.825
        assert scores[(0, 15)] == 0.727
        assert scores[(14, 15)] == 0.774

    def test_simmilarity_matrix_tanimoto(self):
        # Tanimoto
        self.genObj.metric = self.genObj.Tanimoto
        m = mapgen.MapGen(network_obj=self.genObj)

        scores = m.ligands[0]["Scores"]
        assert len(scores) == sum(range(1, 16))
        # Lets check some precalculated values
        assert scores[(0, 1)] == 0.877
        assert scores[(0, 15)] == 0.881
        assert scores[(14, 15)] == 0.943

    def test_simmilarity_matrix_smiles(self):
        """This test is like the above "test_simmilarity_matrix" but it takes
        forever."""
        # Smiles
        self.genObj.metric = self.genObj.SMILES
        m = mapgen.MapGen(network_obj=self.genObj)

        scores = m.ligands[0]["Scores"]
        assert len(scores) == sum(range(1, 16))
        # Lets check some precalculated values
        assert scores[(0, 1)] == 48.5
        assert scores[(0, 15)] == 38.05
        assert scores[(14, 15)] == 36.45

    def test_simmilarity_matrix_mcs(self):
        """This test is like the above "test_simmilarity_matrix" but it takes
        forever."""
        # MCS
        self.genObj.metric = self.genObj.MCS
        m = mapgen.MapGen(network_obj=self.genObj)

        scores = m.ligands[0]["Scores"]
        assert len(scores) == sum(range(1, 16))
        # Lets check some precalculated values
        assert scores[(0, 1)] == 59
        assert scores[(0, 15)] == 51
        assert scores[(14, 15)] == 51

    def test_mcs_calculation(self):
        m = mapgen.MapGen(network_obj=self.genObj)

        assert m.pool.mcs.smartsString == "[#6&R]1(:&@[#7&R]:&@[#6&R](-&!@[" +\
            "#8&!R]-&!@[#6&!R]-&!@[#6&R]2-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&" +\
            "R]-&@[#6&R]-&@2):&@[#6&R]2:&@[#6&R](:&@[#7&R]:&@1):&@[#7&R]:&@" +\
            "[#6&R]:&@[#7&R]:&@2)-&!@[#7&!R]-&!@[#6&R]1:&@[#6&R]:&@[#6&R]:&" +\
            "@[#6&R]:&@[#6&R]:&@[#6&R]:&@1"

    def test_network_loading(self):
        """The network should be loaded at the beginning, because when the
        self.suppl gets consumed it becomes a pain to workwith."""
        m = mapgen.MapGen(network_obj=self.genObj)

        assert len(m.pool) == 16
        assert m.pool[0].GetProp("_Name") == "30"
        assert [_.GetProp("_Name") for _ in m.pool] == [
            '30', '28', '1oiy', '1oi9', '32', '1oiu', '29', '1h1r', '21', '26',
            '1h1s', '31', '20', '22', '17', '1h1q']

    def test_work_with_raw_objects_no_db(self):
        self.sdfIo.seek(0)
        class DummyNetwork:
            MCS = self.genObj.MCS
            Tanimoto = self.genObj.Tanimoto
            MFP = self.genObj.MFP
            SMILES = self.genObj.SMILES

            in_sdf = self.sdfIo
            metric = self.genObj.SMILES

        obj = DummyNetwork()

        m = mapgen.MapGen(network_obj=obj)

        assert m.metric == obj.SMILES

        assert len(m.ligands[0]["Scores"]) == sum(range(1, 16))

class ImageGenerator(TestCase):
    """A class for creating molecule images."""
    def setUp(self):
        self.test_files = Path(__file__).parent / "test_files"
        with open(self.test_files / "CDK2_ligands.sdf", "rb") as m:
            mol_data = io.BytesIO(m.read())

        self.molecules = list(Chem.ForwardSDMolSupplier(mol_data))
        self.pool = mapgen.MoleculePool(self.molecules)

    def test_color_palette(self):
        ## This is the default palette, precalculated
        imgr = mapgen.MoleculeImage(pool_idx=0, pool=self.pool)
        assert imgr.palette == [
            (0, 0, 0),
            (0.9019607843137255, 0.6235294117647059, 0.0),
            (0.33725490196078434, 0.7058823529411765, 0.9137254901960784),
            (0.0, 0.6196078431372549, 0.45098039215686275),
            (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
            (0.0, 0.4470588235294118, 0.6980392156862745),
            (0.8352941176470589, 0.3686274509803922, 0.0),
            (0.8, 0.4745098039215686, 0.6549019607843137)]

        imgr = mapgen.MoleculeImage(pool_idx=0, pool=self.pool,
                                    palette="IBM")
        assert imgr.palette == [
            (0, 0, 0),
            (0.39215686274509803, 0.5607843137254902, 1.0),
            (0.47058823529411764, 0.3686274509803922, 0.9411764705882353),
            (0.8627450980392157, 0.14901960784313725, 0.4980392156862745),
            (0.996078431372549, 0.3803921568627451, 0.0),
            (1.0, 0.6901960784313725, 0.0)]

    def test_needed_properties_for_image_generation(self):
        imgr = mapgen.MoleculeImage(pool_idx=0,
                                    pool=self.pool)

        assert imgr.name == "30"

        assert isinstance(imgr.core, Chem.rdchem.Mol)

    def test_molecules_are_flatten(self):
        imgr = mapgen.MoleculeImage(pool_idx=0,
                                    pool=self.pool)
        assert imgr._flatten_molecule() == None

    def test_substructure_core_finder(self):
        """We have to find the MCS (core) structure in our molecule."""
        imgr = mapgen.MoleculeImage(pool_idx=0,
                                    pool=self.pool)

        assert [_.GetProp("SourceAtomIdx") for _ in imgr.molecule.GetAtoms()]

    def test_png_from_smiles_single_molecule(self):
        molecule = Chem.MolFromSmiles('Cc1nc(C)c(s1)c2ccnc(Nc3ccccc3F)n2')
        pool = mapgen.MoleculePool([molecule])
        imgr = mapgen.MoleculeImage(pool_idx=0, pool=pool)

        with open(self.test_files / "Plain.png", "rb") as r:
            assert imgr.png() == r.read()

    def test_png_generation(self):
        imgr = mapgen.MoleculeImage(pool_idx=0,
                                    pool=self.pool)

        with NamedTemporaryFile() as tmpimg:
            tmpimg.write(imgr.png())
            assert diff(self.test_files / "Sample30.png", tmpimg.name) < 0.01

    def test_png_from_smiles(self):
        mol1 = Chem.MolFromSmiles('Cc1nc(C)c(s1)c2ccnc(Nc3ccccc3F)n2')
        mol2 = Chem.MolFromSmiles('Cc1nc(Nc5ccccc5)c(s1)c2ccnc(Nc3ccccc3F)n2')
        pool = mapgen.MoleculePool([mol1, mol2])

        imgr = mapgen.MoleculeImage(pool_idx=1, pool=pool)
        imgr.name = "Sample"
        with NamedTemporaryFile() as tmpimg:
            tmpimg.write(imgr.png())
            assert diff(self.test_files / "HollowRing.png", tmpimg.name) < 0.01

        imgr = mapgen.MoleculeImage(pool_idx=1, pool=pool)
        imgr.fill_rings = True
        imgr.name = "Sample2"
        with NamedTemporaryFile() as tmpimg:
            tmpimg.write(imgr.png())
            assert diff(self.test_files / "FilledRing.png", tmpimg.name) < 0.01


class PoolGenerator(TestCase):
    def setUp(self):
        self.test_files = Path(__file__).parent / "test_files"
        with open(self.test_files / "CDK2_ligands.sdf", "rb") as m:
            mol_data = io.BytesIO(m.read())

        self.molecules = list(Chem.ForwardSDMolSupplier(mol_data))

    def test_the_cores_can_be_calculated(self):
        pool = mapgen.MoleculePool(self.molecules)

        assert isinstance(pool.mcs, Chem.rdFMCS.MCSResult)
        assert isinstance(pool.core, Chem.rdchem.Mol)
        assert isinstance(pool.query_core, Chem.rdchem.Mol)

    def test_group_decomposition(self):
        pool = mapgen.MoleculePool(self.molecules)

        assert len(pool) == 16
        for molecule in pool.groups:
            assert list(molecule.keys()) == ["Core", "R1", "R2"]

    def test_can_add_molecules_to_pool(self):
        pool = mapgen.MoleculePool(self.molecules)

        assert len(pool) == 16

        pool.append(self.molecules[0])

        assert len(pool) == 17

    def test_can_grow_pool(self):
        pool = mapgen.MoleculePool()

        assert len(pool) == 0
        assert len(pool.groups) == 0

        for m in self.molecules:
            pool.append(m)

        assert len(pool) == 16
        assert len(pool.groups) == 16
