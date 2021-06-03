import io
from pathlib import Path

from django.test import override_settings, TestCase
from model_bakery import baker
import pytest
from rdkit import Chem
from rdkit.DataStructs import FingerprintSimilarity

from networkgen import mapgen
from networkgen.models import Generator as g, Ligand


class MapGenerator(TestCase):
    @staticmethod
    def loadSdf(path):
        """Load a SDF file and return a io.Bytes object."""
        with open(path, "rb") as f:
            return io.BytesIO(f.read())

    def setUp(self):
        self.sdfPath = Path(__file__).parent / "test_files" / "CDK2_ligands.sdf"
        self.sdfIo = self.loadSdf(self.sdfPath)
        self.genObj = baker.make("Generator")
        self.img_dir = Path(__file__).parent / "test_files" / "media"
        self.img_dir.mkdir(exist_ok=True)

    def tearDown(self):
        for f in Path(self.img_dir).iterdir():
            try:
                f.unlink()
            except IsADirectoryError:
                for img in f.iterdir():
                    img.unlink()
                f.rmdir()

        Path(self.img_dir).rmdir()

    def test_mapgen_can_be_inited_with_io_stream(self):
        m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)

        assert m.metric == g.MFP

    def test_mapgen_can_be_inited_with_string(self):
        """Passing a file path looks more human than passing the io around."""
        m = mapgen.MapGen(in_sdf=self.sdfPath, network_obj=self.genObj)

        assert m.metric == g.MFP

    def test_mapgen_setup_of_ligands(self):
        """After MapGen loads the SDF data, the ligand setup must be called."""
        m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)

        assert list(m.ligands.keys()) == [0]
        assert list(m.ligands[0].keys()) == ["Name", "Mol", "FP", "df"]
        assert len(list(m.ligands[0]["Name"])) == 16  # Items in the file

    def test_set_the_similarity_function(self):
        self.genObj.metric = g.MCS
        m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)

        m._set_similarity_function()
        assert m.simF.__name__ == "FindMCS"

        for metric in [g.Tanimoto, g.MFP]:
            m.metric = metric
            m._set_similarity_function()
            assert m.simF.__name__ == "FingerprintSimilarity"

        m.metric = g.SMILES
        m._set_similarity_function()
        assert m.simF.__name__ == "globalms"

    def test_simmilarity_matrix(self):
        """Compute a similarity matrix for each ligand against all the others.

        This method is the meat of the class, so test thorougly. Any fails here
        and the whole web is doomed with bad input.
        """
        ## MFP
        m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)

        matrix = m.ligands[0]["df"]
        assert matrix.shape == (16, 16)  # A matrix lig x lig
        # Only the down triangle is calculated
        assert matrix.iloc[1, 0] > 0
        assert matrix.iloc[1, 0] != matrix.iloc[0, 1]
        # Lets check some precalculated values
        assert matrix.loc["28", "30"] == 0.8245614035087719
        assert matrix.loc["1h1q", "30"] == 0.7272727272727273
        assert matrix.loc["1h1q", "17"] == 0.7735849056603774

    def test_simmilarity_matrix_tanimoto(self):
        # Tanimoto
        self.sdfIo.seek(0)
        self.genObj.metric = g.Tanimoto
        m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)

        matrix = m.ligands[0]["df"]
        assert matrix.shape == (16, 16)
        # Only the down triangle is calculated
        assert matrix.iloc[1, 0] > 0
        assert matrix.iloc[1, 0] != matrix.iloc[0, 1]
        # Lets check some precalculated values
        assert matrix.loc["28", "30"] == 0.8771367521367521
        assert matrix.loc["1h1q", "30"] == 0.8807870370370371
        assert matrix.loc["1h1q", "17"] == 0.942998760842627

    def test_simmilarity_matrix_smiles(self):
        """This test is like the above "test_simmilarity_matrix" but it takes
        forever."""
        # Smiles
        self.sdfIo.seek(0)
        self.genObj.metric = g.SMILES
        m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)

        matrix = m.ligands[0]["df"]
        assert matrix.shape == (16, 16)
        # Only the down triangle is calculated
        assert matrix.iloc[1, 0] > 0
        assert matrix.iloc[1, 0] != matrix.iloc[0, 1]
        # Diagonal for this method is 100.0, not 1.0
        assert matrix.iloc[0, 0] == 100.0
        # Lets check some precalculated values
        assert matrix.loc["28", "30"] == 48.5
        assert matrix.loc["1h1q", "30"] == 38.05
        assert matrix.loc["1h1q", "17"] == 36.45

    def test_simmilarity_matrix_mcs(self):
        """This test is like the above "test_simmilarity_matrix" but it takes
        forever."""
        # MCS
        self.sdfIo.seek(0)
        self.genObj.metric = g.MCS
        m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)

        matrix = m.ligands[0]["df"]
        assert matrix.shape == (16, 16)
        # Only the down triangle is calculated
        assert matrix.iloc[1, 0] > 0
        assert matrix.iloc[1, 0] != matrix.iloc[0, 1]
        # Diagonal for this method is 100.0, not 1.0
        assert matrix.iloc[0, 0] == 100.0
        # Lets check some precalculated values
        assert matrix.loc["28", "30"] == 59.0
        assert matrix.loc["1h1q", "30"] == 51.0
        assert matrix.loc["1h1q", "17"] == 51.0

    def test_mcs_calculation(self):
        self.sdfIo.seek(0)
        m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)

        assert m.mcs.smartsString == "[#6&R]1(:&@[#7&R]:&@[#6&R](-&!@[#8&!R]" +\
            "-&!@[#6&!R]-&!@[#6&R]2-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&R]-" +\
            "&@[#6&R]-&@2):&@[#6&R]2:&@[#6&R](:&@[#7&R]:&@1):&@[#7&R]:" +\
            "&@[#6&R]:&@[#7&R]:&@2)-&!@[#7&!R]-&!@[#6&R]1:&@[#6&R]:&@[#6&R]:" +\
            "&@[#6&R]:&@[#6&R]:&@[#6&R]:&@1"

    def test_network_loading(self):
        """The network should be loaded at the beginning, because when the
        self.suppl gets consumed it becomes a pain to workwith."""
        self.sdfIo.seek(0)
        m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)

        assert len(m.molecules) == 16
        assert m.molecules[0].GetProp("_Name") == "30"
        assert [_.GetProp("_Name") for _ in m.molecules] == [
            '30', '28', '1oiy', '1oi9', '32', '1oiu', '29', '1h1r', '21', '26',
            '1h1s', '31', '20', '22', '17', '1h1q']

    @override_settings(MEDIA_ROOT=Path(__file__).parent / "test_files" / "media")
    def test_network_image_builder(self):
        """Test that the class can build the images for the objects."""
        m = mapgen.MapGen(in_sdf=self.sdfIo, network_obj=self.genObj)

        assert Ligand.objects.count() == 0
        ligands = m.save_ligands()
        assert Ligand.objects.count() == 16

        assert ligands[0].image.width == 400

class ImageGenerator(TestCase):
    """A class for creating molecule images."""
    def setUp(self):
        self.test_files = Path(__file__).parent / "test_files"
        with open(self.test_files / "CDK2_ligands.sdf", "rb") as m:
            mol_data = io.BytesIO(m.read())

        self.molecules = list(Chem.ForwardSDMolSupplier(mol_data))
        self.core = Chem.rdFMCS.FindMCS(self.molecules,
                                        matchValences=False,
                                        ringMatchesRingOnly=True,
                                        completeRingsOnly=True,
                                        matchChiralTag=False)

    def test_color_palette(self):
        ## This is the default palette, precalculated
        imgr = mapgen.MoleculeImage(self.molecules[0])
        assert imgr.palette == [
            (0, 0, 0),
            (0.9019607843137255, 0.6235294117647059, 0.0),
            (0.33725490196078434, 0.7058823529411765, 0.9137254901960784),
            (0.0, 0.6196078431372549, 0.45098039215686275),
            (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
            (0.0, 0.4470588235294118, 0.6980392156862745),
            (0.8352941176470589, 0.3686274509803922, 0.0),
            (0.8, 0.4745098039215686, 0.6549019607843137)]

        imgr = mapgen.MoleculeImage(self.molecules[0], palette="IBM")
        assert imgr.palette == [
            (0, 0, 0),
            (0.39215686274509803, 0.5607843137254902, 1.0),
            (0.47058823529411764, 0.3686274509803922, 0.9411764705882353),
            (0.8627450980392157, 0.14901960784313725, 0.4980392156862745),
            (0.996078431372549, 0.3803921568627451, 0.0),
            (1.0, 0.6901960784313725, 0.0)]

    def test_needed_properties_for_image_generation(self):
        imgr = mapgen.MoleculeImage(molecule=self.molecules[0],
                                    core=self.core)

        assert imgr.name == "30"

        assert isinstance(imgr.core, Chem.rdchem.Mol)

    def test_molecules_are_flatten(self):
        imgr = mapgen.MoleculeImage(molecule=self.molecules[0],
                                    core=self.core)
        assert imgr._flatten_molecule() == None

    def test_substructure_core_finder(self):
        """We have to find the MCS (core) structure in our molecule."""
        imgr = mapgen.MoleculeImage(molecule=self.molecules[0],
                                    core=self.core)

        with self.assertRaises(KeyError):
            # The key as not been set yet
            [_.GetProp("SourceAtomIdx") for _ in imgr.molecule.GetAtoms()]

        assert imgr._find_core() == None

        assert [_.GetProp("SourceAtomIdx") for _ in imgr.molecule.GetAtoms()]

    def test_png_from_smiles_single_molecule(self):
        molecule = Chem.MolFromSmiles('Cc1nc(C)c(s1)c2ccnc(Nc3ccccc3F)n2')
        imgr = mapgen.MoleculeImage(molecule=molecule)

        with open(self.test_files / "Plain.png", "rb") as r:
            assert imgr.png() == r.read()

    def test_png_generation(self):
        imgr = mapgen.MoleculeImage(molecule=self.molecules[0],
                                    core=self.core)

        with open(self.test_files / "Sample30.png", "rb") as r:
            assert imgr.png() == r.read()

    def test_png_from_smiles(self):
        mol1 = Chem.MolFromSmiles('Cc1nc(C)c(s1)c2ccnc(Nc3ccccc3F)n2')
        mol2 = Chem.MolFromSmiles('Cc1nc(Nc5ccccc5)c(s1)c2ccnc(Nc3ccccc3F)n2')

        core = Chem.rdFMCS.FindMCS([mol1, mol2],
                                   matchValences=False,
                                   ringMatchesRingOnly=True,
                                   completeRingsOnly=True,
                                   matchChiralTag=False)

        imgr = mapgen.MoleculeImage(molecule=mol2, core=core)
        imgr.name = "Sample"
        with open(self.test_files / "HollowRing.png", "rb") as r:
            assert imgr.png() == r.read()

        imgr = mapgen.MoleculeImage(molecule=mol2, core=core)
        imgr.fill_rings = True
        imgr.name = "Sample2"
        with open(self.test_files / "FilledRing.png", "rb") as r:
            assert imgr.png() == r.read()
