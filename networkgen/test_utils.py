import io
from pathlib import Path

from django.test import TestCase
from rdkit.DataStructs import FingerprintSimilarity
from model_bakery import baker

from networkgen import mapgen
from networkgen.models import Generator as g


class MapGenerator(TestCase):
    @staticmethod
    def loadSdf(path):
        """Load a SDF file and return a io.Bytes object."""
        with open(path, "rb") as f:
            return io.BytesIO(f.read())

    def setUp(self):
        self.sdfPath = Path(__file__).parent / "test_files/CDK2_ligands.sdf"
        self.sdfIo = self.loadSdf(self.sdfPath)
        self.genObj = baker.make("Generator")

    def test_mapgen_can_be_inited_with_io_stream(self):
        m = mapgen.MapGen(in_sdf=self.sdfIo, metric="mfp")

        assert m.metric == "mfp"

    def test_mapgen_can_be_inited_with_string(self):
        """Passing a file path looks more human than passing the io around."""
        m = mapgen.MapGen(in_sdf=self.sdfPath, metric=g.MFP)

        assert m.metric == g.MFP

    def test_mapgen_setup_of_ligands(self):
        """After MapGen loads the SDF data, the ligand setup must be called."""
        m = mapgen.MapGen(in_sdf=self.sdfIo, metric=g.MFP,
                          network_obj=self.genObj)
        assert m.lig_dict == {}

        m.set_ligdict()

        assert list(m.lig_dict.keys()) == [0]
        assert list(m.lig_dict[0].keys()) == ["Name", "Mol", "FP"]
        assert len(list(m.lig_dict[0]["Name"])) == 16  # Items in the file

    def test_set_the_similarity_function(self):
        m = mapgen.MapGen(in_sdf=self.sdfIo, metric=g.MCS,
                          network_obj=self.genObj)

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
        m = mapgen.MapGen(in_sdf=self.sdfIo, metric=g.MFP,
                          network_obj=self.genObj)

        assert m.lig_dict == {}
        m.sim_mx()

        matrix = m.lig_dict[0]["df"]
        assert matrix.shape == (16, 16)  # A matrix lig x lig
        # Only the down triangle is calculated
        assert matrix.iloc[1, 0] > 0
        assert matrix.iloc[1, 0] != matrix.iloc[0, 1]
        # Lets check some precalculated values
        assert matrix.loc["28", "30"] == 0.8245614035087719
        assert matrix.loc["1h1q", "30"] == 0.7272727272727273
        assert matrix.loc["1h1q", "17"] == 0.7735849056603774

        # Tanimoto
        self.sdfIo.seek(0)
        m = mapgen.MapGen(in_sdf=self.sdfIo, metric=g.Tanimoto,
                          network_obj=self.genObj)

        m.sim_mx()

        matrix = m.lig_dict[0]["df"]
        assert matrix.shape == (16, 16)
        # Only the down triangle is calculated
        assert matrix.iloc[1, 0] > 0
        assert matrix.iloc[1, 0] != matrix.iloc[0, 1]
        # Lets check some precalculated values
        assert matrix.loc["28", "30"] == 0.8771367521367521
        assert matrix.loc["1h1q", "30"] == 0.8807870370370371
        assert matrix.loc["1h1q", "17"] == 0.942998760842627

        # Smiles
        self.sdfIo.seek(0)
        m = mapgen.MapGen(in_sdf=self.sdfIo, metric=g.SMILES,
                          network_obj=self.genObj)

        m.sim_mx()

        matrix = m.lig_dict[0]["df"]
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

        # MCS
        self.sdfIo.seek(0)
        m = mapgen.MapGen(in_sdf=self.sdfIo, metric=g.MCS,
                          network_obj=self.genObj)

        m.sim_mx()

        matrix = m.lig_dict[0]["df"]
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
