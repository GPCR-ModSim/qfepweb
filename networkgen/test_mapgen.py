from pathlib import Path
from unittest import TestCase
import pytest

from . import mapgen


class MapGenTests(TestCase):
    @pytest.mark.timeout(2)
    def test_make_map_doesnt_infinite_loop(self):
        # This method enters a infinite loop in some cases.
        class Network:
            MCS = "MSC"
            SMILES = "SMILES"
            MFP = "MFP"
            Tanimoto = "Tanimoto"
            metric = SMILES

        net = Network()
        sdf = Path(__file__).parent / "test_files" / "Two_Ligands.sdf"
        obj = mapgen.MapGen(network_obj=net, in_sdf=sdf)

        obj.make_map()
