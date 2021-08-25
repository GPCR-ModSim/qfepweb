from pathlib import Path
from unittest import TestCase
import pytest

from . import mapgen

class Network:
    MCS = "MSC"
    SMILES = "SMILES"
    MFP = "MFP"
    Tanimoto = "Tanimoto"
    metric = SMILES


class MapGenTests(TestCase):
    @pytest.mark.timeout(2)
    def test_make_map_doesnt_infinite_loop(self):
        # This method enters a infinite loop in some cases.
        net = Network()
        sdf = Path(__file__).parent / "test_files" / "Two_Ligands.sdf"
        obj = mapgen.MapGen(network_obj=net, in_sdf=sdf)

        obj.make_map()

    def test_setup_of_ligands(self):
        net = Network()
        sdf = Path(__file__).parent / "test_files" / "CDK2_ligands.sdf"
        m = mapgen.MapGen(network_obj=net, in_sdf=sdf)

        m.make_map()

        assert list(m.ligands[0].keys()) == ["Name", "PoolIdx", "FP", "Scores", "Graph"]
        expected = {
            0: {1: {'weight': 48.5}, 11: {'weight': 44.35}},
            1: {6: {'weight': 49.4}, 0: {'weight': 48.5},
                10: {'weight': 48.5}, 5: {'weight': 46.5}},
            2: {6: {'weight': 43.25}, 11: {'weight': 43.0}},
            3: {9: {'weight': 39.5}, 15: {'weight': 38.5}},
            4: {9: {'weight': 38.05}, 8: {'weight': 36.95}},
            5: {10: {'weight': 47.0}, 1: {'weight': 46.5}},
            6: {1: {'weight': 49.4}, 2: {'weight': 43.25}},
            7: {8: {'weight': 39.0}, 12: {'weight': 39.0}, 13: {'weight': 39.0},
                14: {'weight': 37.9}},
            8: {9: {'weight': 39.0}, 7: {'weight': 39.0}, 12: {'weight': 39.0},
                13: {'weight': 39.0}, 14: {'weight': 37.9}, 4: {'weight': 36.95}},
            9: {11: {'weight': 39.9}, 3: {'weight': 39.5}, 8: {'weight': 39.0},
                4: {'weight': 38.05}, 15: {'weight': 38.45}},
            10: {1: {'weight': 48.5}, 5: {'weight': 47.0}},
            11: {0: {'weight': 44.35}, 9: {'weight': 39.9}, 2: {'weight': 43.0}},
            12: {7: {'weight': 39.0}, 8: {'weight': 39.0}},
            13: {7: {'weight': 39.0}, 8: {'weight': 39.0}},
            14: {7: {'weight': 37.9}, 8: {'weight': 37.9}},
            15: {3: {'weight': 38.5}, 9: {'weight': 38.45}}}

        self.assertEqual(m.ligands[0]["Graph"].adj, expected)

    def test_map_with_one_ligand(self):
        net = Network()

        sdf = Path(__file__).parent / "test_files" / "One_ligand.sdf"

        m = mapgen.MapGen(network_obj=net, in_sdf=sdf)

        m.make_map()

        assert m.ligands[0]["Graph"].adj == {'18629-1': {}}

    def test_map_with_two_ligand(self):
        net = Network()

        sdf = Path(__file__).parent / "test_files" / "Two_Ligands.sdf"

        m = mapgen.MapGen(network_obj=net, in_sdf=sdf)

        m.make_map()

        expected = {'30': {'28': {'weight': 48.5}},
                    '28': {'30': {'weight': 48.5}}}

        assert m.ligands[0]["Graph"].adj == expected
