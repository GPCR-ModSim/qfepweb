import io
import json
from pathlib import Path
from uuid import UUID

from django.conf import settings
from django.core.files import File
from django.test import override_settings, TestCase
from model_bakery import baker

from networkgen.models import Generator, Ligand


@override_settings(MEDIA_ROOT=Path(__file__).parent / "test_files" / "media")
class GeneratorModel(TestCase):
    def setUp(self):
        Path(settings.MEDIA_ROOT).mkdir(exist_ok=True)
        self.sdf_path = Path(__file__).parent / "test_files" / "CDK2_ligands.sdf"
        self.in_sdf = File(open(self.sdf_path, "rb"), name="CDK2_ligands.sdf")

    def tearDown(self):
        for l in Ligand.objects.all():
            Path(l.image.path).unlink(missing_ok=True)
        for g in Generator.objects.all():
                Path(g.in_sdf.path).unlink(missing_ok=True)
        try:
            Path(settings.MEDIA_ROOT / "molimages").rmdir()
            Path(settings.MEDIA_ROOT).rmdir()
        except OSError:
            pass

    def test_object_name(self):
        network = baker.make_recipe("networkgen.network")
        assert str(network) == f"Network Generator <{network.uuid}>"

    def test_network_json_is_created_on_saving(self):
        g = Generator(metric=Generator.MFP, in_sdf=self.in_sdf)
        assert g.network is None
        g.save()
        network_json = json.loads(g.network)

        assert list(network_json.keys()) == ["nodes", "edges"]

        ligands = [(_.name, _.uuid, _.image.url) for _ in g.ligand_set.all()]
        for node in network_json["nodes"]:
            assert (node["label"], UUID(node["id"]), node["image"]) in ligands

        uuids = [str(_) for _ in Ligand.objects.values_list("uuid", flat=True)]
        for edge in network_json["edges"]:
            assert edge["from"] in uuids
            assert edge["to"] in uuids

    def test_network_image_builder(self):
        assert Ligand.objects.count() == 0

        network = baker.make_recipe("networkgen.network")

        assert Ligand.objects.count() == 16

        for l in Ligand.objects.all():
            assert l.image.width == 400
            assert l.image.width == 400
