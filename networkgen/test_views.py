"""The basic views for this project."""
from io import BytesIO
from pathlib import Path
import json
from django.template.response import TemplateResponse
from django.test import Client, TestCase
from django.urls import reverse
from model_bakery import baker

from networkgen.forms import GeneratorForm
from networkgen.models import Generator as g, Ligand


class NetworkGen(TestCase):
    """Test the landing of the generator pages.

    It should be available to everyone, without account, and a brief
    introdution to what the server does.

    """
    def setUp(self):
        self.network = baker.make_recipe("networkgen.network")
        self.ligands = baker.make("Ligand",
                                  network=self.network,
                                  _quantity=5)

    def test_the_home_is_available(self):
        page = self.client.get(reverse('networkgen:create'))

        assert page.status_code == 200
        assert "FEP Network" in page.content.decode()  # The title

    def test_the_javascript_and_css_is_included(self):
        page = self.client.get(
            reverse('networkgen:detail', kwargs={"pk": self.network.pk})
        ).content.decode()

        # There are proxies to identify the static assets needed
        requiredJs = ["lodash.min.js", "vis.js", "networkgen.js"]
        requiredCss = ["vis-network.min.css", "use.fontawesome.com",
                       "networkgen.css"]

        for asset in requiredJs + requiredCss:
            assert asset in page

    def test_network_edit_buttons_are_visible(self):
        expected_buttons = ["btn-undo", "btn-redo", "saveButton",
                            "cancelButton"]

        page = self.client.get(
            reverse('networkgen:detail', kwargs={"pk": self.network.pk})
        ).content.decode()

        for btn in expected_buttons:
            assert btn in page

    def test_the_basic_network_generator_formview(self):
        page = self.client.get(reverse('networkgen:create'))
        token = page.context.get("csrf_token")

        expected = TemplateResponse(
            page.wsgi_request,
            "networkgen/generator_form.html",
            context={"form": GeneratorForm(), "csrf_token": token}).render()

        assert page.content.decode() == expected.content.decode()


class NetworkGeneratorViews(TestCase):
    """Test the submission of SDF files to build networks."""

    def setUp(self):
        self.c = Client()

    def test_single_file_submission(self):
        assert Ligand.objects.count() == 0
        with open(Path(__file__).parent
                  / "test_files" / "CDK2_ligands.sdf", "rb") as f:
            sdf = BytesIO(f.read())

        response = self.c.post(reverse("networkgen:create"),
                               {"in_sdf": sdf, "metric": g.MCS},
                               follow=True)
        assert response.status_code == 200

        obj = response.context.get("object")

        # The 16 Ligands has been created
        assert obj.ligand_set.count() == 16

        # The page rendered contains the Json for the network. This is a proxy
        #  of the network view page.
        assert obj.network in response.content.decode()
