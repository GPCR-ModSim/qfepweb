"""The basic views for this project."""
import json
from django.template.response import TemplateResponse
from django.test import TestCase
from django.urls import reverse
from networkgen.forms import GeneratorForm


class NetworkGen(TestCase):
    """Test the landing of the generator pages.

    It should be available to everyone, without account, and a brief
    introdution to what the server does.

    """
    def setUp(self):
        pass

    def test_the_home_is_available(self):
        page = self.client.get(reverse('networkgen:index'))

        assert page.status_code == 200
        assert "FEP Network" in page.content.decode()  # The title

    def test_the_javascript_and_css_is_included(self):
        page = self.client.get(reverse('networkgen:index')).content.decode()

        # There are proxies to identify the static assets needed
        requiredJs = ["lodash.min.js", "vis.js", "networkgen.js"]
        requiredCss = ["vis-network.min.css", "use.fontawesome.com",
                       "networkgen.css"]

        for asset in requiredJs + requiredCss:
         assert asset in page

    def test_the_serving_of_Json_data(self):
        """Right now the json data is a plain file. In the future it might be
        a file generated from some dynamic data."""

        resp = self.client.get(reverse('networkgen:data'))
        page = resp.content.decode()

        assert resp["Content-Type"] == "application/json; charset=utf-8"

        jsonData = json.loads(page)

        assert list(jsonData.keys()) == ['nodes', 'edges']

    def test_network_edit_buttons_are_visible(self):
        expected_buttons = ["btn-undo", "btn-redo", "saveButton",
                            "cancelButton"]

        page = self.client.get(reverse('networkgen:index')).content.decode()

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
