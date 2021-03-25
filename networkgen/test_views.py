"""The basic views for this project."""
from django.test import TestCase
from django.urls import reverse


class NetworkGen(TestCase):
    """Test the landing of the runner pages.

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
