"""The basic views for this project."""
from django.test import TestCase
from django.urls import reverse


class Home(TestCase):
    """Test the landing page.

    It should be available to everyone, without account, and a brief
    introdution to what the server does.

    """
    def setUp(self):
        pass

    def test_the_home_is_available(self):
        page = self.client.get(reverse('runner:index'))

        assert page.status_code == 200
        assert "QligFEP" in page.content.decode()  # The title

    def test_the_form_is_shown(self):
        """The basic input form is rendered in the home page."""
        page = self.client.get(reverse('runner:index')).content.decode()

        form_fields = ["mutation", "forcefield", "sampling", "windows",
                       "system", "temperatures", "replicates", "dual", "start"]

        for name in form_fields:
            assert f'div id="div_id_{name}"' in page
