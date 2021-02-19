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
        page = self.client.get(reverse('home:index'))

        assert page.status_code == 200
        assert "QwebFEP" in page.content.decode()  # The title
