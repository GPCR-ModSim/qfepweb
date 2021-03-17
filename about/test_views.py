"""The basic views for the about page(s)."""
from django.test import TestCase
from django.urls import reverse


class Home(TestCase):
    """Test the about page"""

    def setUp(self):
        pass

    def test_about_is_available(self):
        page = self.client.get(reverse("about:index"))
        assert page.status_code == 200

        assert "<title>About</title>" in page.content.decode()

    def test_about_is_linked_from_navigation(self):
        page = self.client.get(reverse("home:index"))

        assert reverse("about:index") in page.content.decode()

        # Assert one link in the title, another in the footer
        self.assertContains(page, reverse("about:index"), count=2)
