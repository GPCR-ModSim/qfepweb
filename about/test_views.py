"""The basic views for the about page(s)."""
from pathlib import Path
from django.test import TestCase
from django.urls import reverse
from model_bakery import baker

from .models import Person


class Home(TestCase):
    """Test the about page"""

    def setUp(self):
        self.people = []

        self.people.append(baker.make("Person", _create_files=True))

    def tearDown(self):
        for p in self.people:
            to_remove = Path(p.picture.path)
            print(to_remove)
            if to_remove.is_file():
                to_remove.unlink()

    def test_about_is_available(self):
        page = self.client.get(reverse("about:index"))
        assert page.status_code == 200

        assert "<title>About</title>" in page.content.decode()

    def test_about_is_linked_from_navigation(self):
        page = self.client.get(reverse("home:index"))

        assert reverse("about:index") in page.content.decode()

        # Assert one link in the title, another in the footer
        self.assertContains(page, reverse("about:index"), count=2)

    def test_about_includes_skeleton_page(self):
        with self.assertTemplateUsed("upper_navbar.html"):
            self.client.get(reverse("about:index"))

        with self.assertTemplateUsed("footer.html"):
            self.client.get(reverse("about:index"))

    def test_persons_in_the_project_appears(self):
        page = self.client.get(reverse("about:index"))

        content = page.content.decode()

        assert self.people[0].name in content
        assert self.people[0].get_absolute_url() in content
        assert self.people[0].picture.url in content
        assert self.people[0].position in content
