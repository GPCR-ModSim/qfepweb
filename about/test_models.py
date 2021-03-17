from pathlib import Path
from django.urls import reverse
from django.test import TestCase
from model_bakery import baker


class RunnerModel(TestCase):
    def setUp(self):
        self.person = baker.make("Person", _create_files=True)

    def tearDown(self):
        to_remove = Path(self.person.picture.path)
        if to_remove.is_file():
            to_remove.unlink()

    def test_model_fields(self):
        assert self.person.name
        assert self.person.position
        assert self.person.description
        assert self.person.url
        assert self.person.picture

    def test_model_friendly_string(self):
        assert str(self.person) == f"{self.person.name} - {self.person.position}"

    def test_model_absolute_url(self):
        assert self.person.get_absolute_url() == \
            reverse("about:person_detail",
                    kwargs={"pk": self.person.pk, "name": self.person.name})
