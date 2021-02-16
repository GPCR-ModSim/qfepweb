from django.test import TestCase
from runner import forms


class RunnerFormT(TestCase):
    """Test the Runner Input Form."""
    def test_form_is_a_modelform(self):
        """The form autoloads fields from the model."""
        form = forms.RunnerForm()

        for field in ["forcefield", "sampling", "cysbond", "windows", "system",
                      "temperatures", "replicates", "start", "sphere_radius"]:
            assert field in form.fields.keys()

    def test_form_help_texts(self):
        form = forms.RunnerForm()

        assert "The Forcefield to be used" in \
            form.as_p()
        self.fail("Keep working here, adding help texts.")
