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
        renderedForm = form.as_p()

        assert "The Forcefield to be used" in renderedForm
        assert "What type of system we are setting up" in renderedForm
        assert "Size of the simulation sphere, in UNITS" in renderedForm
        assert "Add CYS-bonds as at1:at2,at3:at4" in renderedForm
        assert "Starting FEP in the middle or endpoint" in renderedForm
        assert "Temperature(s), separated by commas" in renderedForm
        assert "The number of repeats to run" in renderedForm
        assert "Lambda spacing type to be used" in renderedForm
        assert "Total number of windows to run" in renderedForm

    def test_form_help_textholders(self):
        form = forms.RunnerForm()
        renderedForm = form.as_p()

        assert "e.g. 1:99,35:150" in renderedForm
        assert "e.g. 298,300" in renderedForm
