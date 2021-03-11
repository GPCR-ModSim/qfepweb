from django.core.exceptions import ValidationError
from django.test import TestCase
from runner import forms


class RunnerFormT(TestCase):
    """Test the Runner Input Form."""
    def setUp(self):
        self.form_data = {"forcefield": "O15",
                          "sampling": "LIN",
                          "windows": 10,
                          "system": "WAT",
                          "replicates": 1,
                          "start": "E",
                          "sphere_radius": 15}

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

    def test_minimal_valid_form(self):
        form = forms.RunnerForm(data=self.form_data)

        assert form.is_valid(), form.errors

    def test_cys_validation(self):
        self.form_data["cysbond"] = "1:99,35:150"

        form = forms.RunnerForm(data=self.form_data)

        assert form.is_valid(), form.errors

        invalid_inputs = ["0", "-1", "1", "1:99,", "A", "1:99,35",
                          "1:99,35:150,", "0:39", "35:15", "15,15"]

        for cys_value in invalid_inputs:
            self.form_data["cysbond"] = cys_value
            form = forms.RunnerForm(data=self.form_data)
            assert not form.is_valid(), "Should be invalid: {cys_value}"

    def test_temperatures_validation(self):
        self.form_data["temperatures"] = "298,300"
        form = forms.RunnerForm(data=self.form_data)

        assert form.is_valid(), form.errors

        invalid_inputs = ["A", "-1", "298:300", "-1,298"]

        for temp_value in invalid_inputs:
            self.form_data["temperatures"] = temp_value
            form = forms.RunnerForm(data=self.form_data)
            assert not form.is_valid(), f"Should be invalid: {temp_value}"
