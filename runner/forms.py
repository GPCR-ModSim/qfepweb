from crispy_forms.helper import FormHelper
from crispy_forms.layout import (ButtonHolder, Div, Fieldset, Layout, Submit,
                                 Reset)
from django.forms import ModelForm

from runner.models import Runner
from runner import validators


class RunnerForm(ModelForm):
    """A base form for the input view, based on the model Runner."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.form_action = ''

        self.fields["cysbond"].validators.append(validators.validate_cysbond)
        self.fields["temperatures"].validators.append(validators.validate_temps)

        self.fields["cysbond"].widget.attrs["placeholder"] = "e.g. 1:99,35:150"
        self.fields["temperatures"].widget.attrs["placeholder"] = "e.g. 298,300"

        self.helper.layout = Layout(
            Fieldset(
                "Q Ligand FEP run parameters",
                Div(
                    Div("forcefield", css_class="col-sm-6"),
                    Div("sampling", css_class="col-sm-6"),
                    css_class="row"),
                Div(
                   Div("cysbond",  css_class="col-sm-6"),
                   Div("windows", css_class="col-sm-6"),
                   css_class="row"),
                Div(
                   Div("system", css_class="col-sm-6"),
                   Div("temperatures", css_class="col-sm-6"),
                   css_class="row"),
                Div(
                   Div("replicates", css_class="col-sm-4"),
                   Div("start", css_class="col-sm-4"),
                   Div("sphere_radius", css_class="col-sm-4"),
                   css_class="row")),
            ButtonHolder(
                Submit("submit", "Submit", css_class="bg-cp1"),
                Reset("reset", "Reset")))

    class Meta:
        """Set the basic config for the form."""
        model = Runner
        fields = ["forcefield", "sampling", "cysbond", "windows", "system",
                  "temperatures", "replicates", "start", "sphere_radius"]

        help_texts = {
            "forcefield": "The Forcefield to be used",
            "system": "What type of system we are setting up",
            "sphere_radius": "Size of the simulation sphere, in UNITS",  # FIXME
            "cysbond": "Add CYS-bonds as at1:at2,at3:at4",
            "start": "Starting FEP in the middle or endpoint",
            "temperatures": "Temperature(s), separated by commas",
            "replicates": "The number of repeats to run",
            "sampling": "Lambda spacing type to be used",
            "windows": "Total number of windows to run",
        }
