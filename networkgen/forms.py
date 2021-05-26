from crispy_forms.helper import FormHelper
from crispy_forms.layout import (Div, Fieldset, Layout)
from django.forms import ModelForm
from networkgen.models import Generator
from django import forms


class GeneratorForm(ModelForm):
    """A base form for the FEP network generator model."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_action = ''
        self.helper.layout = Layout(
            Fieldset(
                "FEP network generator parameters",
                Div(
                    Div("in_sdf", css_class="col-sm-12"),
                    Div("metric", css_class="col-sm-12"),
                    css_class="row")))

    class Meta:
        """Set the basic config for the FEP network form."""
        model = Generator
        fields = ["metric", "in_sdf"]
        help_texts = {
            "metric": "Distance metric to compare chemical similarity " +\
                "between ligand pairs.",
            "in_sdf": "Input .sdf file containing ligands."
        }
