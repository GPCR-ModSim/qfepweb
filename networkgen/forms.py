from crispy_forms.helper import FormHelper
from crispy_forms.layout import (ButtonHolder, Div, Fieldset, Layout, Submit,
                                 Reset)
from django import forms
from django.forms import ModelForm
from networkgen.models import Generator
from networkgen import validators as v


class GeneratorForm(ModelForm):
    """A base form for the FEP network generator model."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # TODO: The following line allows the uploading of multiple files,
        #       but some work in the backend is needed before it's possible
        #self.fields["in_sdf"].widget.attrs.update({"multiple": True})
        self.helper = FormHelper()
        self.helper.form_action = ''
        self.helper.layout = Layout(
            Fieldset(
                "Create a network from a Sdf file",
                Div(
                    Div("in_sdf", css_class="col-sm-12"),
                    Div("in_pdb", css_class="col-sm-12"),
                    Div("metric", css_class="col-sm-12"),
                    css_class="row")),
            ButtonHolder(
                Submit("submit", "Submit", css_class="bg-cp1"),
                Reset("reset", "Reset")))

    class Meta:
        """Set the basic config for the FEP network form."""
        model = Generator
        fields = ["in_pdb", "in_sdf", "metric"]
        help_texts = {
            "metric": "Distance metric to compare chemical similarity " +\
                "between ligand pairs.",
            "in_sdf": "Input .sdf file containing ligands.",
            "in_pdb": "Input [optional] .pdb file containing proteins."
        }
