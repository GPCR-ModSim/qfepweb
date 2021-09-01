from crispy_forms.helper import FormHelper
from crispy_forms.layout import (ButtonHolder, Div, Layout, Fieldset, Reset,
                                 Submit)
from django import forms


class ConfigForm(forms.Form):

    forcefield_dir = forms.CharField(
        help_text="The directory to the input Force Field")
    nodes = forms.IntegerField(
        initial=1,
        min_value=1,
        max_value=32,
        help_text="The number of nodes in your Slurm cluster")
    tasks = forms.IntegerField(
        initial=8,
        min_value=1,
        max_value=32,
        help_text="The number of tasks per Slurm script.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.form_action = ''

        self.helper.layout = Layout(
            Fieldset(
                "Cluster Config",
                Div(
                    Div("forcefield_dir", css_class="col-sm-12"),
                    css_class="row"),
                Div(
                    Div("nodes", css_class="col-sm-4"),
                    Div("tasks", css_class="col-sm-8"),
                    css_class="row")),
            ButtonHolder(
                Submit("submit", "Submit", css_class="bg-cp1"),
                Reset("reset", "Reset")))
