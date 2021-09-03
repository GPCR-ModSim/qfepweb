from crispy_forms.helper import FormHelper
from crispy_forms.layout import (ButtonHolder, Div, Layout, Fieldset, Reset,
                                 Submit)
from django import forms

# 'QDYN'       : 'qdyn=/home/x_aledi/software/q6/bin/qdynp', #fix qdyn= !!!!!
# 'QPREP'      : '/home/diazalej/software/q6/bin/qprep', # NOTE: change to where you are setting up, not where you are running!
# 'QFEP'       : '/home/x_aledi/software/q6/bin/qfep',
class ConfigForm(forms.Form):

    username = forms.CharField(
        help_text="The SLURM username used for submitting jobs on the cluster.")
    QDYN_path = forms.CharField(
        help_text="The path to the QDYN binary. For example: /path/to/software/q/bin/qdyn")
    QPREP_path = forms.CharField(
        help_text="The path to the QPREP binary. For example: /path/to/software/q/bin/qprep")
    QFEP_path = forms.CharField(
        help_text="The path to the QFEP binary. For example: /path/to/software/q/bin/qfep")
    root_directory = forms.CharField(
        help_text="The directory to run FEP in.")
    forcefield_directory = forms.CharField(
        help_text="The directory that contains the required Force Field. For example: /path/to/forcefields/")
    nodes = forms.IntegerField(
        initial=1,
        min_value=1,
        max_value=32,
        help_text="The number of nodes required per job.")
    tasks = forms.IntegerField(
        initial=8,
        min_value=1,
        max_value=32,
        help_text="The number of tasks per job.")
    runtime = forms.IntegerField(
        initial=4,
        min_value=1,
        max_value=32,
        help_text="The maximum job runtime. Input in d-hh:mm:ss.")
    modules = forms.CharField(
        help_text="Additional SLURM modules required per job. Input as comma-separated list.",
        required=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.form_action = ''

        self.helper.layout = Layout(
            Fieldset(
                "Cluster Configuration file generation",
                Div(
                    Div("username", css_class="col-sm-4"),
                    Div("root_directory", css_class="col-sm-8"),
                    css_class="row"),
                Div(
                    Div("QDYN_path", css_class="col-sm-12"),
                    css_class="row"),
                Div(
                    Div("QPREP_path", css_class="col-sm-12"),
                    css_class="row"),
                Div(
                    Div("QFEP_path", css_class="col-sm-12"),
                    css_class="row"),
                Div(
                    Div("forcefield_directory", css_class="col-sm-12"),
                    css_class="row"),
                Div(
                    Div("nodes", css_class="col-sm-4"),
                    Div("tasks", css_class="col-sm-4"),
                    css_class="row")),
                Div(
                    Div("runtime", css_class="col-sm-4"),
                    Div("modules", css_class="col-sm-4"),
                    css_class="row"),
            ButtonHolder(
                Submit("submit", "Download configuration file", css_class="bg-cp1"),
                Reset("reset", "Reset")))
