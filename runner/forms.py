from crispy_forms.helper import FormHelper
from django.forms import ModelForm

from runner.models import Runner


class RunnerForm(ModelForm):
    """A base form for the input view, based on the model Runner."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.helper = FormHelper()
        self.helper.form_action = ''

    class Meta:
        """Set the basic config for the form."""
        model = Runner
        fields = ["mutation", "forcefield", "sampling", "windows", "system",
                  "temperatures", "replicates", "dual", "start"]

        help_texts = {
            "mutation": "The desired mutation given as Wild-ResN-Mut, "
                "e.g. Y197A. ResN is taken from original .pdb"}

