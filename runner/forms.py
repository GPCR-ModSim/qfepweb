from crispy_forms.helper import FormHelper
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
