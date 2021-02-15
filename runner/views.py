from django.views.generic import CreateView

from runner.forms import RunnerForm
from runner.models import Runner


class Runner(CreateView):
    """The main entry.

    Create a new Database entry with all the parameters required for the run.
    """
    form_class = RunnerForm
    model = Runner
