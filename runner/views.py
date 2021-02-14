from django.views.generic import CreateView

from runner.models import Runner


class Runner(CreateView):
    """The main entry.

    Create a new Database entry with all the parameters required for the run.
    """
    model = Runner

    # TODO: Add here the fields needed in the view, or better yet move it to
    #  a dedicated Form
    fields = ["mutation"]
