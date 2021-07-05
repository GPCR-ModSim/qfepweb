from django.views.generic import TemplateView, CreateView
from django.urls import reverse
from networkgen.forms import GeneratorForm
from networkgen.models import Generator
from django.views.generic.detail import DetailView
from networkgen.mapgen import MapGen


class NetworkGen(TemplateView):
    template_name = "networkgen/index.html"


class NetworkData(TemplateView):
    content_type = "application/json; charset=utf-8"
    template_name = "networkgen/graph.json"


class GeneratorView(CreateView):
    """Create a new Database entry with all the parameters required for the FEP network generator.
    """
    form_class = GeneratorForm
    model = Generator

    def get_success_url(self):
        return reverse("networkgen:detail", kwargs={"pk": self.object.uuid})


class GeneratorDetailView(DetailView):
    """Create a new Database entry with all the parameters required for the FEP network generator.
    """
    model = Generator
