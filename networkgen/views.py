from django.views.generic import TemplateView
from django.views.generic import CreateView
from networkgen.forms import GeneratorForm
from networkgen.models import Generator
from django.views.generic.detail import DetailView
from networkgen.mapgen import MapGen
import subprocess


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
    def form_valid(self, form):
        gen = MapGen(self.request.FILES["in_sdf"], form.data.get("metric"))
        gen.process_map()
        print(gen.as_json())
        return super().form_valid(form)
    def get_success_url(self):
        return self.request.path
        


class GeneratorDetailView(DetailView):
    """Create a new Database entry with all the parameters required for the FEP network generator.
    """
    model = Generator
