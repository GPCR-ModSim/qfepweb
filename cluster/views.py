from django.views.generic import FormView, TemplateView
from cluster import forms


class ConfigView(FormView):
    form_class = forms.ConfigForm
    template_name = "cluster/config_create.html"
