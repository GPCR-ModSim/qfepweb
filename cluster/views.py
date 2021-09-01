from django.views.generic import TemplateView


class ConfigView(TemplateView):
    template_name = "cluster/config_create.html"
