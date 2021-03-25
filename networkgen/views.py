from django.views.generic import TemplateView, View


class NetworkGen(TemplateView):
    template_name = "networkgen/index.html"


class NetworkData(TemplateView):
    content_type = "application/json; charset=utf-8"
    template_name = "networkgen/graph.json"
