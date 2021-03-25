from django.views.generic import TemplateView, View


class NetworkGen(TemplateView):
    template_name = "networkgen/index.html"


class NetworkData(View):
    pass
