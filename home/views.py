# from django.views.decorators.cache import cache_page
# from django.utils.decorators import method_decorator
from django.views.generic import TemplateView


# @method_decorator(cache_page(60 * 60), name='dispatch')
class Home(TemplateView):

    template_name = "home/index.html"
