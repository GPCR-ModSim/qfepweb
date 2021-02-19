from django.shortcuts import render
from django.views.decorators.cache import cache_page


#@cache_page(60 * 60 * 24)
def index(request):

    context = {}
    return render(request, 'home/index.html', context)



