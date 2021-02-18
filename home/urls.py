from . import views
from django.urls import path
from django.views.decorators.cache import cache_page

urlpatterns = [
    path('', views.index, name='index')
    ]
