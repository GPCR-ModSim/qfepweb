from django.urls import path
from . import views as v
app_name = "networkgen"
urlpatterns = [
    path("", v.GeneratorView.as_view(), name="home"),
    path("detail/<uuid:pk>", v.GeneratorDetailView.as_view(), name="detail"),
    path("display/", v.NetworkGen.as_view(), name="index"),
    path('display/data/graph.json', v.NetworkData.as_view(), name='data')
]

