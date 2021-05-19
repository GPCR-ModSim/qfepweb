from django.urls import path
from . import views as v
app_name = "networkgen"
urlpatterns = [
    path("", v.GeneratorView.as_view(), name='home'),
    path("detail/<uuid:pk>", v.GeneratorDetailView.as_view(), name="detail")
]