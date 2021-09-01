from django.urls import path
from . import views as v


app_name = "cluster"
urlpatterns = [
    path("", v.ConfigView.as_view(), name="create"),
]
