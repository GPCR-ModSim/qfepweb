from django.urls import path
from . import views as v


app_name = "cluster"
urlpatterns = [
    path("", v.ConfigCreate.as_view(), name="create"),
    path("detail", v.ConfigDetail.as_view(), name="detail"),
]
