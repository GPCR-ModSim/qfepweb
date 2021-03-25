from django.urls import path
from . import views as v


app_name = "networkgen"

urlpatterns = [
    path('', v.NetworkGen.as_view(), name='index'),
    path('', v.NetworkData.as_view(), name='data')
]
