from django.urls import path

from . import views as v

app_name = "utils"

urlpatterns = [
    path('', v.About.as_view(), name='about'),
    path('', v.Help.as_view(), name='help')
]
