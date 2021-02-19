from django.urls import path
from . import views as v

app_name = "runner"

urlpatterns = [
    path('', v.Runner.as_view(), name='qligfep')
]
