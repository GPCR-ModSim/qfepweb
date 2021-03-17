from django.views.generic import DetailView, ListView
from about.models import Person


class Index(ListView):
    model = Person


class PersonDetails(DetailView):
    model = Person
