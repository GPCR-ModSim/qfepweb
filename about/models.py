from django.db import models
from django.urls import reverse


class Person(models.Model):
    """A model to hold data about people on the team."""
    name = models.CharField(max_length=64)
    position = models.CharField(max_length=64)
    description = models.TextField(max_length=1024)
    picture = models.ImageField()
    url = models.URLField()

    def get_absolute_url(self):
        """Return the absolute url for the instance."""
        kwargs = {"pk": self.id, "name": self.name}

        return reverse("about:person_detail", kwargs=kwargs)

    def __str__(self):
        return f"{self.name} - {self.position}"
