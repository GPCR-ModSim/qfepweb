from os import name
import uuid
from django.db import models
from django.urls import reverse
from django_extensions.db.models import TimeStampedModel


class Generator(TimeStampedModel):
    """A model that holds parameters for the FEP network generator."""

    SMILES = "SMILES"
    MFP = "MFP"
    Tanimoto = "Tanimoto"
    MCS = "MCS"
    METRICS_CHOICES = (
        (SMILES, "SMILES"),
        (MFP, "Morgan fingerprints (MFP)"),
        (Tanimoto, "Tanimoto fingerprints (TFP)"),
        (MCS, "Maximum common substructure (MCS)"))

    uuid = models.UUIDField(
        primary_key=True, default=uuid.uuid4, editable=False)
    metric = models.CharField(
        max_length=155, choices=METRICS_CHOICES, default=MFP)
    in_sdf = models.FileField(max_length=255, null=True, help_text='max. 20 Mbs')
    network = models.JSONField(null=True)

    def get_absolute_url(self):
        """ Get absolute url. """
        return reverse("networkgen:detail", kwargs={"pk": self.pk})


class Ligand(models.Model):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    charge = models.IntegerField()
    atom_number = models.IntegerField()
    name = models.CharField(max_length=255)
    smiles = models.CharField(max_length=255)
    image = models.ImageField(max_length=255)
    network = models.ForeignKey("Generator", on_delete=models.CASCADE)
