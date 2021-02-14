import uuid
from django.db import models


class Runner(models.Model):
    """A model to hold each launch/run parameters."""
    O15 = "O15"
    O05 = "O05"
    S_C = "S_C"
    A14 = "A14"
    C36 = "C36"
    FORCEFIELD_CHOICES = (
        (O15, "OPLS 2015"),
        (O05, "OPLS 2005"),
        (S_C, "SIDECHAIN"),
        (A14, "AMBER14sb"),
        (C36, "CHARMM36"))

    LIN = "LIN"
    SIG = "SIG"
    EXP = "EXP"
    REX = "REX"
    SAMPLING_CHOICES = (
        (LIN, "Linear"),
        (SIG, "Sigmoidal"),
        (EXP, "Exponential"),
        (REX, "Reverse exponential"))

    PRO = "PRO"
    WAT = "WAT"
    VAC = "VAC"
    SYSTEM_CHOICES = (
        (PRO, "Protein"),
        (WAT, "Water"),
        (VAC, "Vaccum"))

    ONE = "E"
    HAL = "H"
    START_CHOICES = (
        (ONE, "Start in the endpoint"),
        (HAL, "Start in the middle, recommended for dual"))

    uuid = models.UUIDField(
        primary_key=True, default=uuid.uuid4, editable=False)
    mutation = models.CharField(max_length=10, blank=True)
    forcefield = models.CharField(
        max_length=3, choices=FORCEFIELD_CHOICES, default=O15)
    sampling = models.CharField(
        max_length=3, choices=SAMPLING_CHOICES, default=LIN)
    windows = models.PositiveIntegerField(default=1)
    system = models.CharField(
        max_length=3, choices=SYSTEM_CHOICES, default=PRO)
    # WARNING: the following field is a list field and should include
    # validators as such
    temperatures = models.CharField(max_length=255, blank=True)
    replicates = models.PositiveIntegerField(default=1)
    dual = models.BooleanField(default=False)
    start = models.CharField(
        max_length=3, choices=START_CHOICES, default=ONE)

    def __str__(self):
        """Return a proper string for the model."""
        return f'Run {self.uuid}'
