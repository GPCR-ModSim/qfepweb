import uuid
from django.db import models

from . import validators


class Runner(models.Model):
    """A model to hold each launch/run parameters."""
    O15 = "O15"
    O05 = "O05"
    S_C = "S_C"
    A14 = "A14"
    C36 = "C36"
    C22 = "C22"
    C_T = "C_T"
    FORCEFIELD_CHOICES = (
        (O15, "OPLS 2015"),
        (O05, "OPLS 2005"),
        (S_C, "SIDECHAIN"),
        (A14, "AMBER14sb"),
        (C36, "CHARMM36"),
        (C22, "CHARMM Test"))

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
        (HAL, "Start in the middle"))

    uuid = models.UUIDField(
        primary_key=True, default=uuid.uuid4, editable=False)
    forcefield = models.CharField(
        max_length=3, choices=FORCEFIELD_CHOICES, default=O15)
    sampling = models.CharField(
        max_length=3, choices=SAMPLING_CHOICES, default=LIN)
    # Cysbond is something like 20:150,45:187
    cysbond = models.CharField(max_length=128, blank=True,
                               validators=[validators.validate_cysbond])
    windows = models.PositiveIntegerField(default=1)
    system = models.CharField(
        max_length=3, choices=SYSTEM_CHOICES, default=PRO)
    temperatures = models.CharField(max_length=255, blank=True,
                                    validators=[validators.validate_temps])
    replicates = models.PositiveIntegerField(default=1)
    start = models.CharField(
        max_length=3, choices=START_CHOICES, default=ONE)
    sphere_radius = models.PositiveIntegerField(default=15)

    def __str__(self):
        """Return a proper string for the model."""
        return f'Run {self.uuid}'
