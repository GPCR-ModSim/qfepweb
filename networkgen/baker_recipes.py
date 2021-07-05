import io
from pathlib import Path
from django.core.files import File
from model_bakery.recipe import Recipe

from networkgen.models import Generator

def void():
    pass

sdf = File(
    open(Path(__file__).parent / "test_files" / "CDK2_ligands.sdf", "rb"),
    name="CDK2_ligands.sdf")

network = Recipe(
    Generator,
    in_sdf=sdf)
