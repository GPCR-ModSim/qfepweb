from django.core.exceptions import ValidationError
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier


def two_ligands(sdf_io):
    """Check that the file at "sdf_io" contains at least two molecules."""
    suppl = list(ForwardSDMolSupplier(sdf_io))
    if len(suppl) == 1:
        raise ValidationError(
            "SDF with only one ligand. Needs two at least")
    sdf_io.seek(0)
