from django.core.exceptions import ValidationError
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier


def two_ligands(sdf_io):
    """Check that the file at "sdf_io" contains at least two molecules."""
    suppl = list(ForwardSDMolSupplier(sdf_io))

    if len(suppl) == 0:
        raise ValidationError(
            "SDF file contains no ligands. Needs two at least")

    if len(suppl) == 1:
        raise ValidationError(
            "SDF file contains a single ligand. Needs two at least")

    sdf_io.seek(0)
