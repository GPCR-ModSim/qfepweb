from django.core.exceptions import ValidationError
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier


def valid_sdf(sdf_io):
    """Check that the file at "sdf_io" contains at least two molecules
    that RDKit is able to read."""
    suppl = list(ForwardSDMolSupplier(sdf_io))

    # Check molecule validity. RDKit returns a None molecule when reading fails.
    # see https://github.com/rdkit/rdkit/issues/642 if we want to extend this test.
    for i, mol in enumerate(suppl):
        if mol == None:
            raise ValidationError(
    			f"RDKit is unable to read molecule with index {i} in the input SDF file.")

    nligands = len(suppl)
    sdf_io.seek(0)
    # check that input SDF contains multiple ligands.
    if nligands == 0:
        raise ValidationError(
            "SDF file contains no ligands. Needs two at least.")

    if nligands == 1:
        raise ValidationError(
            "SDF file contains a single ligand. Needs two at least.")

