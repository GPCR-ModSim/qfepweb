import json
from os import name
from pathlib import Path
import uuid
from django.conf import settings
#from django.core.files import File
from django.db import models
from django.urls import reverse
from django_extensions.db.models import TimeStampedModel
from rdkit import Chem

from networkgen import mapgen


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

    def __str__(self):
        return f"Network Generator <{self.uuid}>"

    def _build_json(self, ligands, db_ligands):
        """Build a Json of the RDkit ligands paired with the DB data."""
        # For future ref:
        # edge_keys = ["label", "freenrg", "sem", "crashes", "from", "to"]
        # node_keys = ["shape", "label", "image", "id"]
        nodes = set([])
        result = {"nodes": [], "edges": []}
        key_ligands = {_.name: _ for _ in db_ligands}

        for ligand in ligands.values():
            nodes = nodes | set(
                [node for edge in ligand['Graph'].edges for node in edge])
            for edge in ligand['Graph'].edges:
                result["edges"].append(
                    {"from": str(key_ligands.get(edge[0]).uuid),
                     "to": str(key_ligands.get(edge[1]).uuid)})

        for node in nodes:
            node_ligand = key_ligands.get(node)
            result["nodes"].append(
                {"label": node_ligand.name,
                 "image": node_ligand.image.path,
                 "id": str(node_ligand.uuid)})

        return json.dumps(result)

    def build_network(self):
        """Calculate the network at this point. This includes

            - Do all the calculations for nodes and endges.
            - Build the images for the Ligands
            - Save those Ligands.
        """
        m = mapgen.MapGen(in_sdf=self.in_sdf, network_obj=self)
        m.make_map()

        if len(m.molecules) < 2:
            # FIXME: How does it work throwing an exception while serving?
            raise Exception(f"Number of ligands ({len(lignames)}) must be > 1")

        img_dir = "molimages"
        Path(settings.MEDIA_ROOT / "molimages").mkdir(exist_ok=True)

        ligands = []
        for molecule in m.molecules:
            moleculeImage = mapgen.MoleculeImage(molecule=molecule, core=m.mcs)
            ligand = Ligand(
                charge=Chem.rdmolops.GetFormalCharge(molecule),
                atom_number=len(molecule.GetAtoms()),
                name=moleculeImage.name,
                smiles=Chem.MolToSmiles(molecule, isomericSmiles=True),
                network=self)
            ligand.image = str(Path(img_dir) / f"{ligand.uuid}.png")
            with open(ligand.image.path, "wb") as png:
                png.write(moleculeImage.png())
            ligands.append(ligand)

        db_ligands = Ligand.objects.bulk_create(ligands)

        self.network = self._build_json(m.ligands, db_ligands)

    def save(self, *args, **kwargs):
        if not self.network:
            # Do this only the first time is saved
            self.build_network()
        super().save(*args, **kwargs)

    def get_absolute_url(self):
        """ Get absolute url. """
        return reverse("networkgen:detail", kwargs={"pk": self.pk})


class Ligand(models.Model):
    uuid = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    charge = models.IntegerField()
    atom_number = models.IntegerField()
    name = models.CharField(max_length=255)
    smiles = models.CharField(max_length=255)
    image = models.ImageField(max_length=255, upload_to="molimages")
    network = models.ForeignKey("Generator", on_delete=models.CASCADE)
