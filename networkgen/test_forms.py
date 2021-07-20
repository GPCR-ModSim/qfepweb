from pathlib import Path
from django.core.exceptions import ValidationError
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from networkgen import forms
from networkgen.models import Generator as g
from networkgen import validators as v


class GeneratorForm(TestCase):
    """Test the Generator Input Form"""
    def setUp(self):
        in_sdf = Path(__file__).parent / "test_files" / "Two_Ligands.sdf"
        with open(in_sdf, "rb") as inFile:
            content = inFile.read()

        self.multipart = {"in_sdf": SimpleUploadedFile(in_sdf, content)}
        self.form_data = {"metric": g.SMILES}

    def test_form_is_a_model_form(self):
        form = forms.GeneratorForm()

        for field in self.form_data.keys():
            assert field in form.fields.keys()

    def test_form_help_texts(self):
        form = forms.GeneratorForm()
        renderedForm = form.as_p()

        assert "Distance metric to compare chemical similarity between " +\
            "ligand pairs." in renderedForm
        assert "Input .sdf file containing ligands." in renderedForm

    def test_form_rendering(self):
        form = forms.GeneratorForm()
        renderedForm = form.as_p()

        assert '<label for="id_in_sdf">In sdf:</label>' in renderedForm
        assert '<input type="file" name="in_sdf" required id="id_in_sdf">' in \
            renderedForm

    def test_minimal_valid_form(self):
        form = forms.GeneratorForm(data=self.form_data,
                                   files=self.multipart)

        assert form.is_valid(), form.errors

    def test_one_file_validator(self):
        in_sdf = Path(__file__).parent / "test_files" / "One_ligand.sdf"
        with open(in_sdf, "rb") as inFile:
            content = inFile.read()
        sdf = SimpleUploadedFile(in_sdf, content)

        with self.assertRaisesMessage(
            ValidationError, "SDF file contains a single ligand. Needs two at least."):
            v.valid_sdf(sdf)

    def test_some_invalid_sdfs(self):
        in_sdf = Path(__file__).parent / "test_files" / "One_ligand.sdf"
        with open(in_sdf, "rb") as inFile:
            content = inFile.read()

        self.multipart = {"in_sdf": SimpleUploadedFile(in_sdf, content)}
        self.form_data = {"metric": g.SMILES}

        form = forms.GeneratorForm(data=self.form_data,
                                   files=self.multipart)
        assert not form.is_valid()

        assert form.errors.as_data()["in_sdf"][0].message == \
            "SDF file contains a single ligand. Needs two at least."

    def test_empty_file_validator(self):
        in_sdf = Path(__file__).parent / "test_files" / "Empty_ligand.sdf"
        with open(in_sdf, "rb") as inFile:
            content = inFile.read()
        sdf = SimpleUploadedFile(in_sdf, content)

        with self.assertRaisesMessage(
            ValidationError, "SDF file contains no ligands. Needs two at least."):
            v.valid_sdf(sdf)

    def test_invalid_sdf(self):
        in_sdf = Path(__file__).parent / "test_files" / "Invalid.sdf"
        with open(in_sdf, "rb") as inFile:
            content = inFile.read()
        sdf = SimpleUploadedFile(in_sdf, content)

        with self.assertRaisesMessage(
            ValidationError,
            "RDKit is unable to read molecule with index 0 in the input SDF file."):
            v.valid_sdf(sdf)
