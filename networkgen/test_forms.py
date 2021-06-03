from pathlib import Path
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from networkgen import forms
from networkgen.models import Generator as g


class GeneratorForm(TestCase):
    """Test the Generator Input Form"""
    def setUp(self):
        in_sdf = Path(__file__).parent / "test_files" / "One_ligand.sdf"
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

    def test_some_invalid_sdfs(self):
        self.fail("Continue testing here")
