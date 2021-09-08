from datetime import timedelta
from django.template import Context
from django.test import TestCase

from cluster import forms


class ConfigForm(TestCase):
    def test_form_help_text(self):
        form = forms.ConfigForm()
        renderedForm = form.helper.render_layout(form, Context({}))

        print(renderedForm)

        assert "Cluster Config" in renderedForm
        assert 'name="forcefield_directory"' in renderedForm

    def test_form_buttons(self):
        form = forms.ConfigForm()
        renderedForm = form.helper.render_layout(form, Context({}))

        assert "Download configuration file" in renderedForm
        assert "Reset" in renderedForm

    def test_runtime_cleaning(self):
        data = {"forcefield_directory": "FF",
                "root_directory": "ROOT",
                "qdyn_path": "qdyn/path",
                "qprep_path": "qprep/path",
                "qfep_path": "qfep/path",
                "username": "userX",
                "runtime": timedelta(hours=23, minutes=59, seconds=59),
                "modules": "module1, module2, module3",
                "nodes": 1,
                "tasks": 16}

        form = forms.ConfigForm(data=data)
        if form.is_valid():
            assert form.clean_runtime() == "23:59:59"

        data["runtime"] = timedelta(days=2, hours=10)
        form = forms.ConfigForm(data=data)
        if form.is_valid():
            assert form.clean_runtime() == "2d-10:00:00"

        data["runtime"] = timedelta(hours=24)
        form = forms.ConfigForm(data=data)
        if form.is_valid():
            assert form.clean_runtime() == "1d-0:00:00"
