from django.template import Context
from django.test import TestCase

from cluster import forms


class ConfigForm(TestCase):
    def test_form_help_text(self):
        form = forms.ConfigForm()
        renderedForm = form.helper.render_layout(form, Context({}))

        print(renderedForm)

        assert "Cluster Config" in renderedForm
        assert 'name="forcefield_dir"' in renderedForm

    def test_form_buttons(self):
        form = forms.ConfigForm()
        renderedForm = form.helper.render_layout(form, Context({}))

        assert "Submit" in renderedForm
        assert "Reset" in renderedForm
