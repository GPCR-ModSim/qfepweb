from django.test import Client, TestCase
from django.urls import reverse

from cluster import forms


class Config(TestCase):
    def test_the_config_page_is_available(self):
        page = self.client.get(reverse('cluster:create'))

        assert page.status_code == 200
        assert "Cluster Config" in page.content.decode()

    def test_helper_texts_in_page(self):
        page = self.client.get(reverse('cluster:create'))
        helper_text = forms.ConfigForm().fields["forcefield_dir"].help_text

        assert helper_text in page.content.decode()
