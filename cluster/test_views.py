from datetime import time
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
        helper_text = forms.ConfigForm().fields["forcefield_directory"].help_text

        assert helper_text in page.content.decode()

    def test_labels_in_page(self):
        page = self.client.get(reverse('cluster:create'))
        labels = ["QDYN path", "QPREP path", "QFEP path"]

        for label in labels:
            assert label in page.content.decode()

    def test_placeholders_in_page(self):
        page = self.client.get(reverse('cluster:create'))

        for ph in ["/path/to/software/q/bin/qdyn", "DD HH:MM:SS"]:
            assert ph in page.content.decode()

    def test_form_submission(self):
        response = self.client.post(reverse('cluster:create'),
                                    {"forcefield_directory": "FF",
                                     "root_directory": "ROOT",
                                     "qdyn_path": "qdyn/path",
                                     "qprep_path": "qprep/path",
                                     "qfep_path": "qfep/path",
                                     "username": "userX",
                                     "runtime": time(23, 59, 59),
                                     "modules": "module1, module2, module3",
                                     "nodes": 1,
                                     "tasks": 16},
                                    follow=True)

        assert response.status_code == 200
        assert response.headers.get("Content-Type") == "text/plain; charset=utf-8"

        assert response.content.decode() == \
            """import os

# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")

# CLUSTER INPUTS. To add your own cluster, use the same input as below
CSB = {'NODES'        : '1',
       'NTASKS'       : '16'}
"""
