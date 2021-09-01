from django.test import Client, TestCase
from django.urls import reverse


class Config(TestCase):
    def test_the_config_page_is_available(self):
        page = self.client.get(reverse('cluster:create'))

        assert page.status_code == 200
        assert "Cluster Config" in page.content.decode()
