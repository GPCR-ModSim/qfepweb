from django.test import TestCase
from model_bakery import baker


class RunnerModel(TestCase):
    def setUp(self):
        self.runner = baker.make("Runner")

    def test_model_has_a_string(self):
        assert str(self.runner) == f'Run {self.runner.uuid}'

    def test_models_fields(self):
        assert self.runner.forcefield
        assert self.runner.sampling
        assert self.runner.cysbond == ""
        assert self.runner.windows
        assert self.runner.system
        assert self.runner.temperatures == ""
        assert self.runner.replicates
        assert self.runner.start
        assert self.runner.sphere_radius
