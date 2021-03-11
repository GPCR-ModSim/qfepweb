from django.core.exceptions import ValidationError
from django.test import TestCase
from model_bakery import baker

from . import validators


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


class RunnerValidators(TestCase):
    def test_cysbond_validator(self):
        #"0", "1", "1:99,", "A", "1:99,35", "1:99,35:150,"
        assert validators.validate_cysbond("1:35") is None

        invalid_values = ["A:23", "1:99,35", "1:99,35:", "0:2", "1",
                          "1:99,35:150,", "1:99,"]
        for values in invalid_values:
            with self.assertRaises(ValidationError):
                validators.validate_cysbond(values)

    def test_temperatures_validation(self):
        assert validators.validate_temps("298,300") is None

        invalid_values = ["A", "-1", "298:300", "-1,298"]

        for value in invalid_values:
            with self.assertRaises(ValidationError):
                validators.validate_temps(value)
