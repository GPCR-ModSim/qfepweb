# Generated by Django 3.1.5 on 2021-02-16 12:39

from django.db import migrations, models
import uuid


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Runner',
            fields=[
                ('uuid', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('forcefield', models.CharField(choices=[('O15', 'OPLS 2015'), ('O05', 'OPLS 2005'), ('S_C', 'SIDECHAIN'), ('A14', 'AMBER14sb'), ('C36', 'CHARMM36'), ('C22', 'CHARMM Test')], default='O15', max_length=3)),
                ('sampling', models.CharField(choices=[('LIN', 'Linear'), ('SIG', 'Sigmoidal'), ('EXP', 'Exponential'), ('REX', 'Reverse exponential')], default='LIN', max_length=3)),
                ('cysbond', models.CharField(blank=True, max_length=16)),
                ('windows', models.PositiveIntegerField(default=1)),
                ('system', models.CharField(choices=[('PRO', 'Protein'), ('WAT', 'Water'), ('VAC', 'Vaccum')], default='PRO', max_length=3)),
                ('temperatures', models.CharField(blank=True, max_length=255)),
                ('replicates', models.PositiveIntegerField(default=1)),
                ('start', models.CharField(choices=[('E', 'Start in the endpoint'), ('H', 'Start in the middle')], default='E', max_length=3)),
                ('sphere_radius', models.PositiveIntegerField(default=15)),
            ],
        ),
    ]
