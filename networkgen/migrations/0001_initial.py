# Generated by Django 3.2.2 on 2021-05-26 17:30

from django.db import migrations, models
import django.db.models.deletion
import django_extensions.db.fields
import uuid


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Generator',
            fields=[
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('uuid', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('metric', models.CharField(choices=[('SMILES', 'SMILES'), ('MFP', 'Morgan fingerprints (MFP)'), ('Tanimoto', 'Tanimoto fingerprints (TFP)'), ('MCS', 'Maximum common substructure (MCS)')], default='MFP', max_length=155)),
                ('in_sdf', models.FileField(help_text='max. 20 Mbs', max_length=255, null=True, upload_to='')),
                ('network', models.JSONField(null=True)),
            ],
            options={
                'get_latest_by': 'modified',
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Ligand',
            fields=[
                ('uuid', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
                ('charge', models.IntegerField()),
                ('atom_number', models.IntegerField()),
                ('name', models.CharField(max_length=255)),
                ('smiles', models.CharField(max_length=255)),
                ('image', models.ImageField(max_length=255, upload_to='')),
                ('network', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='networkgen.generator')),
            ],
        ),
    ]
