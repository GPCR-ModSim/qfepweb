# Generated by Django 3.2.6 on 2021-08-10 12:13

from django.db import migrations, models
import networkgen.validators


class Migration(migrations.Migration):

    dependencies = [
        ('networkgen', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='generator',
            name='in_sdf',
            field=models.FileField(help_text='max. 20 Mbs', max_length=255, null=True, upload_to='', validators=[networkgen.validators.valid_sdf]),
        ),
        migrations.AlterField(
            model_name='ligand',
            name='image',
            field=models.ImageField(max_length=255, upload_to='molimages'),
        ),
    ]