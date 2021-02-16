# qfepweb
Development of the web interface to QligFEP.

# Install for development

    $ python -m venv qfepweb
    $ cd qfepweb
    $ source bin/activate
    $ python -m pip install --upgrade pip
    $ git clone git@github.com:GPCR-ModSim/qfepweb.git

At `bin` there is a file called `activate`. Create a new file called
`postactivate` with the following contents:

    export PYTHONPATH=$VIRTUAL_ENV/qfepweb
    export DJANGO_SETTINGS_MODULE=qfepweb.settings

Load it with `source bin/postactivate`

    $ cd qfepweb
    $ pip install -r requirements_devel.txt

Create a file called `.env` and add the following lines. **Never commit this
file**

    DEBUG=on  # or off, for production
    SECRET_KEY=somegoodrandomstring_atleast50length
    # e.g. )Ql?ddAC3jMTybSHe^(vk@Sd=nTF#q&S5m:f$:U8O27^cr5S

If everything went smoothly, we could check Django testing page with:

    $ django-admin runserver

# First migration

You should have seen a warning like:

    You have 18 unapplied migration(s)

The current database is the default `db.sqlite3` (ignored by .gitignore). The
pending migrations are the Django default tables. Apply it with:

    $ python manage.py migrate