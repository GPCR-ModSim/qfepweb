# qfepweb
Development of the web interface to QligFEP.

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

If everything went smoothly, we could check Django testing page with:

  $ django-admin runserver

