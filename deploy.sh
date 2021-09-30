#!/bin/env sh
source ../bin/activate

RED='\033[0;31m'
NOCOLOR='\033[0m'

echo -e "${RED}Syncing with Git${NOCOLOR}"
git pull

echo -e "${RED} Installing any new requirements${NOCOLOR}"
pip install -r requirements.txt

echo -e "${RED} Syncing migrations${NOCOLOR}"
python manage.py migrate

echo -e "${RED} Updating static assets${NOCOLOR}"
python manage.py collectstatic --noinput -v 0

echo -e "${RED} Restarting the server${NOCOLOR}"
systemctl --user restart qfepweb
