#!/bin/sh

ENV=~/gdspy-env

[ ! -d "$ENV" ] && virtualenv --system-site-packages "$ENV"
[ -d "$ENV/lib/python/site-packages/gdspy" ] && rm -r "$ENV/lib/python/site-packages/gdspy"

. "$ENV/bin/activate"

python setup.py build && python setup.py install

echo
cd docs/_static || exit
for i in *.py; do
  python "$i"
done
cd - || exit

deactivate
