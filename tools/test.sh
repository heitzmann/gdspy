#!/bin/bash
VERSION=$1
if [ $(which virtualenv$VERSION) ]; then
  ENV=~/gdspy-env-$VERSION

  [ ! -d "$ENV" ] && virtualenv$VERSION --system-site-packages "$ENV"
  [ -d "$ENV/lib/python$VERSION/site-packages/gdspy" ] && rm -r "$ENV/lib/python$VERSION/site-packages/gdspy"

  source "$ENV/bin/activate"

  python setup.py build && python setup.py install

  echo
  echo Testing $VERSION
  cd docs/_static
  for i in *.py; do
    python "$i"
  done
  cd -

  deactivate
else
  echo "Usage: ./test.sh version"
fi
