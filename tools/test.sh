#!/bin/sh

VERSION=$1
if [ "$(command -v "virtualenv$VERSION")" ]; then
  ENV=~/gdspy-env-$VERSION

  [ ! -d "$ENV" ] && "virtualenv$VERSION" --system-site-packages "$ENV"
  [ -d "$ENV/lib/python$VERSION/site-packages/gdspy" ] && rm -r "$ENV/lib/python$VERSION/site-packages/gdspy"

  . "$ENV/bin/activate"

  python setup.py build && python setup.py install

  echo
  echo "Testing $VERSION"
  cd docs/_static || exit
  for i in *.py; do
    python "$i"
  done
  cd - || exit

  deactivate
else
  echo "Usage: $0 version"
fi
