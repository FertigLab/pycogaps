#! /bin/sh
set -e

python3 ./setup.py install --user
python3 ./vignette.py --user
