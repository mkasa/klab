#!/bin/bash

# symlink is not allowed to install
rm -f bin/qtextasdata.py
# six==1.11.0 version requirement is too hard
if [[ "$($PYTHON -V | cut -d' ' -f2)" =~ "^3.8" ]]; then
    perl -ple 's/six==([\d\.]+)/six==1.14.0/' -i.orig setup.py
elif [[ "$($PYTHON -V | cut -d' ' -f2)" =~ "^3.9" ]]; then
    perl -ple 's/six==([\d\.]+)/six==1.15.0/' -i.orig setup.py
fi

# install as usual
$PYTHON setup.py install --single-version-externally-managed --record record.txt


