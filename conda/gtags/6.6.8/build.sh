#!/bin/bash

set -e
set -o pipefail

./configure --prefix=${PREFIX}
make
make install
