#! /bin/bash

2to3 --no-diffs -w -n landlab;
2to3 --no-diffs -d -w -n landlab;
2to3 --no-diffs -w -n ez_setup.py;
2to3 --no-diffs -w -n setup.py;
