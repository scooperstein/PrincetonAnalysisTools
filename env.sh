#/bin/bash

PYTHONPATH=$PWD/python:$PYTHONPATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$PWD:$PWD/HelperClasses:$PWD/plugins:$PWD/python
