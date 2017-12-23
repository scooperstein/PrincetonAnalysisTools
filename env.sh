#!/usr/bin/env bash

readonly ANALYSIS_TOOLS_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


main() {
  export PYTHONPATH="$ANALYSIS_TOOLS_PATH/python:$PYTHONPATH"
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ANALYSIS_TOOLS_PATH:$ANALYSIS_TOOLS_PATH/HelperClasses:$ANALYSIS_TOOLS_PATH/plugins:$ANALYSIS_TOOLS_PATH/python"
}
main

