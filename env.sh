#!/usr/bin/env bash

readonly ANALYSIS_TOOLS_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

main() {
  export PYTHONPATH="$ANALYSIS_TOOLS_PATH/python:$PYTHONPATH"
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ANALYSIS_TOOLS_PATH:$ANALYSIS_TOOLS_PATH/HelperClasses:$ANALYSIS_TOOLS_PATH/plugins:$ANALYSIS_TOOLS_PATH/python"
}
main

python -c "import varial_ext" &> /dev/null
varial_nonexisting=$?

if [ $varial_nonexisting != 0 ]; then
    if [ -f Varial/setup.py ]; then
        echo "Updating Varial."
        cd Varial
        git pull
        cd -
    else
        echo "Installing Varial."
        git clone https://github.com/HeinerTholen/Varial
    fi

    export PYTHONPATH="$ANALYSIS_TOOLS_PATH/Varial:$PYTHONPATH"
fi
