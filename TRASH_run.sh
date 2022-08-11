#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

#check arguments, if --def present then pass to Rscript to use default R libpath
echo $SCRIPT_DIR
echo "TRASH will be run with $@ arguments"


Rscript $SCRIPT_DIR/src/TRASH_run.R $@
