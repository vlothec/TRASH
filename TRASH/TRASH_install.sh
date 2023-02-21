#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

#check arguments, if --def present then pass to Rscript to use default R libpath

tar xfvz $SCRIPT_DIR/src/mafft-7.490-linux.tgz -C $SCRIPT_DIR/src

echo $SCRIPT_DIR
echo $1
VAR1="--def"

if [ "$1" = "$VAR1" ]
then
	echo "Installing TRASH with default library path"
	Rscript $SCRIPT_DIR/src/TRASH_install.R $VAR1
else
	echo "Installing TRASH with its own library path"
	Rscript $SCRIPT_DIR/src/TRASH_install.R
fi