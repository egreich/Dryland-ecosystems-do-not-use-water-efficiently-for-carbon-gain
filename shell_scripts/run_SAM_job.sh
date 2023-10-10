#!/bin/bash
# analysis.sh - an analysis program
# $1 (input) and $2 (output) are the first and second arguments to this script

# strip off the directory paths to get just the filename
# BASE=`basename $1`
site=$1
seed=$2
modelv=$3

echo 
echo "run('$1');"
date

# Run multicomp script
#R --no-save < $1
./scripts/02_script_SAM.R $site $seed $modelv
echo "run('$1'); done"
date
