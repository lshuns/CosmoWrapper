#!/usr/bin/env bash

PIPEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
WDIR=$1
if [ "$WDIR" == "" ]
then
    echo "ERROR: ARG1 must be the path to the new working directory"
    exit 1
elif [[ -e $WDIR ]]
then
    echo "ERROR: target dir exsts"
    exit 1
elif [[ ! -d $WDIR ]]
then
    mkdir -pv $WDIR
fi

for script in ${PIPEDIR}/src/*{.R,.py,.sh}
do
    ln -sfv $script $WDIR/$(basename $script)
done
cp -v $PIPEDIR/src/vandenBusch2022.param $WDIR/
source ${WDIR}/vandenBusch2022.param

mkdir -pv ${WDIR}/${INPUTDIR}
mkdir -pv ${WDIR}/${OUTPUTDIR}
