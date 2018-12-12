#!/usr/bin/env bash

#module load coltron
module load bamliquidator
module load R

if [ $# -ne 6 ] && [ $# -ne 9 ]
then
	echo "Usage: `basename $0` {coltron path} {enhancer file} {bam} {genome} {output folder} {job name} {distance} {exp file} {exp cutoff}"
	exit 65
fi

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export PATH=$PATH:$script_dir

if [ $# -eq 6 ]
then
	echo "$1 -e $2 -b $3 -g $4 -o $5 -n $6"
	$1 -e $2 -b $3 -g $4 -o $5 -n $6
fi

if [ $# -eq 9 ]
then
	echo "$1 -e $2 -b $3 -g $4 -o $5 -n $6 -d $7 -a $8 -x $9"
	$1 -e $2 -b $3 -g $4 -o $5 -n $6 -d $7 -a $8 -x $9
fi


