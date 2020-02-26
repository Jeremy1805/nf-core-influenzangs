#!/usr/bin/env bash

Rchosenlist=$1
tupleid=$2

listcount=`echo "$Rchosenlist" | wc -w`

if [ $listcount -gt 1 ]; then

	if [ -f "${tupleid}.aligninput" ]; then
		rm "${tupleid}.aligninput"
	fi

	for refseg in ${Rchosenlist}
	do
				cat ${refseg}.fasta >> "${tupleid}.aligninput"
	done

	clustalw -INFILE=${tupleid}.aligninput -ALIGN -OUTPUT=FASTA

fi
