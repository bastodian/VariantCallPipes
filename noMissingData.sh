#!/bin/bash

### Script to parse vcf file and not allow for missing data.

INFILE=$1
OUTFILE=$2

grep -v -m 1 '##' $INFILE >> $OUTFILE

while read line
do
    if [[ "$line" == *"PASS"* ]] && [[ "$line" == *"FS=0.000"* ]] && [[ "$line" != *"./."* ]]
    then
        printf "%s\n" "$line" >> $OUTFILE
    fi
done < $INFILE
