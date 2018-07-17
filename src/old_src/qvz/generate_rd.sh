#!/bin/bash

idx=0
STATSFILE=`mktemp rd_statsXXXXXX`
while [ $idx -lt 20 ]; do
	# If doing ratiod encoding
	comp=$(echo "$idx*0.05" | bc -l)
	bin/qvz -c 1 -f $comp -s $1 $2 | tee -a $STATSFILE

	# Fixed rate encoding
	#comp=$(echo "$idx*0.10" | bc -l)
	#bin/qvz -c 3 -r $comp -s $1 $2 | tee -a $STATSFILE

	idx=$((idx+1))
done
awk '{print $2 $4 $6}' $STATSFILE > $3
rm -f $STATSFILE
