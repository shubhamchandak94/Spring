#!/bin/bash

rm config.py

for f in config*.py
do
	cp $f config.py
	stdbuf -oL echo $f |& tee -a logs/"$f".log
	stdbuf -oL cat $f |& tee -a logs/"$f".log
	stdbuf -oL /usr/bin/time -v ./run.sh -pgc |& tee -a logs/"$f"_thresh4.log
#	stdbuf -oL /usr/bin/time -v ./run.sh -d |& tee -a logs/"$f".log
#	stdbuf -oL ./run.sh -o |& tee -a logs/"$f".log
	rm config.py
done
