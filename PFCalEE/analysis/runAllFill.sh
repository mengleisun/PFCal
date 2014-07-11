#!/bin/sh

INPATH=root://eoscms//eos/cms/store/user/amagnan/HGCalHEGeant4/gitV00-02-01/

for v in 23
do
    for p in e- pi-
    do
	MYDIR=PLOTS/gitV00-02-03/version$v/$p/
	mkdir -p $MYDIR
	for r in 0
	do
	    #for e in 10 15 18 20 25 30 35 40 45 50 60 80
	    for e in 10 18 20 35 40 45 50 60 80
	    #for e in 50 60 80 100 150 200 300 400 500
	    #for e in 5 10 15 20 25 30 40 50 60 #80 100 150 200 300 400 500
	    do
		echo "Processing v${v}_p${p}_e${e}_run${r}"
		./bin/validation 0 $INPATH/$p/ HGcal_version${v}_model3_BOFF_e${e}.root Digi_version${v}_model3_BOFF_e${e}.root $MYDIR/validation_e${e}.root 4 >& v${v}_p${p}_e${e}.log &
	    done
	done
    done
done
