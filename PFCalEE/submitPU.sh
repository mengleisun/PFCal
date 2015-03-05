#!/bin/bash

for file in {100..105}
do
  for i in {0..4}
  do
     j=$(($file-1))*5
     pu=`echo $i+$j+1|bc`
     ./submitProd.py -q 2nd -t V00-03-00 -v 12 -m 2 -n 1000 -f /afs/cern.ch/user/m/msun/HEPMC/Pythia140305_000${file}_${pu}.dat -o /afs/cern.ch/user/m/msun/HGCal_02_14/PFCal/PFCalEE/V12/MinBias/pile_${pu}/ -e /store/cmst3/group/hgcal/Standalone/V12/MinBias/pile_${pu}
  done
done

