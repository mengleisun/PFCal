#!/bin/bash

for en in {3,5,7,10,30,50,70,100}
do
   mkdir /afs/cern.ch/work/m/msun/HGCAL_03_08/PFCal/PFCalEE/analysis/PLOTS/gamma_et${en}_eta1.6
   ./bin/egammaResolution -c scripts/DefaultConfig.cfg -n 1000 -i /afs/cern.ch/work/m/msun/HGCAL_03_08/PFCal/PFCalEE/V34/gitV00-03-08/gamma/  --digifilePath=/afs/cern.ch/work/m/msun/HGCAL_03_08/PFCal/PFCalEE/V34/gitV00-03-08/gamma/ -s HGcal_version34_model2_BOFF_et${en}_eta1.600  -r DigiIC2_version33_model2_BOFF_et${en}_eta1.600 -o ./PLOTS/gamma_et${en}_eta1.6.root --redoStep=2
done
