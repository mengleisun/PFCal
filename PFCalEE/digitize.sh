#/bin/bash

for pu in {18..25}
do
./userlib/bin/digiSingle 10 root://eoscms//eos/cms/store/user/msun/V12/electron/gitV00-02-09/e-/HGcal_version12_model2_BOFF_et20_alpha0.361.root ~/HGCal_02_10/PFCal/PFCalEE/analysis/PLOTS/ 0-23:6,24-33:8 0-33:0.12 0-33:2 50 root://eoscms//eos/cms/store/user/msun/V12/MinBias/HGcal_version12_model2_BOFF_MinBias_${pu}.root > digi_${pu}.log
done
