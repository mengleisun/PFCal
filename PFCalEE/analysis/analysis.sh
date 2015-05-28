#!/bin/bash

en=5
#for en in {5}
#do
  mkdir /afs/cern.ch/work/m/msun/HGCAL_03_08/PFCal/PFCalEE/digi/git_V00-03-03/version_30/model_2/gamma/BON/et_50/eta_2.000/a_2.000/run_0/eta2.0_et${en}_pu0_IC2
  ./bin/egammaResolution -c scripts/DefaultConfig.cfg -n 1000 --nRuns=0 -i root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V33//gitV00-03-07/gamma/ --digifilePath=root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V33//gitV00-03-07/gamma/ -s HGcal_version33_model2_BOFF_et${en}_eta1.600.root  -r DigiIC2_version33_model2_BOFF_et${en}_eta1.600.root  -o /afs/cern.ch/work/m/msun/HGCAL_03_08/PFCal/PFCalEE/digi/git_V00-03-03/version_30/model_2/gamma/BON/et_50/eta_2.000/a_2.000/run_0/eta2.0_et${en}_pu0_IC2.root --redoStep=2
#done

