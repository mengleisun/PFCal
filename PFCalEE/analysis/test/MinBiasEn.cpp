#define PI 3.14159265

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <cmath>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TF1.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

int main(int argc, char** argv){//main  


  std::ostringstream inFilePath("root://eoscms//eos/cms/store/user/msun/V12/gamma//gitV00-03-00/gamma/HGcal_version12_model2_BOFF_et60_eta1.900_run1.root");
  ///////////// input file ////////////////
  TFile *simFile = TFile::Open(inFilePath.str().c_str());

  TTree *lSimTree = (TTree*)simFile->Get("HGCSSTree");
  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  
  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  
  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
  //for HGCAL, true means only 12 FHCAL layers considered (24 are simulated)
  bool concept = true;

  bool selectEarlyDecays = false;

  double minX=-1700,maxX=1700;
  double minY=-1700,maxY=1700;
  double minZ=3170,maxZ=5070;

  unsigned nX=(maxX-minX)/10,nY=(maxY-minY)/10;
  unsigned nZ=maxZ-minZ;

  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  unsigned debug = 0;

  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  HGCSSInfo * info=(HGCSSInfo*)simFile->Get("Info");
  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  //models 0,1 or 3.
  bool isTBsetup = (model != 2);
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;

  //extract input energy

  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();

  myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

  //initialise calibration class
  HGCSSCalibration mycalib(inFilePath.str().c_str());

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;



  TFile *outputFile = TFile::Open("MinbiasEnergy.root","RECREATE");
  outputFile->cd();

  TH1F *p_pdgID = new TH1F("p_pdgID","pdgID",
			      2000,-1000,1000);
  TH1F *p_en = new TH1F("p_en","energy",
                              150,1.5,3.0);
  TProfile *p_longpro = new TProfile("p_longpro","longitudinal profile",
                              30,1,30);


  double radLength[30];

  for (unsigned ievt(0); ievt<1000; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

    lSimTree->GetEntry(ievt);

    if(ievt==0){
      for(unsigned iL(0); iL<30; iL++){
      double w = (*ssvec)[iL].volX0trans()/(*ssvec)[1].volX0trans();
      std::cout << "layer " <<iL+1 << " " << w <<std::endl;
      if(iL==0)radLength[iL]=w;
      else radLength[iL]=radLength[iL-1] + w;
      }
    }
   

    for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles
      p_pdgID->Fill((*genvec)[iP].pdgid());
      p_en->Fill((*genvec)[iP].eta());
    }//loop on gen particles

    double en_profile[30];
    for(unsigned iL(0); iL<30; iL++){
       en_profile[iL]=0;
    }

    unsigned prevLayer = 10000;
    DetectorEnum type = DetectorEnum::FECAL;
    unsigned subdetLayer=0;

    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*simhitvec)[iH];

      //discard some si layers...
      if (lHit.silayer() >= 2) continue; 

      unsigned layer = lHit.layer();

      if (layer >= nLayers) {
	//std::cout << " WARNING! SimHits with layer " << layer << " outside of detector's definition range ! Please fix the digitiser or the detector definition used here. Ignoring..." << std::endl;
	continue;
      }
      if (layer != prevLayer){
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	type = subdet.type;
	subdetLayer = layer-subdet.layerIdMin;
	prevLayer = layer;
	if (debug > 1) std::cout << " - layer " << layer << " " << subdet.name << " " << subdetLayer << std::endl;
      }     

      unsigned sec =  myDetector.getSection(layer);

      double energy = lHit.energy()*mycalib.MeVToMip(layer);
      en_profile[layer] += energy;
 }
    for(unsigned iL(0); iL<30; iL++){
       p_longpro->Fill(radLength[iL], en_profile[iL]);
    }
}

  outputFile->Write();
  outputFile->Close();

 
  return 0;


}//main
