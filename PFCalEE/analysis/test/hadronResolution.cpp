#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TChain.h"

#include "HGCSSInfo.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"

#include "HadEnergy.hh"

#include "TRandom3.h"

bool testInputFile(std::string input, TFile* & file){
  file = TFile::Open(input.c_str());
  
  if (!file) {
    std::cout << " -- Error, input file " << input.c_str() << " cannot be opened. Skipping..." << std::endl;
    return false;
  }
  else std::cout << " -- input file " << file->GetName() << " successfully opened." << std::endl;
  return true;
};


int main(int argc, char** argv){//main  

  if (argc < 7) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <path to input files>"
	      << " <name of input sim file>"
	      << " <name of input reco file>"
	      << " <full path to output file>"
	      << " <number of si layers to consider: 1,2 or 3>" 
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  unsigned genEn[]={5,10,30,50};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  double eta = 2.00;
  unsigned nRuns = 0;

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
//  bool isCalibed = true;
  bool isCalibed = false;

  bool concept = true;

  bool selectEarlyDecay = true;

  //////////////////////////////////////////////////////////
  //// Hardcoded factor ////////////////////////////////////
  //////////////////////////////////////////////////////////
//  double FHtoEslope = 11.085; 
//  double FHtoEoffset = -4.34;
//  double BHtoEslope = 4.079;
//  double BHtoEoffset = 293.55;
//  double BHtoEoffset = 0;
//  double ECALslope = 116.897;
//  double ECALoffset = -101.775;
  
//  double HcalPionOffset = 0.; 
//  double HcalPionCalib = 0.92;

  double FHtoEslope = 1; 
  double FHtoEoffset = 0;
  double BHtoEslope = 1;
//  double BHtoEoffset = 293.55;
  double BHtoEoffset = 0;
  double ECALslope = 1;
  double ECALoffset = 0;
  
  double HcalPionOffset = 0.; 
  double HcalPionCalib = 1;
  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(1234);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string simFileName = argv[3];
  std::string recoFileName = argv[4];

  std::string inFilePath = filePath+simFileName;

  std::string outPath = argv[5];
  unsigned nSiLayers = 2;
  nSiLayers = atoi(argv[6]);

  unsigned debug = 0;
  if (argc >7) debug = atoi(argv[7]);


  bool isEM = false;

  if (inFilePath.find("e-")!=inFilePath.npos || 
      inFilePath.find("e+")!=inFilePath.npos) isEM = true;

  if (selectEarlyDecay && isEM) {
    selectEarlyDecay = false;
  }

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  /////////////////////////////////////////////////////////////
  //input file format
  ///////////////////////////////////////////////////////////////

  std::size_t begin = simFileName.find("_et")+3;
  std::size_t middle = simFileName.find("_eta");
  std::size_t end = simFileName.find_last_of(".root")+1;
  std::size_t run(0);
  if(nRuns>0)run = simFileName.find("_run");
  std::string simHeader = simFileName.substr(0,begin);
  std::string simAppend("");
  if(nRuns>0)simAppend = simFileName.substr(middle,run-middle);
  else simAppend = simFileName.substr(middle,end-middle);

  std::size_t begin_rec = recoFileName.find("_et")+3;
  std::size_t middle_rec = recoFileName.find("_eta");
  std::size_t end_rec = recoFileName.find_last_of(".root")+1;
  std::size_t run_rec(0);
  if(nRuns>0)run_rec = recoFileName.find("_run");
  std::string recHeader = recoFileName.substr(0,begin_rec);
  std::string recAppend;
  if(nRuns>0)recAppend = recoFileName.substr(middle_rec,run_rec-middle_rec);
  else recAppend = recoFileName.substr(middle_rec,end_rec-middle_rec);


  /////////////////////////////////////////////////////////////
  //output
  ///////////////////////////////////////////////////////////////

  TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }

  std::vector<double> EE[nGenEn];
  std::vector<double> EFHCAL[nGenEn];
  std::vector<double> EBHCAL[nGenEn];
  std::vector<double> GlobalC[nGenEn];
  double HvsEslope[nGenEn];
  double FvsBslope[nGenEn];

  for(unsigned iGen(0); iGen < nGenEn; iGen++){

    /////////////////////////////////////////////////////////////
    //input
    /////////////////////////////////////////////////////////////
    TFile * simFile = 0;
    TFile * recFile = 0;
 
    TChain  *lSimTree = new TChain("HGCSSTree");
    TChain  *lRecTree = new TChain("RecoTree");
  
    std::ostringstream inputsim, inputrec;
    if (nRuns == 0){
      inputsim << filePath << "/" << simHeader << genEn[iGen] << simAppend;
      inputrec << filePath << "/" << recHeader << genEn[iGen] << recAppend;;
      if (!testInputFile(inputsim.str(),simFile)) return 1;
      lSimTree->AddFile(inputsim.str().c_str());
      if (!testInputFile(inputrec.str(),recFile)) return 1;
      lRecTree->AddFile(inputrec.str().c_str());
    }
    else {
      for (unsigned i(0);i<nRuns;++i){
        inputsim.str("");
        inputsim << filePath << "/" << simHeader << genEn[iGen] << simAppend << "_run" << i << ".root";
        inputrec.str("");
        inputrec << filePath << "/" << recHeader << genEn[iGen] << recAppend << "_run" << i << ".root";
        if (!testInputFile(inputsim.str(),simFile)) return 1;
        lSimTree->AddFile(inputsim.str().c_str());
        if (!testInputFile(inputrec.str(),recFile)) return 1;
        lRecTree->AddFile(inputrec.str().c_str());
      }
    }
 
  
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
  
  
    std::cout << " -- Version number is : " << versionNumber 
  	    << ", model = " << model
  	    << ", cellSize = " << cellSize
  	    << std::endl;
  
  
    //initialise detector
    HGCSSDetector & myDetector = theDetector();
  
    myDetector.buildDetector(versionNumber,concept,isCaliceHcal);
  
    //initialise calibration class
    HGCSSCalibration mycalib(inFilePath);
    HGCSSDigitisation myDigitiser;
    myDigitiser.setRandomSeed(lRndm->GetSeed());
  
    const unsigned nLayers = myDetector.nLayers();
    const unsigned nSections = myDetector.nSections();
  
    std::cout << " -- N layers = " << nLayers << std::endl
  	    << " -- N sections = " << nSections << std::endl;
    
    HadEnergy myhadReso(genEn[iGen],eta, myDetector, versionNumber, isCalibed);
    myhadReso.addLimMIP(5.5);
    myhadReso.bookHist(outputFile);
    myhadReso.setFHtoE(FHtoEslope, FHtoEoffset);
    myhadReso.setBHtoE(BHtoEslope, BHtoEoffset);
    myhadReso.setEEcalib(ECALslope, ECALoffset);
 
    HGCSSEvent * event = 0;
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;
    
    lSimTree->SetBranchAddress("HGCSSEvent",&event);
    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    
    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  
    const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;
    
    std::cout << "- Processing = " << nEvts  << " events out of " << lSimTree->GetEntries() << std::endl;
  
    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
      if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
      else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
      lSimTree->GetEntry(ievt);
      lRecTree->GetEntry(ievt);
  
      if (debug){
        std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
      }

      myhadReso.fillEnergies(ssvec, rechitvec, myDetector); 

      EE[iGen].push_back(myhadReso.getEE());
      EFHCAL[iGen].push_back(myhadReso.getEFHCAL());
      EBHCAL[iGen].push_back(myhadReso.getEBHCAL());
      GlobalC[iGen].push_back(myhadReso.getGlobalC(5.5));

    }//loop on events
    
    FvsBslope[iGen] =  myhadReso.getSlope();
  }//loop over GenEn 

  outputFile->cd();
 
  TGraph * slopesummary = new TGraph();
  for(unsigned iG(0); iG<nGenEn; iG++){
     if(genEn[iG] > 9){
       Int_t np=slopesummary->GetN();  
       slopesummary->SetPoint(np,genEn[iG],FvsBslope[iG]);
     }
  }
  double slopeMean = fabs(slopesummary->GetMean(2));

  TH1F* p_totalE[nGenEn];
  std::ostringstream lName;
  for(unsigned iG(0); iG < nGenEn; iG++){
     lName.str("");
     lName << "p_totalE_e" << genEn[iG];
     p_totalE[iG] = new TH1F(lName.str().c_str(),"p_totalE",1000,0,genEn[iG]*200);
     
     unsigned nevt = EFHCAL[iG].size();
     double showerE(0);
     double calibedF(0),calibedB(0);
     for(unsigned ievt(0); ievt < nevt; ievt++){
         showerE = GlobalC[iG][ievt]*(EE[iG][ievt]+(EFHCAL[iG][ievt]+(1./slopeMean *EBHCAL[iG][ievt])-HcalPionOffset)/HcalPionCalib);
         p_totalE[iG]->Fill(showerE);
      }
  }
  
  slopesummary->Write();
  outputFile->Write();
  //outputFile->Close();
  
  return 0;
  
  
}//main
