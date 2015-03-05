#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
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

#include "HGCSSInfo.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"

#include "HadEnergy.hh"

#include "TRandom3.h"

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

  unsigned genEn[]={30,40,50,60};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  bool concept = true;

  bool selectEarlyDecay = true;

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
	    << " -- Input energy is: " << genEn << " GeV." << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  /////////////////////////////////////////////////////////////
  //input file format
  ///////////////////////////////////////////////////////////////

  std::size_t begin = simFileName.find("_e")+2;
  std::size_t middle = simFileName.find(".root");
  std::size_t end = simFileName.find_last_of(".root")+1;
  std::string simHeader = simFileName.substr(0,begin);
  std::string simAppend = simFileName.substr(middle,end-middle);

  std::size_t begin_rec = recoFileName.find("_e")+2;
  std::size_t middle_rec = recoFileName.find(".root");
  std::size_t end_rec = recoFileName.find_last_of(".root")+1;
  std::string recHeader = recoFileName.substr(0,begin_rec);
  std::string recAppend = recoFileName.substr(middle_rec,end_rec-middle_rec);



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

  std::vector<double> EFHCAL[nGenEn];
  std::vector<double> EBHCAL[nGenEn];
  std::vector<double> GlobalC[nGenEn];
  double FvsBslope[nGenEn];

  for(unsigned iGen(0); iGen < nGenEn; iGen++){

    /////////////////////////////////////////////////////////////
    //input
    /////////////////////////////////////////////////////////////
 
    std::ostringstream input;
    input << filePath << "/" << simHeader << genEn[iGen] << simAppend;
   
    TFile *simFile = TFile::Open(input.str().c_str());
  
    if (!simFile) {
      std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
      return 1;
    }
    else std::cout << " -- input file " << simFile->GetName() << " successfully opened." << std::endl;
    
    TTree *lSimTree = (TTree*)simFile->Get("HGCSSTree");
    if (!lSimTree){
      std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
      return 1;
    }
  
    input.str("");
    input << filePath << "/" << recHeader << genEn[iGen] << recAppend;;
    
    TFile *recFile = TFile::Open(input.str().c_str());
  
    if (!recFile) {
      std::cout << " -- Error, input file " << input.str() << " cannot be opened. Exiting..." << std::endl;
      return 1;
    }
    else std::cout << " -- input file " << recFile->GetName() << " successfully opened." << std::endl;
  
    TTree *lRecTree = (TTree*)recFile->Get("RecoTree");
    if (!lRecTree){
      std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
      return 1;
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
    
    
    HadEnergy myhadReso(genEn[iGen], myDetector, versionNumber);
    myhadReso.addLimMIP(5.5);
    myhadReso.bookHist(outputFile);
  
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

      EFHCAL[iGen].push_back(myhadReso.getEFHCAL());
      EBHCAL[iGen].push_back(myhadReso.getEBHCAL());
      GlobalC[iGen].push_back(myhadReso.getGlobalC(5.5));

    }//loop on events
    
    FvsBslope[iGen] =  myhadReso.getSlope();

  }//loop over GenEn 

  outputFile->cd();
 
  double GenEn[nGenEn];
  for(unsigned iG(0); iG <nGenEn; iG++){
     GenEn[iG] = genEn[iG];
  } 
  TGraph * slopesummary = new TGraph(nGenEn, GenEn, FvsBslope);
  double slopeMean = fabs(slopesummary->GetMean(2));

  TH1F* p_totalE[nGenEn];
  std::ostringstream lName;
  for(unsigned iG(0); iG < nGenEn; iG++){
     lName.str("");
     lName << "p_totalE_e" << genEn[iG];
     p_totalE[iG] = new TH1F(lName.str().c_str(),"p_totalE",1000,0,genEn[iG]*200);
     
     unsigned nevt = EFHCAL[iG].size();
     double showerE(0);
     for(unsigned ievt(0); ievt < nevt; ievt++){
         showerE = GlobalC[iG][ievt]*(EFHCAL[iG][ievt]+(1./slopeMean *EBHCAL[iG][ievt]));
         p_totalE[iG]->Fill(showerE);
      }
  }
  
  slopesummary->Write();
  outputFile->Write();
  //outputFile->Close();
  
  return 0;
  
  
}//main
