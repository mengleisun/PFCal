#ifndef HadEnergy_h
#define HadEnergy_h

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TAxis.h"
#include "TMath.h"
#include "TProfile.h"

#include "HGCSSRecoHit.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDetector.hh"

class HadEnergy{

public:
  HadEnergy(double genEn,
            double eta,
            HGCSSDetector & myDetector,
            const unsigned versionNumber=12,
            bool isCalibed=true);
  ~HadEnergy();

  inline void addLimMIP(double lim){
    if(lim >0)LimMIP_.push_back(lim);
  }

  void bookHist(TFile *outputFile);

  void setabsweight();

  inline void setOutputFile(const std::string outputName){
       outputName_ = outputName;
       TFile *outputFile_ =  TFile::Open(outputName.c_str(),"RECREATE");
       bookHist(outputFile_); 
  };
  
  inline void setFHtoE(const double FHtoEslope,const double FHtoEoffset){
       FHtoEslope_ = FHtoEslope;
       FHtoEoffset_ = FHtoEoffset;
  }

  inline void setBHtoE(const double BHtoEslope,const double BHtoEoffset){
       BHtoEslope_ = BHtoEslope;
       BHtoEoffset_ = BHtoEoffset;
  }
  
  inline void setEEcalib(const double ECALslope,const double ECALoffset){
       ECALslope_ = ECALslope;
       ECALoffset_ = ECALoffset;
  }

  bool fillEnergies(const std::vector<HGCSSSamplingSection> * ssvec,
                    const std::vector<HGCSSRecoHit> * rechitvec,
                    HGCSSDetector & myDetector);
 
  double calcGlobalC(double LimMIP, double EmipMean); 

  inline double getEE(){  return EE_;}

  inline double getEFHCAL(){  return EFHCAL_;}
  
  inline double getEBHCAL(){  return EBHCAL_;}
  
  inline double getGlobalC(unsigned iLim) { 
      if(iLim < Cglobal_.size()) return Cglobal_[iLim];
      else return 0;
  }

  double getSlope();

private:

  std::string outputName_;
  TFile *outputFile_;
  TTree *outtree_;
  
  //for tree
  double nLayers_;
  unsigned version_;
  unsigned nSections_;
  double genEn_; 
  double eta_;
  std::vector<double> absweight_;

  double EE_;
  double EFHCAL_;
  double EBHCAL_;
  std::vector<double> LimMIP_;
  std::vector<double> Cglobal_; 

  TH1F *p_spectrum;
  TH1F *p_spectrum_hightail;
  TH1F *p_spectrum_lowtail;
  TH1F *p_EE;
  TH1F *p_EFHCAL;
  TH1F *p_EBHCAL;
  TH2F *p_FHcalVsBHcal;
  TH2F *p_spectrumByLayer;
  std::vector<TH2F*> p_FHcalVsC; 

  bool isCalibed_;
  double FHtoEslope_; 
  double FHtoEoffset_;
  double BHtoEslope_;
  double BHtoEoffset_;
  double ECALslope_;
  double ECALoffset_;
};

#endif
