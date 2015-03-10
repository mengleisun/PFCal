#include "HadEnergy.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"

HadEnergy::HadEnergy(double genEn,
                     HGCSSDetector & myDetector,
                     const unsigned versionNumber){
  genEn_ = genEn;  
  version_ = versionNumber;
  nSections_ = myDetector.nSections();
  nLayers_ = myDetector.nLayers();
}


HadEnergy::~HadEnergy(){
}


void HadEnergy::bookHist(TFile *outputFile){
     outputFile->cd();
 
     std::ostringstream lName;
     lName.str("");
     lName << "p_spectrum_e" << genEn_;
     p_spectrum = new TH1F(lName.str().c_str(),"p_spectrum",1000,0,500);

     lName.str("");
     lName << "p_spectrum_hightail_e" << genEn_;
     p_spectrum_hightail = new TH1F(lName.str().c_str(),"p_spectrum_hightail",1000,0,500);

     lName.str("");
     lName << "p_spectrum_lowtail_e" << genEn_;
     p_spectrum_lowtail = new TH1F(lName.str().c_str(),"p_spectrum_lowtail",1000,0,500);

     lName.str("");
     lName << "p_FHcalVsBHcal_e" << genEn_;
     p_FHcalVsBHcal = new TH2F(lName.str().c_str(),"p_FHcalVsBHcal",1000,0,30000,1000,0,30000);

     unsigned p_size = LimMIP_.size();
     p_FHcalVsC.reserve(p_size);

     for(unsigned iLim(0); iLim < p_size; iLim++){
         lName.str("");
         lName << "p_FHcalVsC_e" << genEn_ << "_" << iLim;
         p_FHcalVsC[iLim] =  new TH2F(lName.str().c_str(),";C_{global}; EFHcal",1000,0,2,1000,0,30000);
     } 
}


bool HadEnergy::fillEnergies(const std::vector<HGCSSSamplingSection> * ssvec,
			     const std::vector<HGCSSRecoHit> * rechitvec,
			     HGCSSDetector & myDetector){

  p_spectrum->Reset();
  double EE(0), EFHCAL(0), EBHCAL(0);

  int nSec = nSections_;
  double EmipMean[nSec];
  int    nhits[nSec]; 
  double recSum[nSec];
  for(unsigned iS(0); iS < nSec; iS++){
    recSum[iS] = 0;
    EmipMean[iS] = 0;
    nhits[iS] = 0;
  }

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
    const HGCSSRecoHit lHit = (*rechitvec)[iH];
    
    unsigned layer = lHit.layer();
    if (layer >= nLayers_) {
       std::cout << " WARNING! RecoHits with layer " << layer << " outside of detector's definition range ! Please fix the digitiser or the detector definition used here. Ignoring..." << std::endl;
       continue;
    }
      
    double absweight = (*ssvec)[layer].volX0trans()/(*ssvec)[0].volX0trans();
    unsigned sec =  myDetector.getSection(layer);

    double energy = lHit.energy();
 
    recSum[sec] += energy*absweight;
    EmipMean[sec] += energy;
    nhits[sec] += 1;

    p_spectrum->Fill(energy);
  }//loop on hits
  
  for(unsigned iS(0); iS < nSec; iS++){
    if(nhits[iS]!=0)EmipMean[iS] = EmipMean[iS]/nhits[iS];
  }
  EE = recSum[0] + recSum[1] + recSum[2];
  EFHCAL = recSum[3];
  EBHCAL = recSum[4]; 

  double GenMIP = genEn_ * 100;
  if(EFHCAL < GenMIP*0.95)p_spectrum_hightail->Add(p_spectrum);
  else if(EFHCAL > GenMIP*1.05)p_spectrum_lowtail->Add(p_spectrum);
  p_FHcalVsBHcal->Fill(EFHCAL, EBHCAL);  

  EE_ = EE;
  EFHCAL_ = EFHCAL;
  EBHCAL_ = EBHCAL;
  
  for(unsigned iLim(0); iLim < LimMIP_.size(); iLim++){
     double globalC = calcGlobalC(LimMIP_[iLim], EmipMean[0]);
     Cglobal_.push_back(globalC); 
     p_FHcalVsC[iLim]->Fill(globalC, EFHCAL);
  }

}

double HadEnergy::calcGlobalC(double LimMIP, double EmipMean){
  Int_t binLim = p_spectrum->GetXaxis()->FindBin(LimMIP);
  Int_t binAve = p_spectrum->GetXaxis()->FindBin(EmipMean);
  float countLim(0);
  for(unsigned iB(0); iB < binLim; iB++){
     countLim += p_spectrum->GetBinContent(iB);
  }
  float countAve(0);
  for(unsigned iB(0); iB < binAve; iB++){
     countAve += p_spectrum->GetBinContent(iB);
  }
  double globalC(0);
  if(countAve!=0)globalC = countLim/countAve;
  
  return globalC;
}

double HadEnergy::getSlope(){
  TProfile *prof = p_FHcalVsBHcal->ProfileX();
  TF1 *f1 = new TF1("f1","pol1",0, genEn_*100*0.6);
  prof->Fit(f1,"R");
  double slope = f1->GetParameter(1);

  return slope;
}
