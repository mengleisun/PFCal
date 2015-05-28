#include "HadEnergy.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"

HadEnergy::HadEnergy(double genEn,
                     double eta,
                     HGCSSDetector & myDetector,
                     const unsigned versionNumber,
                     bool isCalibed){
  genEn_ = genEn;  
  eta_ = eta;
  version_ = versionNumber;
  nSections_ = myDetector.nSections();
  nLayers_ = myDetector.nLayers();

  FHtoEslope_ = 1.; 
  FHtoEoffset_ = 0.;
  BHtoEslope_ = 1.;
  BHtoEoffset_ = 0.;
  ECALslope_ = 1.;
  ECALoffset_ = 0.;

  isCalibed_ = isCalibed;
  setabsweight();
}


HadEnergy::~HadEnergy(){
}


void HadEnergy::bookHist(TFile *outputFile){
     outputFile->cd();
 
     std::ostringstream lName;
     lName.str("");
     lName << "p_spectrum_e" << genEn_;
     p_spectrum = new TH1F(lName.str().c_str(),"p_spectrum",1000,0,250);

     lName.str("");
     lName << "p_spectrumByLayer_e" << genEn_;
     p_spectrumByLayer = new TH2F(lName.str().c_str(),"p_spectrumByLayer",nLayers_,0,nLayers_-1,100,0,5);

     lName.str("");
     lName << "p_spectrum_hightail_e" << genEn_;
     p_spectrum_hightail = new TH1F(lName.str().c_str(),"p_spectrum_hightail",1000,0,250);

     lName.str("");
     lName << "p_spectrum_lowtail_e" << genEn_;
     p_spectrum_lowtail = new TH1F(lName.str().c_str(),"p_spectrum_lowtail",1000,0,250);

     lName.str("");
     lName << "p_EE_e" << genEn_;
     p_EE = new TH1F(lName.str().c_str(),"p_EE",1000,0, isCalibed_==true? genEn_*20: genEn_*1000);

     lName.str("");
     lName << "p_EFHCAL_e" << genEn_;
     p_EFHCAL = new TH1F(lName.str().c_str(),"p_EFHCAL",1000,0,isCalibed_==true? genEn_*20: genEn_*1000);

     lName.str("");
     lName << "p_EBHCAL_e" << genEn_;
     p_EBHCAL = new TH1F(lName.str().c_str(),"p_EBHCAL",1000,0,isCalibed_==true? genEn_*20: genEn_*1000);

     lName.str("");
     lName << "p_FHcalVsBHcal_e" << genEn_;
     p_FHcalVsBHcal = new TH2F(lName.str().c_str(),"p_FHcalVsBHcal",1000,0,isCalibed_==true? genEn_*20: genEn_*1000,1000,0,isCalibed_==true? genEn_*20: genEn_*1000);

     unsigned p_size = LimMIP_.size();
     p_FHcalVsC.reserve(p_size);

     for(unsigned iLim(0); iLim < p_size; iLim++){
         lName.str("");
         lName << "p_FHcalVsC_e" << genEn_ << "_" << iLim;
         p_FHcalVsC[iLim] =  new TH2F(lName.str().c_str(),";C_{global}; EFHcal",700,0.8,1.5,1000,0,isCalibed_==true? genEn_*10: genEn_*1000);
     } 
}


bool HadEnergy::fillEnergies(const std::vector<HGCSSSamplingSection> * ssvec,
			     const std::vector<HGCSSRecoHit> * rechitvec,
			     HGCSSDetector & myDetector){

  p_spectrum->Reset();
  double EE(0), EFHCAL(0), EBHCAL(0);

  int nSec = nSections_;
  double EmipMean[nSec];
  double EmipMeanFH(0);
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
//    double absweight = 1; 
//    if(version_ == 30 || version_ == 33)absweight = absweight_[layer];
//    else if(version_ == 31)absweight = absweight_[layer + 28];
//    else if(version_ == 32)absweight = absweight_[layer + 40];
    unsigned sec =  myDetector.getSection(layer); 

    double energy = lHit.energy();

    p_spectrumByLayer->Fill(layer,energy);
 
    recSum[sec] += energy*absweight;
//    recSum[sec] += energy;
    EmipMean[sec] += energy;
    nhits[sec] += 1;

    DetectorEnum type = myDetector.detType(sec);
    if(type == DetectorEnum::FHCAL)p_spectrum->Fill(energy);
  }//loop on hits

  int HitFH(0);
  for(unsigned iS(0); iS < nSec; iS++){
    //get total energy for each sub-detector
    DetectorEnum type = myDetector.detType(iS); 
    if(type == DetectorEnum::FECAL || type == DetectorEnum::MECAL || type == DetectorEnum::BECAL)
       EE += (recSum[iS]-ECALoffset_)/ECALslope_ ;
    else if(type == DetectorEnum::FHCAL){
       EFHCAL += (recSum[iS]-FHtoEoffset_)/FHtoEslope_;
       HitFH += nhits[iS];
       EmipMeanFH += EmipMean[iS];
    }
    else if(type == DetectorEnum::BHCAL1 || type == DetectorEnum::BHCAL2)
       EBHCAL += (recSum[iS]-BHtoEoffset_)/BHtoEslope_; 
    else {
       std::cout << "the subdetector type is not defined" << std::endl;
     }
  }
  if(HitFH!=0)EmipMeanFH = EmipMeanFH/HitFH; 

  p_EE->Fill(EE);
  p_EFHCAL->Fill(EFHCAL);
  p_EBHCAL->Fill(EBHCAL);

  EE_ = EE;
  EFHCAL_ = EFHCAL;
  EBHCAL_ = EBHCAL;

  double GenE = genEn_*cosh(eta_);
  if(EFHCAL < 0.9*GenE/1.25)p_spectrum_lowtail->Add(p_spectrum);
  else if(EFHCAL > 1.1*GenE/1.25)p_spectrum_hightail->Add(p_spectrum);
  p_FHcalVsBHcal->Fill(EFHCAL, EBHCAL);  

  for(unsigned iLim(0); iLim < LimMIP_.size(); iLim++){
     double globalC = calcGlobalC(LimMIP_[iLim], EmipMeanFH);
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
  double max = genEn_<5? genEn_*cosh(eta_)*0.6: genEn_<10? genEn_*cosh(eta_)*0.8: genEn_*cosh(eta_)*0.9;
  TF1 *f1 = new TF1("f1","pol1",0, max);
  prof->Fit(f1,"R");
  double slope = f1->GetParameter(1);

  return slope;
}

void HadEnergy::setabsweight(){
   absweight_.clear();
   absweight_.reserve(52);
   absweight_[0]=1;
   absweight_[1]=1.00258;
   absweight_[2]=0.984423;
   absweight_[3]=1.00258;
   absweight_[4]=0.984423;
   absweight_[5]=1.00258;
   absweight_[6]=0.984423;
   absweight_[7]=1.00258;
   absweight_[8]=0.984423;
   absweight_[9]=1.00258;
   absweight_[10]=1.33536;
   absweight_[11]=1.3627;
   absweight_[12]=1.33536;
   absweight_[13]=1.3627;
   absweight_[14]=1.33536;
   absweight_[15]=1.3627;
   absweight_[16]=1.33536;
   absweight_[17]=1.3627;
   absweight_[18]=1.33536;
   absweight_[19]=1.3627;
   absweight_[20]=1.9495;
   absweight_[21]=1.9629;
   absweight_[22]=1.9495;
   absweight_[23]=1.9629;
   absweight_[24]=1.9495;
   absweight_[25]=1.9629;
   absweight_[26]=1.9495;
   absweight_[27]=2.01643;
   absweight_[28]=6.00121;
   absweight_[29]=5.31468;
   absweight_[30]=5.31468;
   absweight_[31]=5.31468;
   absweight_[32]=5.31468;
   absweight_[33]=5.31468;
   absweight_[34]=5.31468;
   absweight_[35]=5.31468;
   absweight_[36]=5.31468;
   absweight_[37]=5.31468;
   absweight_[38]=5.31468;
   absweight_[39]=5.31468;
   absweight_[40]=8.71728;
   absweight_[41]=8.00569;
   absweight_[42]=8.00569;
   absweight_[43]=8.00569;
   absweight_[44]=8.00569;
   absweight_[45]=8.00569;
   absweight_[46]=8.00569;
   absweight_[47]=8.00569;
   absweight_[48]=8.00569;
   absweight_[49]=8.00569;
   absweight_[50]=8.00569;
   absweight_[51]=8.00569;
}
