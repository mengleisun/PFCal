#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
	
bool checkHistPointer(TH1 *pointer, std::string name){
  if (!pointer) {
    std::cout << " -- ERROR, pointer for histocalibam " << name << " is null. Exiting..." << std::endl;
    return false;
  }
  return true;
}

int plotEMCalib(){//main  

  const unsigned nS = 1;
  std::string scenario[nS] = {
    "pi+/"
  };
  
  const unsigned nV = 2;
  TString version[nV] = {"30","31"};
  unsigned versionN[nV] = {30,31};

  TString pSuffix = "";

  unsigned genEn[]={5,10,30};

  unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);

  double SectorE[nV][nGenEn];
  for(unsigned iV(0); iV<nV; iV++){
    for(unsigned iE(0); iE< nGenEn; iE++){
        SectorE[iV][iE] = 0;
    }
  }

  const unsigned limRef = 2;

  std::ostringstream lName;

  const unsigned nCtot = 2;
  TCanvas *myc[nCtot];
  for (unsigned iC(0); iC<nCtot;++iC){
    lName.str("");
    lName << "Canvas_" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1000,700);
  }
  myc[0]->Divide(3,3);

  //draw calibration curves
  TGraphErrors *calib = new TGraphErrors();
  calib->SetName("calib");
  calib->SetMarkerStyle(21);
  calib->SetTitle("");

  std::ostringstream saveName;
  
  for (unsigned iV(0); iV<nV;++iV){//loop on versions
    for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
    
      TString plotDir = "../PLOTS/gitV00-03-03/version"+version[iV]+"/"+scenario[iS]+"/";
      
      TH1F *p_EE[nGenEn];
      TH1F *p_EFHCAL[nGenEn];
      TH1F *p_EBHCAL[nGenEn];

      for (unsigned iE(0); iE<nGenEn; ++iE){

	std::cout << " -- Processing energy " << genEn[iE] << std::endl;

	lName.str("");
	lName << plotDir << "hadron.root";
	TFile *inputFile = TFile::Open(lName.str().c_str());
	if (!inputFile) {
	  std::cout << " -- Error, input file " << lName.str() << " cannot be opened. Going to next..." << std::endl;
	  continue;
	  //return 1;
	}

	lName.str("");
	lName << "p_EE_e" << genEn[iE];
	p_EE[iE] = (TH1F*)gDirectory->Get(lName.str().c_str());
	if (!checkHistPointer(p_EE[iE],lName.str())) return 1;

	lName.str("");
	lName << "p_EFHCAL_e" << genEn[iE];
	p_EFHCAL[iE] = (TH1F*)gDirectory->Get(lName.str().c_str());
	if (!checkHistPointer(p_EFHCAL[iE],lName.str())) return 1;

	lName.str("");
	lName << "p_EBHCAL_e" << genEn[iE];
	p_EBHCAL[iE] = (TH1F*)gDirectory->Get(lName.str().c_str());
	if (!checkHistPointer(p_EBHCAL[iE],lName.str())) return 1;

	myc[0]->cd(iE+1);

        switch(versionN[iV]){
        case 30:
        {
	    p_EE[iE]->Fit("gaus","LR+","",
			p_EE[iE]->GetMean()-2*p_EE[iE]->GetRMS(),
			p_EE[iE]->GetMean()+2*p_EE[iE]->GetRMS());
	    TF1 *fitResult = p_EE[iE]->GetFunction("gaus");
            SectorE[iV][iE] = fitResult->GetParameter(1); 
        }
        case 31:
        {
	    p_EE[iE]->Fit("gaus","LR+","",
			p_EE[iE]->GetMean()-2*p_EE[iE]->GetRMS(),
			p_EE[iE]->GetMean()+2*p_EE[iE]->GetRMS());
	    TF1 *fitResult = p_EE[iE]->GetFunction("gaus");
            SectorE[iV][iE] = fitResult->GetParameter(1);
        }
        case 32:
        {
	    p_EE[iE]->Fit("gaus","LR+","",
			p_EE[iE]->GetMean()-2*p_EE[iE]->GetRMS(),
			p_EE[iE]->GetMean()+2*p_EE[iE]->GetRMS());
	    TF1 *fitResult = p_EE[iE]->GetFunction("gaus");
	    SectorE[iV][iE] = fitResult->GetParameter(1);
        }   
        }

      }//loop on energies
 
    }//loop on scenarios
    
  }//loop on versions
  
  myc[1]->cd();
  for(unsigned iE(0); iE<nGenEn; iE++){
      Int_t np = calib->GetN();
      calib->SetPoint(np,SectorE[0][iE],SectorE[1][iE]);
  }
  calib->GetXaxis()->SetLabelSize(0.06);
  calib->GetXaxis()->SetTitleSize(0.06);
  calib->GetYaxis()->SetLabelSize(0.06);
  calib->GetYaxis()->SetTitleSize(0.06);
  calib->GetXaxis()->SetTitleOffset(0.7);
  calib->GetYaxis()->SetTitleOffset(0.8);
  calib->SetTitle("");
  calib->GetXaxis()->SetTitle("ECAL energy [GeV]");
  calib->GetYaxis()->SetTitle("HCAL energy [GeV]");
  calib->Draw();
  calib->Fit("pol1");
  TF1 *fitFunc = calib->GetFunction("pol1");
  char buf[500];
  TLatex lat;
  lat.SetTextColor(6);
  sprintf(buf,"<E> #propto E(a + b #times E)");
  lat.DrawLatex(genEn[0],calib->GetYaxis()->GetXmax()*0.9,buf);
  sprintf(buf,"a = %3.3f #pm %3.3f ",fitFunc->GetParameter(0),fitFunc->GetParError(0));
  lat.DrawLatex(genEn[0],calib->GetYaxis()->GetXmax()*(0.8),buf);
  sprintf(buf,"b = %3.3f #pm %3.3f",fitFunc->GetParameter(1),fitFunc->GetParError(1));
  lat.DrawLatex(genEn[0],calib->GetYaxis()->GetXmax()*(0.7),buf);
  sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatex(genEn[0],calib->GetYaxis()->GetXmax()*(0.5),buf);

  saveName.str("");
  saveName << "../PLOTS/EbeamvsEshower";
  myc[1]->Update();
  myc[1]->Print((saveName.str()+".png").c_str());
  myc[1]->Print((saveName.str()+".pdf").c_str());

  return 0;


}//main
