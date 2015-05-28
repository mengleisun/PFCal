#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"

double E(const unsigned pT, const double eta){
  return pT*cosh(eta);
};

int plotHvsE(){//main

  unsigned rebinSim = 1;
  const unsigned nV = 2;
  TString version[nV] = {"31","31"};
  double eta = 2.0;
  
  std::string unit = "MIPs";
  const char* unitStr = unit.c_str();

  std::ostringstream saveName;

  bool doCalib = true;
  unsigned genEn[]={5,10,30,50};
  const unsigned nGenEn=sizeof(genEn)/sizeof(unsigned);
  unsigned rebin[20] = {4,4,4,6,6,
			6,6,8,8,10,
			10,10,100,100,100,
			6,6,6,6,6};

  //canvas so they are created only once
  TCanvas *mycE[nV];
  unsigned nx=0,ny=0;
  if (nGenEn>12) {nx=5;ny=3;}
  else if (nGenEn > 10)
    {nx=4;ny=3;}
  else if (nGenEn > 6)
    {nx=5;ny=2;}
  else if (nGenEn > 4)
    {nx=3;ny=2;}
  else if (nGenEn > 2)
    {nx=2;ny=2;}
  else 
    {nx=nGenEn;ny=1;}

  std::ostringstream lName;
  for(unsigned iV(0); iV < nV; iV++){
    lName.str("");
    lName << "mycE_" << iV;
    mycE[iV] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    mycE[iV]->Divide(nx,ny);
  }

  TCanvas *myc[2];
  for(unsigned iC(0); iC < 3; iC++){
    lName.str("");
    lName << "myc" << iC;
    myc[iC]= new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
  }

  TString plotDir = "../PLOTS/Hadron/";
      
  TFile *inputFile[nV];
  std::ostringstream linputStr;
  for(unsigned iV(0); iV < nV; iV++){
    linputStr.str("");
    linputStr << plotDir << "e_version" + version[iV] + "_model2_BON_eta2.000.root";
    inputFile[iV] = TFile::Open(linputStr.str().c_str());
    if (!inputFile[iV]) {
      std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Exiting..." << std::endl;
      return 1;
    }
    else std::cout << " -- File " << inputFile[iV]->GetName() << " sucessfully opened." << std::endl;
  }

  TH1F *p_Etotal[nV][nGenEn];
  for (unsigned iE(0); iE<nGenEn; ++iE){
    std::cout << "- Processing energy : " << genEn[iE] << std::endl;

    for(unsigned iV(0); iV < nV; iV++){	  
	lName.str("");
        if(doCalib){
	  if(version[iV]=="30")lName << "p_EE";
	  else if(version[iV]=="31")lName << "p_EFHCAL";
	  else if(version[iV]=="32")lName << "p_EBHCAL";
        } 
        else lName << "p_totalE";
        lName << "_e" << genEn[iE];
	p_Etotal[iV][iE] = (TH1F*)inputFile[iV]->Get(lName.str().c_str());
	if (!p_Etotal[iV][iE]){
	  std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null. Exiting..." << std::endl;
	  return 1;
	}
	p_Etotal[iV][iE]->Sumw2();

        std::cout << " --- Sim E = entries " << p_Etotal[iV][iE]->GetEntries() 
	    << " mean " << p_Etotal[iV][iE]->GetMean() 
	    << " rms " << p_Etotal[iV][iE]->GetRMS() 
	    << " overflows " << p_Etotal[iV][iE]->GetBinContent(p_Etotal[iV][iE]->GetNbinsX()+1)
	    << std::endl;

        if(version[iV] == "30")p_Etotal[iV][iE]->Rebin(genEn[iE]<100?8:genEn[iE]<50?16:genEn[iE]<300?20:40);
     }

  }//loop on energies

  gStyle->SetOptStat(0);

  //draw calibration curves
  TGraphErrors *calib = new TGraphErrors();
  calib->SetName("calib");
  calib->SetMarkerStyle(21);
  calib->SetTitle("");
  TGraphErrors *calibFit = (TGraphErrors *) calib->Clone("calibFit");
  TGraphErrors *reso = (TGraphErrors *) calib->Clone("reso");
  
	//simhits
  for (unsigned iE(0); iE<nGenEn; ++iE){
      std::cout << "- Processing energy : " << genEn[iE] << std::endl;
      TF1 *fitResult[nV];
      for(unsigned iV(0); iV < nV; ++iV){
  	  //plot total E
	  mycE[iV]->cd(iE+1);
	  gStyle->SetOptFit(0);
	  double eMin = p_Etotal[iV][iE]->GetMean()-10*p_Etotal[iV][iE]->GetRMS();
	  double eMax = p_Etotal[iV][iE]->GetMean()+10*p_Etotal[iV][iE]->GetRMS();
	  p_Etotal[iV][iE]->GetXaxis()->SetRangeUser(eMin,eMax);
	  p_Etotal[iV][iE]->Draw("PE");
	  char buf[500];
	  sprintf(buf,"E=%d GeV",genEn[iE]);
	  p_Etotal[iV][iE]->SetTitle(buf);
	  p_Etotal[iV][iE]->Fit("gaus","LR0","",
			    p_Etotal[iV][iE]->GetMean()-2*p_Etotal[iV][iE]->GetRMS(),
			    p_Etotal[iV][iE]->GetMean()+2*p_Etotal[iV][iE]->GetRMS());

	  fitResult[iV] = p_Etotal[iV][iE]->GetFunction("gaus");

	  p_Etotal[iV][iE]->Fit("gaus","LR+","same",
			    fitResult[iV]->GetParameter(1)-2*fitResult[iV]->GetParameter(2),
			    fitResult[iV]->GetParameter(1)+2*fitResult[iV]->GetParameter(2));
	  fitResult[iV] = p_Etotal[iV][iE]->GetFunction("gaus");
	  TLatex lat;
	  double latx = std::max(0.,p_Etotal[iV][iE]->GetMean()-5*p_Etotal[iV][iE]->GetRMS());
	  double laty = p_Etotal[iV][iE]->GetMaximum();
	  sprintf(buf,"<E> = %3.3f %s",p_Etotal[iV][iE]->GetMean(),unitStr);
	  lat.DrawLatex(latx,laty*0.9,buf);
	  sprintf(buf,"RMS = %3.3f #pm %3.1f %s",p_Etotal[iV][iE]->GetRMS(),p_Etotal[iV][iE]->GetRMSError(),unitStr);
	  lat.DrawLatex(latx,laty*0.8,buf);
	  sprintf(buf,"RMS/mean = %3.3f",p_Etotal[iV][iE]->GetRMS()/p_Etotal[iV][iE]->GetMean());
	  lat.DrawLatex(latx,laty*0.7,buf);
	  sprintf(buf,"<Efit> = %3.3f +/- %3.3f %s",fitResult[iV]->GetParameter(1),fitResult[iV]->GetParError(1),unitStr);
	  lat.DrawLatex(latx,laty*0.6,buf);
	  sprintf(buf,"RMSfit = %3.3f +/- %3.3f %s",fitResult[iV]->GetParameter(2),fitResult[iV]->GetParError(2),unitStr);
	  lat.DrawLatex(latx,laty*0.5,buf);
	  sprintf(buf,"RMS/meanfit = %3.3f",fitResult[iV]->GetParameter(2)/fitResult[iV]->GetParameter(1));
	  lat.DrawLatex(latx,laty*0.4,buf);

	  sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitResult[iV]->GetChisquare(),fitResult[iV]->GetNDF(),fitResult[iV]->GetChisquare()/fitResult[iV]->GetNDF());
          lat.DrawLatex(latx,laty*0.3,buf);
       } //loop on versions

       Int_t np=calib->GetN();
       calib->SetPoint(np,E(genEn[iE],eta),p_Etotal[1][iE]->GetMean());
       calib->SetPointError(np,0,p_Etotal[1][iE]->GetMeanError());
       calibFit->SetPoint(np,E(genEn[iE],eta),fitResult[1]->GetParameter(1));
       calibFit->SetPointError(np,0,fitResult[1]->GetParError(1));
       reso->SetPoint(np, E(genEn[iE],eta),fitResult[1]->GetParameter(2)/fitResult[1]->GetParameter(1));
       double errFit = fitResult[1]->GetParameter(2)/fitResult[1]->GetParameter(1)*sqrt(pow(fitResult[1]->GetParError(2)/fitResult[1]->GetParameter(2),2)+pow(fitResult[1]->GetParError(1)/fitResult[1]->GetParameter(1),2));
       reso->SetPointError(np,0,errFit);
   }//loop on energies


   for(unsigned iV(0); iV<nV; iV++){
       saveName.str("");
       saveName << plotDir << "/SimG4Etotal_" << version[iV];
       mycE[iV]->Update();
       mycE[iV]->Print((saveName.str()+".png").c_str());
       mycE[iV]->Print((saveName.str()+".pdf").c_str());
   }
	//draw calib

   std::string type[3];
   for(unsigned i=0; i<3 ; i++) {
      myc[i]->cd();
      type[i] = i==0 ? "calib" : i==1? "calibFit": "reso"; 

      std::cout << "- Processing type : " << type[i] << std::endl;

      TGraphErrors * gr =( i==0 ? calib : i==1? calibFit : reso);

      gr->GetXaxis()->SetLabelSize(0.06);
      gr->GetXaxis()->SetTitleSize(0.06);
      gr->GetYaxis()->SetLabelSize(0.06);
      gr->GetYaxis()->SetTitleSize(0.06);
      gr->GetXaxis()->SetTitleOffset(0.7);
      gr->GetYaxis()->SetTitleOffset(0.8);

      gr->SetTitle("");
      gr->Draw("ap");
   
      if(i<2){ 
        gr->GetXaxis()->SetTitle("Gen energy [GeV]");
        gr->GetYaxis()->SetTitle("Reco energy [MIP]");

        char buf[500];
        TF1 *fitFunc=new TF1("calib","[0]+[1]*x",0,20000);
        fitFunc->SetLineColor(1);
        gr->Fit(fitFunc,"RME");
        TLatex lat;
        lat.SetTextColor(1);
        sprintf(buf,"<HE> #propto a + b #times <EE> ");
        lat.DrawLatex(genEn[0],gr->GetYaxis()->GetXmax()*0.9,buf);
        sprintf(buf,"a = %3.3f #pm %3.3f %s",fitFunc->GetParameter(0),fitFunc->GetParError(0),unitStr);
        lat.DrawLatex(genEn[0]+i/4*50,gr->GetYaxis()->GetXmax()*(0.8-i/2*0.25),buf);
        sprintf(buf,"b = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(1),fitFunc->GetParError(1),unitStr);
        lat.DrawLatex(genEn[0]+i/4*50,gr->GetYaxis()->GetXmax()*(0.7-i/2*0.25),buf);
        sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
        lat.DrawLatex(genEn[0]+i/4*50,gr->GetYaxis()->GetXmax()*(0.6-i/2*0.25),buf);
      }
      else
      {
		TF1 *fitFunc2;
		fitFunc2 =new TF1("reso","sqrt([0]*[0]/x+[1]*[1])",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());

		fitFunc2->SetParameter(0,0.2);
		fitFunc2->SetParLimits(0,0,1);
		fitFunc2->SetParameter(1,0.01);
		fitFunc2->SetParLimits(1,0,1);

		gr->Fit(fitFunc2,"RME");
		double sigmaStoch = fitFunc2->GetParameter(0);
		double sigmaStochErr = fitFunc2->GetParError(0);
		double sigmaConst = fitFunc2->GetParameter(1);
		double sigmaConstErr = fitFunc2->GetParError(1);

	        TLatex lat;
	        char buf[500];
		lat.SetTextColor(1);
		sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E}} #oplus c");
		double Emin = 40;
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax(),buf);
		sprintf(buf,"s=%3.3f #pm %3.3f",sigmaStoch,sigmaStochErr);
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.9-i/6*0.25),buf);
		sprintf(buf,"c=%3.3f #pm %3.3f",sigmaConst,sigmaConstErr);
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.82-i/6*0.25),buf);
		sprintf(buf,"chi2/NDF = %3.3f/%d = %3.3f",fitFunc2->GetChisquare(),fitFunc2->GetNDF(),fitFunc2->GetChisquare()/fitFunc2->GetNDF());
		lat.DrawLatex(Emin,gr->GetYaxis()->GetXmax()*(0.66-i/6*0.25),buf);
	}
    }//loop on versions

    for(unsigned iC(0); iC<3; iC++){ 
      saveName.str("");
      if(iC==0)saveName << plotDir << "/EHvsEE_raw_" << version[1];
      else if(iC==1)saveName << plotDir << "/EHvsEE_fit_" << version[1];
      else if(iC==2)saveName << plotDir << "/Reso_" << version[1];
      myc[iC]->Update();
      myc[iC]->Print((saveName.str()+".png").c_str());
      myc[iC]->Print((saveName.str()+".pdf").c_str());
   }
  return 0;
  
  
}//main
