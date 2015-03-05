#include <TPaveStats.h>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TInterpreter.h>
#include "TH1F.h"
#include <TLatex.h>
#include <TFile.h>
#include <TLatex.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#define PI 3.14159265358979

void drawGap()
{
 TCanvas *myc;
 myc= new TCanvas("resolution_eta23","Electron Resolution",0,0,600,600);

 myc->cd();
 unsigned mips[3] = {4,8,20};
 unsigned version[2]={12,13};

 double slope[2][3];
 double offset[2][3];
 double noiseterm[2][3];

 for(unsigned iversion(0); iversion<2; iversion++){
 for(unsigned ithres(0); ithres<3; ithres++){

     double noisevalue[2][7][7];
     double noiseerr[2][7][7];
     std::ifstream noiselog;
     std::ostringstream noisename;
     noisename << "../data/Noiseterm_tr" << mips[ithres] << "_v" << version[iversion] << ".dat";
     noiselog.open(noisename.str().c_str());
    std::cout << " successfully open  " << noisename.str().c_str() << std::endl;
    for(unsigned iPu(0); iPu<2; iPu++){
        for(unsigned iEta(0); iEta<7; iEta++){
            for(unsigned iSR(0); iSR<7; iSR++){
               noiselog >> noisevalue[iPu][iEta][iSR] >> noiseerr[iPu][iEta][iSR];
            }
         }
     }

    std::ostringstream filename;
    filename << "fit_tr" << mips[ithres] << "_v" << version[iversion] << ".root";
    TFile *Reso_collection = TFile::Open(filename.str().c_str());
    TGraphErrors *reso;
    TLegend *leg;
    std::ostringstream rootFileName;
    std::ostringstream graphName;
    std::ostringstream legendName;

    graphName.str("");
    graphName << "resoRecoFit2eta23pu1"; 
    reso= (TGraphErrors *)Reso_collection->Get(graphName.str().c_str());
    TF1 *resoModel=new TF1("resomodel","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0,500);
    resoModel->SetParameter(0,0.2);
    resoModel->SetParLimits(0,0,1);
    resoModel->SetParameter(1,0);
    resoModel->SetParLimits(1,0,1);
    resoModel->SetParLimits(2,0,2);
    resoModel->SetParameter(2,noisevalue[1][3][2]);
    resoModel->SetLineColor(iversion+1);

   reso->GetYaxis()->SetRangeUser(0.00,0.10);
   reso->GetYaxis()->SetLabelSize(0.03);
   reso->GetXaxis()->SetLabelSize(0.03);
   reso->SetLineColor(iversion+1);
   reso->SetMarkerStyle(ithres+20);
   reso->SetMarkerSize(1.5);
   reso->SetMarkerColor(iversion+1);
   reso->Fit(resoModel,"MER+");

    slope[iversion][ithres]=resoModel->GetParameter(0);
    offset[iversion][ithres]=resoModel->GetParameter(1);
    noiseterm[iversion][ithres]=resoModel->GetParameter(2);
 
     leg = new TLegend(0.35, 0.8-ithres*0.06-iversion*0.18, 0.9, 0.86-ithres*0.06-iversion*0.18);
     leg->AddEntry(reso,Form("%dmm,%dMips #frac{#sigma}{E} #propto #frac{%3.3f}{#sqrt{E}} + #frac{%3.3f}{E} + %3.3f",version[iversion]==12?2:4,mips[ithres]/4,slope[iversion][ithres],noiseterm[iversion][ithres], offset[iversion][ithres]),"p");
     leg->SetTextSize(0.03);
     leg->SetBorderSize(0);
     leg->SetFillColor(0);

     if(iversion==0 && ithres==0){
	 reso->Draw("ap");
	 leg->Draw("ap");
     } else {
     reso->Draw("psame");
     leg->Draw("psame");
     }
  }
}

          
     myc->SaveAs("resolution_eta23.png");


}
