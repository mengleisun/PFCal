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
#include "TH2F.h"
#include <TLatex.h>
#include <TFile.h>
#include <TLatex.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#define PI 3.14159265358979

void drawAllEta()
{
 TCanvas *myc[5];
 std::ostringstream canvasName;
 for(unsigned iSR(0);iSR<5;iSR++){
     canvasName.str("");
     canvasName << "resolution_SR" << iSR;
     myc[iSR] = new TCanvas(canvasName.str().c_str(),"Electron Resolution",0,0,1200,600);
     myc[iSR]->Divide(2); 
 }

 unsigned pu[2] = {0,140};
 unsigned eta[]={17,19,21,23,25,27,29};
 const unsigned netaEn=sizeof(eta)/sizeof(unsigned);

 double slope[netaEn];
 double slopeErr[netaEn];
 double offset[netaEn];
 double offsetErr[netaEn];
 double noiseterm[netaEn];
 double noiseErr[netaEn];

 double noisevalue[2][7][7];
 double noiseerr[2][7][7];
 std::ifstream noiselog;
 std::ostringstream noisename;
 noisename << "../data/Noiseterm_tr20_v13.dat";
 noiselog.open(noisename.str().c_str());
 for(unsigned iPu(0); iPu<2; iPu++){
    for(unsigned iEta(0); iEta<7; iEta++){
       for(unsigned iSR(0); iSR<7; iSR++){
          noiselog >> noisevalue[iPu][iEta][iSR] >> noiseerr[iPu][iEta][iSR];
       }
    }
 }


 TFile *Reso_collection = TFile::Open("fit_tr20_v13.root");
 TGraphErrors *reso[netaEn];
 TLegend *leg[netaEn];
 std::ostringstream rootFileName;
 std::ostringstream graphName;
 std::ostringstream legendName;
 for(unsigned iSR(0); iSR<5; iSR++){
     for(unsigned ipu(0);ipu<2;ipu++){
             myc[iSR]->cd(ipu+1);
             std::ostringstream mgName("");
             mgName << "(4mm air gap)  PU" << pu[ipu] << " SR" << iSR;
             if(iSR==0)mgName << "(1x1 cm^2)";
             else if(iSR==1)mgName << "(2x2 cm^2);E (GeV); ";
             else if(iSR==2)mgName << "(3x3 cm^2);E (GeV); ";
             else if(iSR==3)mgName << "(4x4 cm^2);E (GeV); ";
             else if(iSR==4)mgName << "(5x5 cm^2);E (GeV); ";
             TH2F *frame = new TH2F(mgName.str().c_str(),mgName.str().c_str(),300,0,1500,100,0,0.1);
             frame->SetStats(0);
             frame->Draw();
	     for(unsigned iEta(0); iEta < netaEn; iEta++){
		 graphName.str("");
		 graphName << "resoRecoFit" << iSR << "eta" << eta[iEta] << "pu" << ipu; 
                 std::cout << "try to get graph:" << graphName.str().c_str() << std::endl;
		 reso[iEta] = (TGraphErrors *)Reso_collection->Get(graphName.str().c_str());
//		 TF1 *resoModel=new TF1("resomodel","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0,500);
		 TF1 *resoModel=new TF1("resomodel","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0,1000);
		 resoModel->SetParameter(0,0.2);
		 resoModel->SetParLimits(0,0,1);
		 resoModel->SetParameter(1,0);
		 resoModel->SetParLimits(1,0,1);
                 resoModel->SetParLimits(2,0,10);
                 resoModel->SetParameter(2,noisevalue[ipu][iEta][iSR]);
                 resoModel->SetLineColor(iEta+1);

	  //       reso[iEta]->SetTitle(mgName.str().c_str());
           //      reso[iEta]->GetXaxis()->SetRangeUser(0.,500.);
                 reso[iEta]->GetYaxis()->SetRangeUser(0.00,0.10);
                 reso[iEta]->GetYaxis()->SetLabelSize(0.03);
                 reso[iEta]->GetXaxis()->SetLabelSize(0.03);
                 reso[iEta]->GetXaxis()->SetRange(0,1000);
                 reso[iEta]->SetLineColor(iEta+1);
		 reso[iEta]->SetMarkerStyle(iEta%4 +20);
		 reso[iEta]->SetMarkerSize(1.5);
		 reso[iEta]->SetMarkerColor(iEta+1);
 //	         reso[iEta]->Fit(resoModel,"MER+");

//		 slope[iEta]=resoModel->GetParameter(0);
//		 slopeErr[iEta]=resoModel->GetParError(0);
//		 offset[iEta]=resoModel->GetParameter(1);
//		 offsetErr[iEta]=resoModel->GetParError(1);
//                 noiseterm[iEta]=resoModel->GetParameter(2);
//                 noiseErr[iEta]=resoModel->GetParError(2);
 
		 leg[iEta] = new TLegend(0.7, 0.8-iEta*0.06, 0.9, 0.86-iEta*0.06);
//		 leg[iEta]->AddEntry(reso[iEta],Form("eta=%1.1f #frac{#sigma}{E} #propto #frac{%3.3f}{#sqrt{E}} + #frac{%3.3f}{E} + %3.3f",(double)eta[iEta]/10.,slope[iEta],noiseterm[iEta], offset[iEta]),"p");
		 leg[iEta]->AddEntry(reso[iEta],Form("eta=%1.1f",(double)eta[iEta]/10.),"p");
		 leg[iEta]->SetTextSize(0.03);
		 leg[iEta]->SetBorderSize(0);
		 leg[iEta]->SetFillColor(0);

	         reso[iEta]->Draw("psame");
                 leg[iEta]->Draw("psame");
	      }
          myc[iSR]->cd(ipu+1);
          
     }
     canvasName.str("");
     canvasName << "resolution_SR" << iSR << "_tr20_v13_nofit" << ".png";
     myc[iSR]->SaveAs(canvasName.str().c_str());
   }


}
