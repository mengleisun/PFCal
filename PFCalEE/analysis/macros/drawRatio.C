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


void drawRatio()
{
 TCanvas *myc;
 myc= new TCanvas("resolution_eta29","Electron Resolution",0,0,600,600);

 myc->cd();

 unsigned mips[3] = {4,8,20};
 unsigned version[2]={12,13};

 double slope[2][3];
 double offset[2][3];
 double noiseterm[2][3];

 const unsigned eta = 29;
 const double theta = 0.110;
 double EtToEfactor = sin(theta);
 unsigned etaIndex = 6;
 unsigned PUIndex=1;
 double Et[10]={2,5,10,20,40,60,80,100,150,200};
 for(unsigned i(0); i<10; i++){
    Et[i] = Et[i]/EtToEfactor;
 }
 double ResoReco[2][3][10];
 double ResoRecoErr[2][3][10];

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
    graphName << "resoRecoFit2eta" << eta <<"pu" << PUIndex; 
    reso= (TGraphErrors *)Reso_collection->Get(graphName.str().c_str());
    unsigned yindex = 0; 
    for(unsigned i(0); i<10; i++){
       if(abs(Et[i] - reso->GetX()[yindex]) < 2){
          ResoReco[iversion][ithres][i]=reso->GetY()[yindex];
          ResoRecoErr[iversion][ithres][i]=reso->GetErrorY(yindex);          
          yindex +=1;
        } else {
          ResoReco[iversion][ithres][i]=0;
          ResoRecoErr[iversion][ithres][i]=0;
        } 
       std::cout << Et[i] <<" " << ResoReco[iversion][ithres][i] << " " << ResoRecoErr[iversion][ithres][i] << std::endl;
    }  

    TF1 *resoModel=new TF1("resomodel","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0,Et[9]);
    resoModel->SetParameter(0,0.2);
    resoModel->SetParLimits(0,0,1);
    resoModel->SetParameter(1,0);
    resoModel->SetParLimits(1,0,1);
  //  resoModel->SetParLimits(2,0,0);
    resoModel->FixParameter(2,0);
    resoModel->SetLineColor(iversion+1);

   reso->GetYaxis()->SetRangeUser(0.00,0.10);
   reso->GetYaxis()->SetLabelSize(0.03);
   reso->GetXaxis()->SetLabelSize(0.03);
   reso->SetLineColor(iversion+1);
   reso->SetMarkerStyle(ithres+20);
   reso->SetMarkerSize(1.5);
   reso->SetMarkerColor(iversion+1);
//   reso->Fit(resoModel,"MER+");

//    slope[iversion][ithres]=resoModel->GetParameter(0);
//    offset[iversion][ithres]=resoModel->GetParameter(1);
//    noiseterm[iversion][ithres]=resoModel->GetParameter(2);
 
     leg = new TLegend(0.6, 0.8-ithres*0.06-iversion*0.18, 0.9, 0.86-ithres*0.06-iversion*0.18);
//     leg->AddEntry(reso,Form("%dmm,%dMips #frac{#sigma}{E} #propto #frac{%3.3f}{#sqrt{E}} + #frac{%3.3f}{E} + %3.3f",version[iversion]==12?2:4,mips[ithres]/4,slope[iversion][ithres],noiseterm[iversion][ithres], offset[iversion][ithres]),"p");
     leg->AddEntry(reso,Form("%dmm,%dMips",version[iversion]==12?2:4,mips[ithres]/4),"p");
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

double ratio[3][10];
double ratioerr[3][10];
for(unsigned ithres(0); ithres < 3; ithres ++){
    for(unsigned i(0); i< 10; i++){
        if(ResoReco[0][ithres][i] > 0 && ResoReco[1][ithres][i] > 0){
           ratio[ithres][i] = ResoReco[1][ithres][i]/ResoReco[0][ithres][i];
           ratioerr[ithres][i] = ratio[ithres][i]*sqrt(ResoRecoErr[1][ithres][i]*ResoRecoErr[1][ithres][i]/(ResoReco[1][ithres][i]*ResoReco[1][ithres][i])+
             ResoRecoErr[0][ithres][i]*ResoRecoErr[0][ithres][i]/(ResoReco[0][ithres][i]*ResoReco[0][ithres][i]) );
        } else {
          ratio[ithres][i] = 0;
        }
    }
}  


TCanvas *myr;
 myr= new TCanvas("ratio_eta29","Electron Resolution",0,0,600,200);

 myr->cd();
TGraphErrors *ratioplot[3];
for(unsigned ithres(0); ithres <3; ithres++){
    ratioplot[ithres] = new TGraphErrors();
    for(unsigned i(0); i<10; i++){
       if(ratio[0][i] > 0.1 && ratio[0][i] < 2)
          { Int_t np=ratioplot[ithres]->GetN();
            ratioplot[ithres]->SetPoint(np,Et[i],ratio[ithres][i]);
            ratioplot[ithres]->SetPointError(np,0,ratioerr[ithres][i]);
          }
    }
     ratioplot[ithres]->SetLineColor(ithres+1);
     ratioplot[ithres]->SetMarkerColor(ithres+1);
     ratioplot[ithres]->SetMarkerStyle(21);
     TLegend *leg = new TLegend(0.5, 0.8-ithres*0.06, 0.9, 0.86-ithres*0.06);
     leg->AddEntry(ratioplot[ithres],Form("4mm / 2mm resolution, %dMips",mips[ithres]/4),"p");
     leg->SetTextSize(0.03);
     leg->SetBorderSize(0);
     leg->SetFillColor(0);

    ratioplot[ithres]->GetYaxis()->SetRangeUser(0.6,1.4);
    if(ithres==0){ratioplot[ithres]->Draw(); leg->Draw();}
    else {ratioplot[ithres]->Draw("same"); leg->Draw("same");}
}
 
        
     myc->SaveAs("resolution_eta29_pu140.eps");
     myr->SaveAs("ratio_eta29_pu140.eps");
     myc->SaveAs("resolution_eta29_pu140.png");
     myr->SaveAs("ratio_eta29_pu140.png");

}
