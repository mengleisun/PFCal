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

void plotReso(){
TFile *inputFile = TFile::Open("CalibReso_vsE_IC2.root");
inputFile->cd("SR2");
gStyle->SetOptStat(0);
TGraphErrors *gr = (TGraphErrors*)gDirectory->Get("resoRecoFit2eta2pu1");
TCanvas *c = new TCanvas("reso","reso",1500,1000);
c->SetFillStyle(4000);
c->SetLogx();
//gPad->SetFillStyle(4000);
//gPad->SetFrameFillColor(0);
//gPad->SetFrameFillStyle(0);
//gPad->SetLogx();
TPad *newpad=new TPad("trans","trans",0,0,1,1);
newpad->SetFillStyle(4000);
newpad->SetFrameFillColor(0);
newpad->SetFrameFillStyle(0);
newpad->SetLogx();
newpad->Draw();
newpad->cd();

TH2F *dummy = new TH2F("reso","",1,7.6,2000,1,0,0.12);
dummy->GetXaxis()->SetLabelSize(0);
dummy->GetYaxis()->SetLabelSize(0);
dummy->Draw();


gr->SetMarkerSize(2);
gr->Draw("same p");

TF1 *fitFunc =new TF1("reso","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",7,600);

fitFunc->SetParameter(0,0.25);
fitFunc->SetParLimits(0,0,1);
fitFunc->SetParameter(1,0.017);
fitFunc->SetParLimits(1,0,1);
fitFunc->SetParameter(2,0);
fitFunc->SetParLimits(2,0,1);

fitFunc->SetLineColor(2);
fitFunc->SetLineWidth(4);
fitFunc->Draw("same");

c->SaveAs("transparent.eps");
}
