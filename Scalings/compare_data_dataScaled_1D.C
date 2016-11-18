#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "TStyle.h"
#include "../utils.h"

using namespace std;
// global vars definition
const char* trigcap("");
const char* evSelCut("");
const char* data2evSelCut("");
const char* evSelCutCap("");
const char* data2evSelCutCap("");
const char* data2Trigg("");
const char* saveLabel("");
const char* scaledD("");

TDatime* date = new TDatime();

void fit(TTree* ts,TTree* tref, TString var, int nBins, double xMin, double xMax, TCut dataCut = "pcollisionEventSelection", TCut dataRefCut = "pcollisionEventSelection");

void compare_data_dataScaled_1D(const char* fname_data="root://eoscms//eos/cms/store/group/phys_heavyions/kjung/ExpressForests/v1/Merged/HiForest_run285090_Express.root", const char* fname_data2="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_HLT_SingleTrack_JSonPPb_forestv84.root")
{
  gStyle->SetOptStat(0);
  
  TFile* f1 = TFile::Open(fname_data);
  TTree* t1 = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
  TTree* t1_skim = (TTree*) f1 -> Get("skimanalysis/HltTree");
  TTree* t1_hlt = (TTree*) f1 -> Get("hltanalysis/HltTree");
  t1->AddFriend(t1_skim);
  t1->AddFriend(t1_hlt);
  TFile* f2 = TFile::Open(fname_data2);
  TTree* t2 = (TTree*) f2 -> Get("hiEvtAnalyzer/HiTree");
  TTree* t2_skim = (TTree*) f2 -> Get("skimanalysis/HltTree");
  TTree* t2_hlt = (TTree*) f2 -> Get("hltanalysis/HltTree");
  t2->AddFriend(t2_skim);
  t2->AddFriend(t2_hlt);
  
  //total
  double hiHFMax = 300;
  double hiHFSSMax = 200;
  double hiHFSSTruncMax = 80;
  double hiBinMax = 200;
  double hiHFhitMax = 8000;
  double hiETMax = 50;
  double hiEEMax = 50;
  double hiEBMax = 100;
  double hiNpixMax = 2000;
  double hiNpixelTracksMax = 500;
  double hiNtracksMax = 400;
  double hiNtracksCutEtaMax = 350;
  double hiNtracksCutMax = 150;
  double hiZDCMax = 40000;
  
  int nBin = 50;
  //  const char* trigcut = "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1)";
  const char* trigcut = "HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1";
  //  trigcap = "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part*";
  trigcap = "HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1";
  evSelCut = "pBeamScrapingFilter && pPAprimaryVertexFilter && phfCoincFilter1 && pVertexFilterCutG";
  const string cap = "pAExpress_5TeV_run285090";
  evSelCutCap = "BS+PV+HFC+PVG";
  TCut selCut = Form("%s && %s",trigcut,evSelCut);
  cout << "Aplying the following event selection to 2015 data: " << selCut.GetTitle() << endl;
  
  data2evSelCut = "pPAcollisionEventSelectionPA && run<211256";
  data2Trigg = "HLT_PAZeroBiasPixel_SingleTrack_v1";
  TCut data2Cut = Form("%s && %s",data2Trigg,data2evSelCut);
  data2evSelCutCap = "pPAcollisionEventSelectionPA ; run<211256";
  cout << "Aplying the following event selection to 2013 data: " << data2Cut.GetTitle() << endl;
  
  saveLabel = "MBFullSel";
  scaledD = "2013 data";
  
  fit(t2, t1, "hiHF",nBin,0,hiHFMax,data2Cut,selCut);
  fit(t2, t1, "hiHFplus",nBin,0,hiHFSSMax,data2Cut,selCut);
  fit(t2, t1, "hiHFminus",nBin,0,hiHFSSMax,data2Cut,selCut);
  fit(t2, t1, "hiHFplusEta4",nBin,0,hiHFSSTruncMax,data2Cut,selCut);
  fit(t2, t1, "hiHFminusEta4",nBin,0,hiHFSSTruncMax,data2Cut,selCut);
//  fit(t2, t1, "hiBin",nBin,0,hiBinMax,data2Cut,selCut);
  fit(t2, t1, "hiHFhit",nBin,0,hiHFhitMax,data2Cut,selCut);
//  fit(t2, t1, "hiET",nBin,0,hiETMax,data2Cut,selCut);
//  fit(t2, t1, "hiEB",nBin,0,hiEBMax,data2Cut,selCut);
//  fit(t2, t1, "hiEE",nBin,0,hiEEMax,data2Cut,selCut);
  fit(t2, t1, "hiNpix",nBin,0,hiNpixMax,data2Cut,selCut);
  fit(t2, t1, "hiNtracks",nBin,0,hiNtracksMax,data2Cut,selCut);
//  fit(t2, t1, "hiNtracksPtCut",nBin,0,hiNtracksCutMax,data2Cut,selCut);
//  fit(t2, t1, "hiNtracksEtaCut",nBin,0,hiNtracksCutEtaMax,data2Cut,selCut);
//  fit(t2, t1, "hiNtracksEtaPtCut",nBin,0,hiNtracksCutMax,data2Cut,selCut);
//  //  fit(t2, t1, "hiNpixelTracks",nBin,0,hiNpixelTracksMax,data2Cut,selCut);
//  fit(t2, t1, "hiZDC",nBin,0,hiZDCMax,data2Cut,selCut);
//  fit(t2, t1, "hiZDCplus",nBin,0,hiZDCMax,data2Cut,selCut);
//  fit(t2, t1, "hiZDCminus",nBin,0,hiZDCMax,data2Cut,selCut);
  
}

double chi2(TH1* h1, TH1* h2, double rangeMin){
  
  
  double c = 0;
  int ndof = 0;
  
  for(int i = h1->GetXaxis()->FindBin(rangeMin); i < h1->GetNbinsX()+1; ++i){
    double y1 = h1->GetBinContent(i);
    double y2 = h2->GetBinContent(i);
    double e1 = h1->GetBinError(i);
    double e2 = h2->GetBinError(i);
    
    double dy = y1-y2;
    double de2 = e1*e1+e2*e2;
    if(de2 > 0) {
      c += dy*dy/de2;
      ndof += 1;
    }
  }
  //return c;
  return c/ndof;
}


TH1* scale(int index, const char* var, TTree* nt, double s = 1, int nBins=100, double xMin=0, double xMax=1000, TCut cut = "", double rangeMin = 0)
{
  TH1D* h = new TH1D(Form("h_%s_%d",var,index),Form(";%s [GeV];event fraction",var),nBins,xMin,xMax);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  
  nt->Draw(Form("%f*%s>>%s",s,var,h->GetName()),cut);
  
  double rangeIntegral = h->Integral(h->GetXaxis()->FindBin(rangeMin),h->GetNbinsX()+1);
  if(rangeIntegral > 0){
    h->Scale(1./rangeIntegral);
  }
  
  return h;
  
}


void fit(TTree* ts,TTree* tref, TString var, int nBins, double xMin, double xMax, TCut dataCut, TCut dataRefCut)
{
  TH1::SetDefaultSumw2();
  
  static double rangeMin = 40;
  double variation = 0.2;
  double scaleMin = 0.95;
  if (!strcmp(var.Data(),"hiNtracks"))
  {
    scaleMin = 0.96;
    variation = 0.06;
    rangeMin = 50;
  }
  if (!strcmp(var.Data(),"hiHFminusEta4") || !strcmp(var.Data(),"hiHFplusEta4"))
  {
    rangeMin = 10;
    variation = 0.1;
  }
  if (!strcmp(var.Data(),"hiHFminus"))
  {
    scaleMin = 0.98;
    variation = 0.14;
    rangeMin = 30;
  }
  if (!strcmp(var.Data(),"hiHFplus")) rangeMin = 30;
  if (!strcmp(var.Data(),"hiHFhit")) rangeMin = 1500;
  if (!strcmp(var.Data(),"hiNpix"))
  {
    rangeMin = 400;
    scaleMin = 0.7;
  }
  
  TH1* href = scale(0,var.Data(),tref,1., nBins, xMin,xMax,dataRefCut,rangeMin);
  
  cout<<"REFERENCE DATA ENTRIES : "<<tref->GetEntries(dataRefCut)<<endl;
  cout<<"SCALED DATA ENTRIES : "<<ts->GetEntries(dataCut)<<endl;
  
  static const int N = 20;
  TH1* h[N];
  double chi2s[N];
  double s[N];
  int ibest = 0;
  TGraphErrors* g = new TGraphErrors(N);
    
  for(int i = 0; i < N; ++i){
    s[i] = scaleMin+(variation/N)*i;
    h[i] = scale(i+1,var.Data(),ts,s[i],nBins, xMin,xMax,dataCut,rangeMin);
    h[i]->SetLineColor(2);
    h[i]->SetMarkerColor(2);
    h[i]->SetMarkerStyle(20);
    chi2s[i] = chi2(h[i], href,rangeMin);
    g->SetPoint(i,s[i],chi2s[i]);
    if(chi2s[i]< chi2s[ibest]) ibest = i;
    cout<<"Scale factor = " << s[i] << " ; chi2 = "<<chi2s[i]<<endl;
    
  }
  double chiMax = chi2s[0]*(1.2);
  double chiMin = chi2s[ibest]*(0.8);
  
  double eff = href->Integral()/(href->Integral(href->GetXaxis()->FindBin(rangeMin)+1,href->GetNbinsX()+1) + h[ibest]->Integral(0,h[ibest]->GetXaxis()->FindBin(rangeMin)));
  //double eff = href->Integral()/h[ibest]->Integral();
  
  cout<<"the best scale is : "<<s[ibest]<<" i = "<<ibest<<endl;
  cout<<"relative efficiency is : "<<eff<<endl;
  
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/scalingFactors/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  
  ofstream outputf(Form("scalingFactors/centEff_fitResults_%s_%s.txt",var.Data(),saveLabel), ios::app);
  outputf << setprecision(4) <<eff*100 << "\t" << chi2s[ibest] << endl;
  outputf.close();
  
  TCanvas* c1 = new TCanvas(Form("c1_%s",var.Data()),Form("c1_%s",var.Data()),600,600);
  g->GetXaxis()->SetTitle("scale factor");
  g->GetYaxis()->SetTitle("reduced #chi^2 (normalized by ndf)");
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();
  g->SetMarkerStyle(20);
  g->GetYaxis()->SetRangeUser(chiMin,chiMax);
  g->Draw("Ap");
  
  
  c1->Print(Form("scalingFactors/figure_chi2_%s_%s_%d.pdf",var.Data(),saveLabel,date->GetDate()));
  
  TCanvas* c2 = new TCanvas(Form("c2_%s",var.Data()),Form("c2_%s",var.Data()),600,600);
  c2->SetLogy();
  
  h[ibest]->SetMaximum(10.*h[ibest]->GetMaximum());
  h[ibest]->Draw();
  href->SetMarkerStyle(20);
  href->SetMarkerSize(0.8);
  href->Draw("same");
  TLine* vertical = new TLine(rangeMin, h[ibest]->GetMinimum(), rangeMin, h[ibest]->GetMaximum());
  vertical->SetLineStyle(2);
  vertical->Draw("same");
  
  TLegend * leg1 = new TLegend(0.4,0.6,0.90,0.85);
  leg1->SetFillStyle(1);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  
  leg1->SetTextSize(0.025);
  leg1->AddEntry(href,"2016 data","p");
  leg1->AddEntry(href,trigcap,"");
  leg1->AddEntry(href,evSelCutCap,"");
  
  leg1->AddEntry(h[ibest],Form("%s based fit",scaledD),"lp");
  leg1->AddEntry(h[ibest],Form("%s scaling : %f",scaledD,s[ibest]),"");
  leg1->AddEntry(h[ibest],Form("Rel. Eff. : %f",eff),"");
  leg1->AddEntry(h[ibest],Form("#chi^{2}/ndf : %f",chi2s[ibest]),"");
  leg1->Draw();
  
  c2->Print(Form("scalingFactors/figure_fit_%s_%s_%d.pdf",var.Data(),saveLabel,date->GetDate()));
  
  TCanvas* c3 = new TCanvas(Form("c3_%s",var.Data()),Form("c3_%s",var.Data()),600,600);
  c3->SetLogx();
  TH1D* hdiff = (TH1D*)href->Clone("hdiff");
  hdiff->SetTitle(Form(";%s [GeV];data - fit",var.Data()));
  hdiff->Add(h[ibest],-1);
  hdiff->Draw();
  TLine* horizontal = new TLine(xMin,0,xMax,0);
  horizontal->SetLineStyle(2);
  horizontal->Draw("same");
  
  
  TLegend * leg2 = new TLegend(0.27,0.15,0.82,0.44);
  leg2->SetFillStyle(1);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  
  leg2->SetTextSize(0.025);
  leg2->AddEntry(hdiff,"data - fit","p");
  leg2->AddEntry(href,trigcap,"");
  leg2->AddEntry(href,evSelCutCap,"");
  leg2->AddEntry(href,data2Trigg,"");
  leg2->AddEntry(href,data2evSelCutCap,"");
  
  leg2->Draw();
  
  c3->Print(Form("scalingFactors/figure_diff_%s_%s_%d.pdf",var.Data(),saveLabel,date->GetDate()));
  
}

