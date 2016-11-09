// plot centrality variables
// Author : Javier Martin
// 8 Nov 2016

//basic c++ header, string ...
#include <iostream>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>
//tree, hist, vector ...
#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>
#include <math.h>
#include <TH1.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TCut.h>
//canvas, legend, latex ... //cosmetic
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TAxis.h>

//private setup
#include "../utils.h"

// global vars definition
const char* trigcap("");
const char* evSelCut("");

void plotVar(TTree* t1=0 ,TString var="hiHF", int nBins=10, double xMin=0, double xMax=10, TCut cut1="(1)", const string cap = "");
void plotCentVarsDistr(const char* fname = "/afs/cern.ch/user/k/kjung/public/forestTests/HiForestAOD_run284755_lumi414.root")
{
  TFile* f1 = TFile::Open(fname);
  TTree* t1 = (TTree*) f1 -> Get("hiEvtAnalyzer/HiTree");
  TTree* t1_skim = (TTree*) f1 -> Get("skimanalysis/HltTree");
  TTree* t1_hlt = (TTree*) f1 -> Get("hltanalysis/HltTree");
  t1->AddFriend(t1_skim);
  t1->AddFriend(t1_hlt);
  
  //total
  double hiHFMax = 9000;
  double hiBinMax = 200;
  double hiHFhitMax = 250000;
  double hiETMax = 1800;
  double hiEEMax = 3500;
  double hiEBMax = 4000;
  double hiNpixMax = 50000;
  double hiNpixelTracksMax = 5000;
  double hiNtracksMax = 3500;
  double hiZDCMax = 1;
  
  int nBin = 50;
  const char* trigcut = "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1)";
  trigcap = "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part*";
  evSelCut = "pBeamScrapingFilter && pPAprimaryVertexFilter && phfCoincFilter1 && pVertexFilterCutG";
  const char* selCut ="" ;//= Form("%s && %s",trigcut,evSelCut);
  const string cap = "pAMinimumBias_runX";
  plotVar(t1, "hiHF",nBin,0,hiHFMax,selCut, cap);
  plotVar(t1, "hiHFplus",nBin,0,hiHFMax,selCut, cap);
  plotVar(t1, "hiHFminus",nBin,0,hiHFMax,selCut, cap);
  plotVar(t1, "hiHFplusEta4",nBin,0,hiHFMax,selCut, cap);
  plotVar(t1, "hiHFminusEta4",nBin,0,hiHFMax,selCut, cap);
  plotVar(t1, "hiBin",nBin,0,hiBinMax,selCut, cap);
  plotVar(t1, "hiHFhit",nBin,0,hiHFhitMax,selCut, cap);
  plotVar(t1, "hiET",nBin,0,hiETMax,selCut, cap);
  plotVar(t1, "hiEB",nBin,0,hiEBMax,selCut, cap);
  plotVar(t1, "hiEE",nBin,0,hiEEMax,selCut, cap);
  plotVar(t1, "hiNpix",nBin,0,hiNpixMax,selCut, cap);
  plotVar(t1, "hiNtracks",nBin,0,hiNtracksMax,selCut, cap);
  plotVar(t1, "hiNtracksPtCut",nBin,0,hiNtracksMax,selCut, cap);
  plotVar(t1, "hiNtracksEtaCut",nBin,0,hiNtracksMax,selCut, cap);
  plotVar(t1, "hiNtracksEtaPtCut",nBin,0,hiNtracksMax,selCut, cap);
  //  plotVar(t1, "hiNpixelTracks",nBin,0,hiNpixelTracksMax,selCut, cap);
  plotVar(t1, "hiZDC",nBin,0,hiZDCMax,selCut, cap);
  plotVar(t1, "hiZDCplus",nBin,0,hiZDCMax,selCut, cap);
  plotVar(t1, "hiZDCminus",nBin,0,hiZDCMax,selCut, cap);
  
} // main function

void plotVar(TTree* t1, TString var, int nBins, double xMin, double xMax, TCut cut1, const string cap)
{
  SetHistTitleStyle();
  SetyjPadStyle();
  gStyle->SetOptStat(0);
  
  static int j = 0;
  TCanvas* c=  new TCanvas(Form("c_%s_%d",var.Data(),j),"", 700,600);
  if (strcmp(var,"hiBin")) gPad->SetLogy();
  
  TH1D* h1 = new TH1D(Form("h1_%s_%d",var.Data(),j), Form(";%s;",var.Data()), nBins,xMin,xMax);
  h1->Sumw2();
  t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), cut1);
  //  h1->Scale( 1. / t1->GetEntries(cut1));
  SetHistColor(h1,2);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.8);
  h1->SetNdivisions(505);
  TLegend* l1 = new TLegend(0.5,0.6,0.90,0.8);
  legStyle(l1);
  l1->AddEntry(h1, "Data", "p");
  
  //  double range = cleverRange(h1,1.5,1.e-4);
  h1->DrawCopy("L");
  l1->Draw("same");
  
  drawText(cap.data(),0.2,0.2+0.7);
  //drawText(cut1.GetTitle(),0.2,0.2);
  drawText(trigcap,0.2,0.2+0.06);
  drawText(evSelCut,0.2,0.2);
  
  // save canvas
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/centralityDistributions/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  c->SaveAs(Form("centralityDistributions/distribution_%s_%s.pdf",var.Data(),cap.data()));
  j++;
}
