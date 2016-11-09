// to draw the correlation plots of centrality variables
//
// Author Javier Martin
// 9 Nov 2016
//
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TPad.h"
#include "stdio.h"
#include "TSystem.h"
#include "../utils.h"

const char* trigcut = "(HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1 || HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1)";
const char* trigcap = "HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part*";
const char* evSelCut = "pBeamScrapingFilter && pPAprimaryVertexFilter && phfCoincFilter1 && pVertexFilterCutG";
const char* selCut ="" ;//= Form("%s && %s",trigcut,evSelCut);
const char* cap = "MBFullEvtSel";

void Get2DPlots(TTree* t_evt=0, TString v1="hiHF", TString v2="hiNpix", int xbin=200, double xmin=0, double xmax=4500, int ybin=200, double ymin=0, double ymax=10000, TCut cut="");

void plotCentVarsCorrel(const char* fname="root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pPb_EposMinBias_5TeV_8022_Forest/EPOS5TeV_GEN_SIM/crab_pPb_EposMinBias_5TeV_8022_Forest_corrCentrTable/161108_122954/0000/HiForestAOD_92.root")
{
  TCut trig(trigcut);
  
  TH1::SetDefaultSumw2();
  gStyle -> SetOptStat(0);
  SetHistTitleStyle(0.06,0.04);
  SetyjPadStyle();
  
  TFile *fin = TFile::Open(fname);
  TTree *t_evt = (TTree*) fin -> Get("hiEvtAnalyzer/HiTree");
  TTree *t_skim = (TTree*) fin -> Get("skimanalysis/HltTree");
  TTree *t_hlt = (TTree*) fin -> Get("hltanalysis/HltTree");
  t_evt -> AddFriend(t_hlt);
  t_evt -> AddFriend(t_skim);
  
  double hiNpixMax = 5000;
  double hiZDCMax = 80000;
  double hiHFMax = 500;
  double hiNTracksMax = 500;
  double hiHFhitMax = 3000;
  int nbin = 100;
  
  Get2DPlots(t_evt, "hiHFplus","hiHFminus",nbin,0,hiHFMax,nbin,0,hiHFMax,selCut);
  Get2DPlots(t_evt, "hiHFhitPlus","hiHFhitMinus",nbin,0,hiHFhitMax,nbin,0,hiHFhitMax,selCut);
  Get2DPlots(t_evt, "hiHF","hiNpix",nbin,0,hiHFMax,nbin,0,hiNpixMax,selCut);
  Get2DPlots(t_evt, "hiHF","hiZDC",nbin,0,hiHFMax,nbin,0,hiZDCMax,selCut);
  Get2DPlots(t_evt, "hiHF","hiNtracks",nbin,0,hiHFMax,nbin,0,hiNTracksMax,selCut);
  Get2DPlots(t_evt, "hiHFplusEta4","hiHFminusEta4",nbin,0,hiHFMax,nbin,0,hiHFMax,selCut);
  Get2DPlots(t_evt, "hiHFplusEta4","hiNtracks",nbin,0,hiHFMax,nbin,0,hiNTracksMax,selCut);
  Get2DPlots(t_evt, "hiHFplusEta4","hiZDC",nbin,0,hiHFMax,nbin,0,hiZDCMax,selCut);
  Get2DPlots(t_evt, "hiHFplusEta4","hiZDCplus",nbin,0,hiHFMax,nbin,0,hiZDCMax,selCut);
  Get2DPlots(t_evt, "hiHFplusEta4","hiZDCminus",nbin,0,hiHFMax,nbin,0,hiHFMax,selCut);
  Get2DPlots(t_evt, "hiNpix","hiZDC",nbin,0,hiNpixMax,nbin,0,hiZDCMax,selCut);
  Get2DPlots(t_evt, "hiZDC","hiNtracks",nbin,0,hiZDCMax,nbin,0,hiNTracksMax,selCut);
  
}

void Get2DPlots(TTree* t_evt, TString v1, TString v2, int xbin, double xmin, double xmax, int ybin, double ymin, double ymax, TCut cut)
{
  TCanvas *c1 = new TCanvas(Form("%sVS%s",v2.Data(), v1.Data()), Form("%sVS%s",v2.Data(), v1.Data()), 300,300);
  c1->cd();
  
  TH2D *h2D = new TH2D(Form("h2D_%s_%s_filter_%s",v1.Data(),v2.Data(),cap), Form("%s;%s;%s",cap, v1.Data(), v2.Data()), xbin, xmin, xmax, ybin, ymin, ymax);
  t_evt->Draw(Form("%s:%s>>+%s",v2.Data(), v1.Data(), h2D->GetName() ), selCut);
  h2D = (TH2D*)gDirectory->Get(h2D->GetName());

  h2D->Draw("colz");
  gPad->SetLogz();
  
  // save canvas
  string mainDIR = gSystem->ExpandPathName(gSystem->pwd());
  string saveDIR = mainDIR + "/2DcentralityDistributions/";
  void * dirp = gSystem->OpenDirectory(saveDIR.c_str());
  if (dirp) gSystem->FreeDirectory(dirp);
  else gSystem->mkdir(saveDIR.c_str(), kTRUE);
  c1->SaveAs(Form("2DcentralityDistributions/centVarsCorrel_%s_%s_%s.pdf", v1.Data(),v2.Data(),cap));
}

