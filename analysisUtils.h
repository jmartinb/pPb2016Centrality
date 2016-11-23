// Javier Martin Blanco
// 23 Nov 2016

#ifndef analysisUtils
#define analysisUtils

#include "TTree.h"
#include "TLorentzVector.h"
#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <string.h>

using namespace std;

const double HFin = 2.87;
const double HFout = 5.02;
const float eThr = 3;
const int nTower= 1;

bool isDSEvent(TTree* t, int jentry)
{
  // This function tells if a MC event is Double Sided or not
  // It checks if there is at least one pion+-, pion0, K+-, K0L, proton, neutron,electron, muon or photon in
  // both sides (3<|eta|<5) with an energy > eThr
  // To use this function you need to: TTree *t_gen = (TTree*) fin -> Get("HiGenParticleAna/hi");
  // and t-> AddFriend(t_gen);
  
  Int_t mult;
  std::vector<float>* pt=0;
  std::vector<float>* eta=0;
  std::vector<float>* phi=0;
  std::vector<int>*   pdg=0;
  t->SetBranchAddress("mult",&mult);
  t->SetBranchAddress("pt",&pt);
  t->SetBranchAddress("eta",&eta);
  t->SetBranchAddress("phi",&phi);
  t->SetBranchAddress("pdg",&pdg);
  
  t->GetEntry(jentry);
  
  int nPosTower = 0;
  int nNegTower = 0;
  for(int i=0; i<mult;i++){
    float m = -1.;
    if(fabs(pdg->at(i)) == 111) m = 0.13498; //pion 0
    else if(fabs(pdg->at(i)) == 211) m = 0.13957; //pion +-
    else if(fabs(pdg->at(i)) == 321) m = 0.49368; //kaon +-
    else if(fabs(pdg->at(i)) == 130) m = 0.49761; //kaon 0 long
    else if(fabs(pdg->at(i)) == 2212) m = 0.93827; //proton
    else if(fabs(pdg->at(i)) == 2112) m = 0.93957; //neutron
    else if(fabs(pdg->at(i)) == 22) m = 0.; //photon
    else if(fabs(pdg->at(i)) == 11) m = 0.000511; //electron
    else if(fabs(pdg->at(i)) == 13) m = 0.10566; //muon
    else continue;
    TLorentzVector * v = new TLorentzVector();
    v->SetPtEtaPhiM(pt->at(i),eta->at(i),phi->at(i),m);
    double e = v->Energy();
    delete v;
    if(e<eThr) continue;
    if( (eta->at(i)>HFin ) && ( (eta->at(i))<HFout ) ) nPosTower++;
    else if( (eta->at(i)>-HFout ) && ( (eta->at(i))<-HFin ) ) nNegTower++;
    else continue;
  }
  
  if(nPosTower>=nTower && nNegTower>=nTower) return true;
  
  return false;
}

#endif
