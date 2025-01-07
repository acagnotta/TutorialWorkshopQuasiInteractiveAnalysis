//#ifndef POST_H
//#define POST_H

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <map>

#include "TDavixFile.h"

using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using rvec_f = const RVec<float> &;
using rvec_i = const RVec<int> &;
using rvec_b = const RVec<bool> &;
using rvec_rvec_i = const RVec<RVec<int>> &;

const float TopRes_trs=  0.5411276;
const float TopMix_trs=  0.7584613561630249;
const float TopMer_trs=  0.8;//0.94; 0.8 for 2022 correspond to fpr 6% on ttbar
const float dR=  0.8;

const float btagDeepB_mediumWP_2018   = 0.2783;
const float btagPNet_mediumWP_2022    = 0.245 ;
const float btagPNet_mediumWP_2022EE  = 0.2605;

const float btagDeepB_looseWP_2018   = 0.0490;
const float btagPNet_looseWP_2022    = 0.047 ;
const float btagPNet_looseWP_2022EE  = 0.0499;

//  Top Resolved threshold {'fpr 10': 0.1334565, 'fpr 5': 0.24193972, 'fpr 1': 0.5411276, 'fpr 01': 0.77197933}
//  Top Mixed threshold {'fpr 10': 0.13067308068275452, 'fpr 5': 0.2957885265350342, 'fpr 1': 0.7584613561630249, 'fpr 01': 0.9129540324211121}

// ########################################################
bool isMC(int SampleFlag){
  if (SampleFlag == 0) return false;
  else return true;
}

// ############### from skimtree_utils

float deltaPhi (float phi1, float phi2){
  float dphi = (phi1-phi2);
  while(dphi >  M_PI) dphi -= 2*M_PI;
  while(dphi < -M_PI) dphi += 2*M_PI;
  return dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2){
  return hypot(eta1 - eta2, deltaPhi(phi1, phi2)); 
}


// ########################################################
// ############## Data2018 ################################
// ########################################################

bool hemveto(rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Electron_eta, rvec_f Electron_phi){
  float hemvetoetaup = -3.05;
  float hemvetoetadown = -1.35;
  float hemvetophiup = -1.62;
  float hemvetophidown = -0.82;
  bool passesMETHEMVeto = true;

  for(size_t i = 0; i < Jet_eta.size(); i++){
    if(Jet_eta[i]>hemvetoetaup && Jet_eta[i]<hemvetoetadown && Jet_phi[i]>hemvetophiup && Jet_phi[i]<hemvetophidown){
      passesMETHEMVeto = false;
    }
  }
  for(size_t i = 0; i < Electron_eta.size(); i++){
    if(Electron_eta[i]>hemvetoetaup && Electron_eta[i]<hemvetoetadown && Electron_phi[i]>hemvetophiup && Electron_phi[i]<hemvetophidown){
      passesMETHEMVeto = false;
    }
  }
  return passesMETHEMVeto; 
}

int w_nominalhemveto(float w_nominal, bool HEMVeto){
  int w_nominal_update = w_nominal;
  if (HEMVeto == false){
    w_nominal_update = w_nominal * 0.354;
    // w_nominal_update = w_nominal * 0.932;
  }
  return w_nominal_update;
}


// ########################################################
// ########## PRESELECTION ################################
// ########################################################
int nVetoElectron(rvec_f Electron_pt, rvec_f Electron_cutBased, rvec_f Electron_eta)
{
  int n=0;
  for(int i = 0; i<Electron_pt.size(); i++)
  {
    if(Electron_cutBased[i]>=1 && Electron_pt[i] > 30 && abs(Electron_eta[i])<2.5) n+=1;
  }
  return n;
}

int nVetoMuon(rvec_f Muon_pt, rvec_f Muon_eta, rvec_f Muon_looseId)
{
  int n=0;
  for(int i = 0; i<Muon_pt.size(); i++)
  {
    if(Muon_looseId[i]==1 && Muon_pt[i] > 30 && abs(Muon_eta[i])<2.4) n+=1;
  }
  return n;
}

RVec<int> GetGoodJet(rvec_f Jet_pt, rvec_f Jet_eta, rvec_i Jet_jetId)
{
  RVec<int> ids;
  for(int i = 0; i<Jet_pt.size(); i++)
  {
      if (Jet_pt[i]>30 && abs(Jet_eta[i])<2.4 && Jet_jetId[i]>1)
      {
        ids.emplace_back(i);
      }
  }
  return ids;
}

RVec<int> GetGoodFatJet(rvec_f Jet_pt, rvec_f Jet_eta, rvec_i Jet_jetId, rvec_f FatJet_particleNetWithMass_TvsQCD)
{
  RVec<int> ids;
  for(int i = 0; i<Jet_pt.size(); i++)
  {
    // taglio in eta portato da 2.7 a 2.4 -> per definizione forward jets
      if (Jet_pt[i]>150 && abs(Jet_eta[i])<2.4 && Jet_jetId[i]>1 && FatJet_particleNetWithMass_TvsQCD[i]>0.)
      {
        ids.emplace_back(i);
      }
  }
  return ids;
}

int nGoodJet(rvec_i GoodJet_idx)
{
  return GoodJet_idx.size();
}


// ########################################################
// ######## Event variables (deltaphi, etajet) ############
// ########################################################

float min_DeltaPhi(float MET_phi, rvec_f Jet_phi, rvec_i GoodJet_idx)
{
  float min_dphi = 4;
  for(int i = 0; i < GoodJet_idx.size(); i++)
  {
    float dphi = deltaPhi(MET_phi, Jet_phi[GoodJet_idx[i]]);
    if (dphi < min_dphi) min_dphi = dphi;
  }
  return min_dphi;
}


Int_t nForwardJet(rvec_f Jet_pt, rvec_f Jet_jetId, rvec_f Jet_eta)
{
  int nfwdjet = 0;
  for(int i = 0; i < Jet_pt.size(); i++)
  {
    if (Jet_pt[i]>30 && Jet_jetId[i] && abs(Jet_eta[i])>2.4)
    {
      nfwdjet += 1;
    }
  }
  return nfwdjet;
}
// ########################################################
// ############## TopSelection ############################
// ########################################################

Int_t select_bestTop(rvec_f FatJet_particleNet_TvsQCD, rvec_i GoodFatJet_idx){
  RVec<float> scores;
  int idx;
  //topselect = GoodTopMer_idx;
  for(int i = 0; i < GoodFatJet_idx.size(); i++)
  {
      scores.emplace_back(FatJet_particleNet_TvsQCD[GoodFatJet_idx[i]]);
  }
  idx = ArgMax(scores);
  return idx;
}