////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef _PHENOANALYZER_H_
#define _PHENOANALYZER_H_

#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TEnv.h"
#include "TRandom.h"

#include <iostream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"

#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TApplication.h"
#include "DelphesFunctions.h"
#include "TDirectory.h"
#include "TFile.h"
#include <fstream>

using namespace std;

class PhenoAnalysis {
public :
   PhenoAnalysis(TChain&, TFile*, TDirectory* dir[], int nDir);
   ~PhenoAnalysis();
   void crateHistoMasps (int);
   bool overlapingObjects(double, double, double, double, double);
   double calculateE(double, double, double);
   double normalizedDphi(double);

   // For Jets
   std::map<unsigned int, TH1*> _hmap_Nevents;
   std::map<unsigned int, TH1*> _hmap_lead_jet_pT;
   std::map<unsigned int, TH1*> _hmap_lead_jet_eta;
   std::map<unsigned int, TH1*> _hmap_lead_jet_phi;
   std::map<unsigned int, TH1*> _hmap_n_jets;
   std::map<unsigned int, TH1*> _hmap_slead_jet_pT;
   std::map<unsigned int, TH1*> _hmap_slead_jet_eta;
   std::map<unsigned int, TH1*> _hmap_slead_jet_phi;
   // For Taus
   std::map<unsigned int, TH1*> _hmap_tau1_pT;
   std::map<unsigned int, TH1*> _hmap_tau1_eta;
   std::map<unsigned int, TH1*> _hmap_tau1_phi;
   std::map<unsigned int, TH1*> _hmap_tau2_pT;
   std::map<unsigned int, TH1*> _hmap_tau2_eta;
   std::map<unsigned int, TH1*> _hmap_tau2_phi;
   std::map<unsigned int, TH1*> _hmap_n_tau;
   // Topology
   std::map<unsigned int, TH1*> _hmap_ht;
   std::map<unsigned int, TH1*> _hmap_st;
   std::map<unsigned int, TH1*> _hmap_dijet_mass;
   std::map<unsigned int, TH1*> _hmap_dijet_deltaEta;

   //2D plots
   std::map<unsigned int, TH2*> _hmap_jet_met_metDphi;
   std::map<unsigned int, TH2*> _hmap_Dphi_tau_met_with_met;
   std::map<unsigned int, TH2*> _hmap_Dphi_tau_jet_with_jPt;
   std::map<unsigned int, TH2*> _hmap_MT_with_MET;
   std::map<unsigned int, TH2*> _hmap_MT_with_Dphi_tau_met;
   std::map<unsigned int, TH2*> _hmap_MT_with_Dphi_jet_met;
   std::map<unsigned int, TH2*> _hmap_Meff_with_MET;
};

#endif
