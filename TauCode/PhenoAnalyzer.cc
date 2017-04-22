/*
@file PhenoAnalyzer.cc
@author Andres Florez
@author Carlos Miguel Patino
@date April 2, 2017

Code used to perform phenomenological analysis of Heavy Neutrinos in the tau channel
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include "ROOTFunctions.h"
#include "DelphesFunctions.h"
#include "PhenoAnalyzer.h"
#include <time.h>

int main(int argc, char *argv[]) {

  //TApplication app("App",&argc, argv);
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile * HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 11;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("Taus_pT_min");
  theDirectory[2]  = HistoOutputFile->mkdir("Taus_eta_min");
  theDirectory[3]  = HistoOutputFile->mkdir("Taus_mass_min");
  theDirectory[4]  = HistoOutputFile->mkdir("MET");
  theDirectory[5]  = HistoOutputFile->mkdir("N_bjets");
  theDirectory[6]  = HistoOutputFile->mkdir("Jets_pT_min");
  theDirectory[7]  = HistoOutputFile->mkdir("Transmass_min");
  theDirectory[8]  = HistoOutputFile->mkdir("VBF_jets_opposite_hemispheres");
  theDirectory[9]  = HistoOutputFile->mkdir("VBF_jets_delta");
  theDirectory[10]  = HistoOutputFile->mkdir("VBF_diJetMass");
  PhenoAnalysis BSM_analysis(chain, HistoOutputFile, theDirectory, nDir);

}

using namespace std;
PhenoAnalysis::PhenoAnalysis(TChain& chain, TFile* theFile, TDirectory *cdDir[], int nDir)
{
  ifstream inFile;
  inFile.open ("config.in", ios::in);

  if (!inFile)
    {
      cerr << "ERROR: Can't open input file: " << endl;
      exit (1);
    }

  string inputType = "";

  //This set of lines are used to open and read the "config.in" file.
  ///////////////////////////////////////////////////////////////////////
  TEnv *params = new TEnv ("config_file");
  params->ReadFile ("config.in", kEnvChange);

  double tau_pt_min      = params->GetValue ("tau_pt_min", 10.);
  double b_jet_pt_min     = params->GetValue ("b_jet_pt_min", 30.0);
  double DR_jet_lep_max  = params->GetValue ("DR_jet_lep_max", 0.3);
  double jet_min_pt       = params->GetValue ("jet_min_pt", 15);
  double VBF_jetPt_min = params->GetValue("VBF_jetPt_min", 40.0);
  double tau_pt_cut       = params->GetValue("tau_pt_cut", 20.);
  double tau_eta_cut      = params->GetValue ("tau_eta_cut", 2.1);
  double deltaEta_diJet_cut = params->GetValue("deltaEta_diJet_cut", 3.8);
  double diJetmass_cut    = params->GetValue("diJetmass_cut", 500.0);
  double tauMass_cut = params->GetValue("tauMass_cut", 100.0);
  double MET_cut = params->GetValue("MET_cut", 50.0);
  double transmass_cut = params->GetValue("transmass_cut", 50.0);
  crateHistoMasps(nDir);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  MissingET *METpointer;

  for(Int_t entry = 0; entry < numberOfEntries; ++entry){

    treeReader->ReadEntry(entry);
    int pass_cuts[nDir] = {0};
    TLorentzVector jetLeadingVec(0., 0., 0., 0.); //tc stands for tau channel
    TLorentzVector jetSleadingVec(0., 0., 0., 0.);

    TLorentzVector Jet_1(0., 0., 0., 0.);
    TLorentzVector Jet_2(0., 0., 0., 0.);

    TLorentzVector jet_i(0., 0., 0., 0.);
    TLorentzVector elec_i(0., 0., 0., 0.);
    TLorentzVector muon_i(0., 0., 0., 0.);
    TLorentzVector tmp_tlv (0., 0., 0., 0.); //Vector used to swap vectors

    TLorentzVector Tau1HadTLV (0., 0., 0., 0.);
    TLorentzVector Tau2HadTLV (0., 0., 0., 0.);
    TLorentzVector Tau3HadTLV (0., 0., 0., 0.);

    Track *track_tau1;
    vector<TLorentzVector> jetsList;

    bool fill_tau1 = false;
    bool fill_tau2 = false;
    bool fill_tau3 = false;

    METpointer = (MissingET*) branchMissingET->At(0);
    double MET = METpointer->MET;
    double Met_phi = METpointer->Phi;
    double transmass = 0.;
    int njets_counter = 0;
    int ntau_counter = 0;
    int nmuon_counter = 0;
    double DiJetMass_final = 100.;
    int nBJets = 0;

    //////////////////Tau Channel///////////
    //Search for taus and bjets
    for (int j = 0; j < branchJet->GetEntriesFast(); j++){

      Jet *jet = (Jet*) branchJet->At(j);

      if((jet->BTag == 1) && (jet->PT > b_jet_pt_min)){nBJets++;}

      //Tau search
      if((jet->TauTag == 1) && (jet->PT > tau_pt_min)){

        ntau_counter++;
        double tau_energy = calculateE(jet->Eta, jet->PT, jet->Mass);

        if((fill_tau1 == false) && (fill_tau2 == false) && (fill_tau3 == false)){

          Tau1HadTLV.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
          fill_tau1 = true;
          continue;
        }

        if((fill_tau1 == true) && (fill_tau2 == false) && (fill_tau3 == false)){
          Tau2HadTLV.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
          fill_tau2 = true;
          continue;
        }

        if((fill_tau1 == true) && (fill_tau2 == true) && (fill_tau3 == false)){
          Tau3HadTLV.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
          fill_tau3 = true;
          continue;
        }
      }
    }

    bool tau1_muon_overlap = false;
    bool tau2_muon_overlap = false;
    bool tau3_muon_overlap = false;
    bool tau1_elec_overlap = false;
    bool tau2_elec_overlap = false;
    bool tau3_elec_overlap = false;

    //Check if taus overlap with muons
    for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){

      if(tau1_muon_overlap && tau2_muon_overlap && tau3_muon_overlap){break;}

      Muon *muon = (Muon*) branchMuon->At(muo);
      double muon_energy = calculateE(muon->Eta, muon->PT, 0.1056583715);
      muon_i.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);

      double DR_tau1_muon = Tau1HadTLV.DeltaR(muon_i);
      double DR_tau2_muon = Tau2HadTLV.DeltaR(muon_i);
      double DR_tau3_muon = Tau3HadTLV.DeltaR(muon_i);

      if (DR_tau3_muon < DR_jet_lep_max){
        tau3_muon_overlap = true;
        Tau3HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }

      if(DR_tau2_muon < DR_jet_lep_max){
        tau2_muon_overlap = true;
        Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }

      if(DR_tau1_muon < DR_jet_lep_max){
        tau1_muon_overlap = true;
        Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }
    } // for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++)

    //Check if taus overlap with electrons
    for (int el = 0; el < branchElectron->GetEntriesFast(); el++){

      if(tau1_elec_overlap && tau2_elec_overlap && tau3_elec_overlap){break;}

      Electron *elec = (Electron*) branchElectron->At(el);
      double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
      elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);

      double DR_tau1_elec = Tau1HadTLV.DeltaR(elec_i);
      double DR_tau2_elec = Tau2HadTLV.DeltaR(elec_i);
      double DR_tau3_elec = Tau3HadTLV.DeltaR(elec_i);

      if (DR_tau3_elec < DR_jet_lep_max){
        tau3_elec_overlap = true;
        Tau3HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }

      if(DR_tau2_elec < DR_jet_lep_max){
        tau2_elec_overlap = true;
        Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }

      if(DR_tau1_elec < DR_jet_lep_max){
        tau1_elec_overlap = true;
        Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }
    }


    //Check if jets overlap with taus
    int n_jets = 0;

    for (int l = 0; l < branchJet->GetEntriesFast(); l++){

      Jet *jet = (Jet*) branchJet->At(l);

      if((jet->PT > jet_min_pt) && (jet->TauTag == 0) && (jet->BTag == 0)){

        double jet_i_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
        jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_i_energy);

        //Check if jet overlaps with any of the taus
        double DR_tau1_jet = Tau1HadTLV.DeltaR(jet_i);
        double DR_tau2_jet = Tau2HadTLV.DeltaR(jet_i);
        double DR_tau3_jet = Tau3HadTLV.DeltaR(jet_i);

        if((DR_tau1_jet < DR_jet_lep_max)){
          Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          ntau_counter--;
        }

        if(DR_tau2_jet < DR_jet_lep_max){
          Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          ntau_counter--;
        }

        if(DR_tau3_jet < DR_jet_lep_max){
          Tau3HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          ntau_counter--;
        }

        if((n_jets <= 6) && (jet->TauTag == 0) && (jet->BTag == 0) && (abs(jet->Eta) < 5.0)){
          jetsList.push_back(jet_i);
          n_jets++;
        }
      }
    }


    //Order taus by pt
    for(int pt_order = 0; pt_order < 3; pt_order++){
      if(Tau2HadTLV.Pt() < Tau3HadTLV.Pt()){
        tmp_tlv = Tau2HadTLV;
        Tau2HadTLV = Tau3HadTLV;
        Tau3HadTLV = tmp_tlv;
      }
      if(Tau1HadTLV.Pt() < Tau2HadTLV.Pt()){
        tmp_tlv = Tau1HadTLV;
        Tau1HadTLV = Tau2HadTLV;
        Tau2HadTLV = tmp_tlv;
      }
    }

    double tau1_track_DR_min = 999.;
    double tau1_track_DR;
    //Search tau track
    for(int i = 0; i < branchTrack->GetEntriesFast(); i++){
      Track *track = (Track*) branchTrack->At(i);
      tau1_track_DR = calculate_deltaR(Tau1HadTLV, track);

      if(tau1_track_DR < tau1_track_DR_min){
        tau1_track_DR_min = tau1_track_DR;
        track_tau1 = track;
      }
    }


    int dijet_index1 = 0;
    int dijet_index2 = 0;

    //Search DiJetMass
    for(int k = 0; k < jetsList.size(); k++){

      if (jetsList.size() < 4){break;}
      Jet_1 = jetsList[k];

      if ((Jet_1.Pt() < VBF_jetPt_min) || (abs(Jet_1.Eta()) > 5.0)){continue;}

      for (int sj = k + 1; sj < jetsList.size(); sj++){

        if (sj != k){

          Jet_2 = jetsList[sj];
          if ((Jet_2.Pt() < VBF_jetPt_min) || (abs(Jet_2.Eta()) > 5.0)){continue;}
          double DiJetMass = (Jet_1+Jet_2).M();
          if (DiJetMass > DiJetMass_final){
            DiJetMass_final = DiJetMass;
            jetLeadingVec = Jet_1;
            jetSleadingVec = Jet_2;
            dijet_index1 = k;
            dijet_index2 = sj;
          }
        }
      }
    }

    if(jetLeadingVec.Pt() < jetSleadingVec.Pt()){
      tmp_tlv = jetLeadingVec;
      jetLeadingVec = jetSleadingVec;
      jetSleadingVec = tmp_tlv;
    }

    if(jetsList.size() > 3){

      //Remove the diJet pair from the jets list
      tmp_tlv = jetsList[jetsList.size() - 1];
      jetsList[jetsList.size() - 1] = jetsList[dijet_index1];
      jetsList[dijet_index1] = tmp_tlv;

      tmp_tlv = jetsList[jetsList.size() - 2];
      jetsList[jetsList.size() -2] = jetsList[dijet_index2];
      jetsList[dijet_index2] = tmp_tlv;

      jetsList.pop_back();
      jetsList.pop_back();

      //Order jets by pt
      for(int jet_order = 0; jet_order < 4; jet_order++){
        if(jetsList[2].Pt() < jetsList[3].Pt()){
          tmp_tlv = jetsList[2];
          jetsList[2] = jetsList[3];
          jetsList[3] = tmp_tlv;
        }
        if(jetsList[1].Pt() < jetsList[2].Pt()){
          tmp_tlv = jetsList[1];
          jetsList[1] = jetsList[2];
          jetsList[2] = tmp_tlv;
        }
        if(jetsList[0].Pt() < jetsList[1].Pt()){
          tmp_tlv = jetsList[0];
          jetsList[0] = jetsList[1];
          jetsList[1] = tmp_tlv;
        }
      }
    }

    //Check for jets pt condition
    int jet_pt_condition = 0;
    for (int i = 0; i < jetsList.size(); i++) {
      if(jetsList[i].Pt() > jet_min_pt){
        jet_pt_condition++;
      }
    }

    double delta_eta_diJet = abs(jetLeadingVec.Eta()-jetSleadingVec.Eta());
    double tauMass = (Tau1HadTLV + Tau2HadTLV).M();
    transmass = TMath::Sqrt(TMath::Abs(2*Tau1HadTLV.Pt()*MET*(1-TMath::Cos(normalizedDphi(Tau1HadTLV.Phi() - Met_phi)))));
    double ht = 0.;
    double st = 0.;

    if (jet_pt_condition > 1) {

      ht += jetLeadingVec.Pt() + jetSleadingVec.Pt();

      for (int i = 0; i < jetsList.size(); i++) {

        ht += jetsList[i].Pt();

      }
    }

    if (jet_pt_condition > 1) {

      st += Tau1HadTLV.Pt() + Tau2HadTLV.Pt() + Tau3HadTLV.Pt() + jetLeadingVec.Pt() + jetSleadingVec.Pt();

      for (int i = 0; i < jetsList.size(); i++) {
	      st += jetsList[i].Pt();
      }
    }

    //////// Apply cuts /////////
    // Events with no cuts
    pass_cuts[0] = 1;

    // Events with 2 taus with min pt
    if ((Tau1HadTLV.Pt() > tau_pt_cut) && (Tau2HadTLV.Pt() > tau_pt_cut)){
      pass_cuts[1] = 1;
    }
    // Events with 2 taus with max eta
    if((pass_cuts[1] == 1) && (abs(Tau1HadTLV.Eta()) < tau_eta_cut) && (abs(Tau2HadTLV.Eta()) < tau_eta_cut)){
      pass_cuts[2] = 1;
    }
    // Min tau system mass
    if ((pass_cuts[2] == 1) && (tauMass > tauMass_cut)){
      pass_cuts[3] = 1;
    }
    //Min MET cut
    if ((pass_cuts[3] == 1) && (MET > MET_cut)){
      pass_cuts[4] = 1;
    }
    // Number of bjets cut
    if ((pass_cuts[4] == 1) && (nBJets == 0)){
      pass_cuts[5] = 1;
    }
    // Jets with min pt
    if ((pass_cuts[5] == 1) && (jet_pt_condition > 1)){
      pass_cuts[6] = 1;
    }
    //Transverse mass cut
    if((pass_cuts[6] == 1) && (transmass > transmass_cut)){
      pass_cuts[7] = 1;
    }
    // Opposite hemisfere in dijet cut
    if((pass_cuts[7] == 1) && ((jetLeadingVec.Eta()*jetSleadingVec.Eta()) < 0) ){
      pass_cuts[8] = 1;
    }
    // Delta eta in dijet pair cut
    if ((pass_cuts[8] == 1) && (delta_eta_diJet > deltaEta_diJet_cut)){
      pass_cuts[9] = 1;
    }
    //Min DiJetMass cut
    if ((pass_cuts[9] == 1) && (DiJetMass_final > diJetmass_cut)){
      pass_cuts[10] = 1;
    }


    //Fill histograms
    for (int i = 0; i < nDir; i++){

      _hmap_Nevents[i]->Fill(0.0);
      _hmap_n_jets[i]->Fill(n_jets);
      _hmap_n_tau[i]->Fill(ntau_counter);

      if (pass_cuts[i] == 1){

        _hmap_Nevents[i]->Fill(1.0);

        if(jetLeadingVec.Pt() > 1.0){
	  _hmap_lead_jet_pT[i]->Fill(jetLeadingVec.Pt());
	  _hmap_lead_jet_eta[i]->Fill(jetLeadingVec.Eta());
	  _hmap_lead_jet_phi[i]->Fill(jetLeadingVec.Phi());
	  _hmap_slead_jet_pT[i]->Fill(jetSleadingVec.Pt());
	  _hmap_slead_jet_eta[i]->Fill(jetSleadingVec.Eta());
	  _hmap_slead_jet_phi[i]->Fill(jetSleadingVec.Phi());
        }
	if(Tau1HadTLV.Pt() > 1.0){
	  _hmap_tau1_pT[i]->Fill(Tau1HadTLV.Pt());
	  _hmap_tau1_eta[i]->Fill(Tau1HadTLV.Eta());
	  _hmap_tau1_phi[i]->Fill(Tau1HadTLV.Phi());
        }
        if(Tau2HadTLV.Pt() > 1.0){
	  _hmap_tau2_pT[i]->Fill(Tau2HadTLV.Pt());
	  _hmap_tau2_eta[i]->Fill(Tau2HadTLV.Eta());
	  _hmap_tau2_phi[i]->Fill(Tau2HadTLV.Phi());
        }
  if(tauMass > 0){_hmap_tauMass[i]->Fill(tauMass);}
  if(MET > 0){_hmap_MET[i]->Fill(MET);}
  if(transmass > 0){_hmap_transmass[i]->Fill(transmass);}
	if(ht > 0.0){_hmap_ht[i]->Fill(ht);}
        if(st > 0.0){_hmap_st[i]->Fill(st);}
	if(DiJetMass_final > 0.0){_hmap_dijet_mass[i]->Fill(DiJetMass_final);}
	if(delta_eta_diJet > 0){_hmap_dijet_deltaEta[i]->Fill(delta_eta_diJet);}
      }
    }
  }// end entry loop for tau channel

  theFile->cd();
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      _hmap_Nevents[d]->Write();
      _hmap_n_jets[d]->Write();
      _hmap_n_tau[d]->Write();
      _hmap_lead_jet_pT[d]->Write();
      _hmap_lead_jet_eta[d]->Write();
      _hmap_lead_jet_phi[d]->Write();
      _hmap_slead_jet_pT[d]->Write();
      _hmap_slead_jet_eta[d]->Write();
      _hmap_slead_jet_phi[d]->Write();
      _hmap_tau1_pT[d]->Write();
      _hmap_tau1_eta[d]->Write();
      _hmap_tau1_phi[d]->Write();
      _hmap_tau2_pT[d]->Write();
      _hmap_tau2_eta[d]->Write();
      _hmap_tau2_phi[d]->Write();
      _hmap_tauMass[d]->Write();
      _hmap_MET[d]->Write();
      _hmap_transmass[d]->Write();
      _hmap_ht[d]->Write();
      _hmap_st[d]->Write();
      _hmap_dijet_mass[d]->Write();
      _hmap_dijet_deltaEta[d]->Write();

    }
  theFile->Close();
}

PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

double PhenoAnalysis::calculateE(double eta, double pt, double mass){

  double theta = TMath::ATan(TMath::Exp(-eta));
  double sin_theta = TMath::Sin(2*theta);
  double p= pt/sin_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));

  return e;

}

double PhenoAnalysis::calculate_deltaR(TLorentzVector vector,  Track* track){

  double eta1 = vector.Eta();
  double phi1 = vector.Phi();
  double eta2 = track->Eta;
  double phi2 = track->Phi;
  double deltaR = sqrt(pow(eta1-eta2,2) + pow(phi1-phi2,2));
  return deltaR;
}

double PhenoAnalysis::normalizedDphi(double phi){
  const double PI  = 3.141592653589793238463;
  double twoPI = 2.0*PI;
  if ( phi < -PI ){phi += twoPI;}
  if ( phi > PI ){phi = twoPI-phi;}
  else phi = TMath::Abs(phi);
  return phi;
}
void PhenoAnalysis::crateHistoMasps (int directories)
{
  for (int i = 0; i < directories; i++)
    {
      _hmap_Nevents[i]       = new TH1F("Nevents", "Nevents", 3,0,3);
      _hmap_lead_jet_pT[i]   = new TH1F("jetLeadPt",    "p_{T} Leading Jet", 200, 0., 2000.);
      _hmap_lead_jet_eta[i]  = new TH1F("jetLeadEta",   "#eta Leading Jet", 50, -5.0, 5.0);
      _hmap_lead_jet_phi[i]  = new TH1F("jetLeadPhi",   "#phi Leading Jet", 70, -3.6, 3.6);
      _hmap_slead_jet_pT[i]   = new TH1F("jetSleadPt",    "p_{T} Sub-leading Jet", 200, 0., 2000.);
      _hmap_slead_jet_eta[i]  = new TH1F("jetSleadEta",   "#eta Sub-leading Jet", 50, -5.0, 5.0);
      _hmap_slead_jet_phi[i]  = new TH1F("jetSleadPhi",   "#phi Sub-leading Jet", 70, -3.6, 3.6);
      _hmap_n_jets[i]        = new TH1F("nJets",         "N(jet)", 6, 0, 6);
      _hmap_n_tau[i]         = new TH1F("nTaus",          "N(#tau)", 4, 0, 4);
      _hmap_tau1_pT[i]       = new TH1F("tau1Pt",        "p_{T}(#tau_{1})", 200, 0., 2000.);
      _hmap_tau1_eta[i]      = new TH1F("tau1Eta",       "#eta(#tau_{1})", 50, -3.5, 3.5);
      _hmap_tau1_phi[i]      = new TH1F("tau1Phi",       "#phi(#tau_{1})", 70, -3.6, 3.6);
      _hmap_tau2_pT[i]       = new TH1F("tau2Pt",        "p_{T}(#tau_{2})", 200, 0., 2000.);
      _hmap_tau2_eta[i]      = new TH1F("tau2Eta",       "#eta(#tau_{2})", 50, -3.5, 3.5);
      _hmap_tau2_phi[i]      = new TH1F("tau2Phi",       "#phi(#tau_{2})", 70, -3.6, 3.6);
      _hmap_tauMass[i]       = new TH1F("tauMass", "m(#t)", 100, 0, 1000);
      _hmap_MET[i]           = new TH1F("MET", "MET", 100, 0, 1000);
      _hmap_transmass[i]     = new TH1F("transMass", "Transverse Mass", 100, 0, 1000);
      _hmap_ht[i]            = new TH1F("HT", "H_{T}", 100, 0, 5000);
      _hmap_st[i]            = new TH1F("ST", "S_{T}", 100, 0, 5000);
      _hmap_dijet_mass[i]     = new TH1F("diJetMass", "diJet_Mass", 100, 0, 5000);
      _hmap_dijet_deltaEta[i] = new TH1F("diJetDeltaEta", "diJet_deltaEta", 160, 0, 8);
    }
}
