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
  int nDir = 12;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("Taus_pT_min");
  theDirectory[2]  = HistoOutputFile->mkdir("Taus_eta_min");
  theDirectory[3]  = HistoOutputFile->mkdir("Impact_parameter_track1_min");
  theDirectory[4]  = HistoOutputFile->mkdir("Taus_mass_min");
  theDirectory[5]  = HistoOutputFile->mkdir("MET_min");
  theDirectory[6]  = HistoOutputFile->mkdir("N_bjets");
  theDirectory[7]  = HistoOutputFile->mkdir("Jets_pT_min");
  theDirectory[8]  = HistoOutputFile->mkdir("Transmass_min");
  theDirectory[9]  = HistoOutputFile->mkdir("VBF_jets_opposite_hemispheres");
  theDirectory[10]  = HistoOutputFile->mkdir("VBF_jets_delta");
  theDirectory[11]  = HistoOutputFile->mkdir("VBF_diJetMass");
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

  double muon_pt_min      = params->GetValue ("muon_pt_min", 10.);
  double b_jet_pt_min     = params->GetValue ("b_jet_pt_min", 30.0);
  double DR_jet_lep_max  = params->GetValue ("DR_jet_lep_max", 0.3);
  double jet_min_pt       = params->GetValue ("jet_min_pt", 15);
  double muon_pt_cut       = params->GetValue("muon_pt_cut", 20.);
  double muon_eta_cut      = params->GetValue ("muon_eta_cut", 2.1);
  double deltaEta_diJet_cut = params->GetValue("deltaEta_diJet_cut", 3.8);
  double diJetmass_cut    = params->GetValue("diJetmass_cut", 500.0);
  double MET_cut = params->GetValue("MET_cut", 50.0);
  double transmass_cut = params->GetValue("transmass_cut", 50.0);
  crateHistoMasps(nDir);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  MissingET *METpointer;

  for(Int_t entry = 0; entry < numberOfEntries; ++entry){

    ////////////////Muon Channel///////////
    treeReader->ReadEntry(entry);
    int pass_cuts[nDir] = {0};
    TLorentzVector Jet_leading_vec(0., 0., 0., 0.); //mc stands for muon channel
    TLorentzVector Jet_sleading_vec(0., 0., 0., 0.);

    TLorentzVector Jet_1(0., 0., 0., 0.);
    TLorentzVector Jet_2(0., 0., 0., 0.);

    TLorentzVector jet_i(0., 0., 0., 0.);
    TLorentzVector elec_i(0., 0., 0., 0.);
    TLorentzVector tmp_tlv (0., 0., 0., 0.); //Vector used to swap vectors

    TLorentzVector Muon1TLV(0., 0., 0., 0.);
    TLorentzVector Muon2TLV(0., 0., 0., 0.);
    TLorentzVector Muon3TLV(0., 0., 0., 0.);

    vector<TLorentzVector> jetsList;

    bool fill_muon1 = false;
    bool fill_muon2 = false;
    bool fill_muon3 = false;

    METpointer = (MissingET*) branchMissingET->At(0);
    double MET = METpointer->MET;
    double Met_phi = METpointer->Phi;
    int njets_counter = 0;
    int nmuon_counter = 0;
    double DiJetMass_final = 100.;
    int n_b_jets = 0;
    bool four_jets_condition = false;

    //Search for muons
    for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){

      nmuon_counter++;

      Muon *muon = (Muon*) branchMuon->At(muo);

      double muon_energy = calculateE(muon->Eta, muon->PT, 0.1056583715);
      if((fill_muon1 == false) && (fill_muon2 == false) && (fill_muon3 == false)){
        Muon1TLV.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
        fill_muon1 = true;
        continue;
      }

      if((fill_muon1 == true) && (fill_muon2 == false) && (fill_muon3 == false)){
        Muon2TLV.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
        fill_muon2 = true;
        continue;
      }

      if((fill_muon1 == true) && (fill_muon2 == true) && (fill_muon3 == false)){
        Muon3TLV.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
        fill_muon3 = true;
        continue;
      }
    }

    bool muon1_elec_overlap = false;
    bool muon2_elec_overlap = false;
    bool muon3_elec_overlap = false;

    //Check if muons overlap with electrons
    for (int el = 0; el < branchElectron->GetEntriesFast(); el++){

      if(muon1_elec_overlap && muon2_elec_overlap && muon3_elec_overlap){break;}

      Electron *elec = (Electron*) branchElectron->At(el);
      double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
      elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);

      double DR_muon1_elec = Muon1TLV.DeltaR(elec_i);
      double DR_muon2_elec = Muon2TLV.DeltaR(elec_i);
      double DR_muon3_elec = Muon3TLV.DeltaR(elec_i);


      if (DR_muon3_elec < DR_jet_lep_max){
        muon3_elec_overlap = true;
        Muon3TLV.SetPtEtaPhiE(0.,0.,0.,0.);
        nmuon_counter--;
      }

      if(DR_muon2_elec < DR_jet_lep_max){
        muon2_elec_overlap = true;
        Muon2TLV.SetPtEtaPhiE(0.,0.,0.,0.);
        nmuon_counter--;
      }

      if(DR_muon1_elec < DR_jet_lep_max){
        muon1_elec_overlap = true;
        Muon1TLV.SetPtEtaPhiE(0.,0.,0.,0.);
        nmuon_counter--;
      }
    }

    //Check if jets overlap with muons

    int n_jets = 0;

    for (int l = 0; l < branchJet->GetEntriesFast(); l++){

      Jet *jet = (Jet*) branchJet->At(l);

      if(jet->BTag == 1){n_b_jets++;}

      double jet_i_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
      jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_i_energy);

      //Check if jet overlaps with any of the taus
      double DR_muon1_jet = Muon1TLV.DeltaR(jet_i);
      double DR_muon2_jet = Muon2TLV.DeltaR(jet_i);
      double DR_muon3_jet = Muon3TLV.DeltaR(jet_i);

      if((DR_muon1_jet < DR_jet_lep_max)){
        Muon1TLV.SetPtEtaPhiE(0.,0.,0.,0.);
        nmuon_counter--;
      }

      if(DR_muon2_jet < DR_jet_lep_max){
        Muon2TLV.SetPtEtaPhiE(0.,0.,0.,0.);
        nmuon_counter--;
      }

      if(DR_muon3_jet < DR_jet_lep_max){
        Muon3TLV.SetPtEtaPhiE(0.,0.,0.,0.);
        nmuon_counter--;
      }

      if((n_jets <= 6) && (jet->TauTag == 0) && (jet->BTag == 0) && (jet_i.Pt() > jet_min_pt)){
        jetsList.push_back(jet_i);
        n_jets++;
      }
    }

    //Order muons by pt
    for(int pt_order = 0; pt_order < 3; pt_order++){
      if(Muon2TLV.Pt() < Muon3TLV.Pt()){
        tmp_tlv = Muon2TLV;
        Muon2TLV = Muon3TLV;
        Muon3TLV = tmp_tlv;
      }
      if(Muon1TLV.Pt() < Muon2TLV.Pt()){
        tmp_tlv = Muon1TLV;
        Muon1TLV = Muon2TLV;
        Muon2TLV = tmp_tlv;
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
          if ((Jet_2.Pt() < 30.0) || (abs(Jet_2.Eta()) > 5.0)){continue;}
          double DiJetMass = (Jet_1+Jet_2).M();
          if (DiJetMass > DiJetMass_final){
            DiJetMass_final = DiJetMass;
            Jet_leading_vec = Jet_1;
            Jet_sleading_vec = Jet_2;
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

      tmp_tlv = jetsList[jetsList.size() - 1];
      jetsList[jetsList.size() - 1] = jetsList[dijet_index1];
      jetsList[dijet_index1] = tmp_tlv;

      tmp_tlv = jetsList[jetsList.size() - 2];
      jetsList[jetsList.size() - 2] = jetsList[dijet_index2];
      jetsList[dijet_index2] = tmp_tlv;

      jetsList.pop_back();
      jetsList.pop_back();



      //Order jets by pt
      for(int jet_order = 0; jet_order < 4; jet_order++){
        if(jetsList[2].Pt() < jetsList[3].Pt()){
          tmp_tlv = jetsList[2];;
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
    bool pass_jet_min_pt = true;
    for (int i = 0; i < jetsList.size(); i++) {
      if(jetsList[i].Pt() > 30.){
        jet_pt_condition++;
      }
    }

    double delta_eta_diJet = abs(Jet_leading_vec.Eta()-Jet_sleading_vec.Eta());
    //double tauMass = (Tau1HadTLV + Tau2HadTLV).M();
    double ht = 0.;
    transmass = TMath::Sqrt(TMath::Abs(2*Muon1TLV.Pt()*MET*(1-TMath::Cos(normalizedDphi(Muon1TLV.Phi() - Met_phi)))));
    double st = 0.;

    if (jet_pt_condition > 1) {
      st += Muon1TLV.Pt() + Muon2TLV.Pt() + Jet_leading_vec.Pt() + Jet_sleading_vec.Pt();

      for (int i = 0; i < jetsList.size(); i++) {
        ht += jetsList[i].Pt();
      }
    }


    //////// Apply cuts /////////
    // Events with no cuts
    pass_cuts[0] = 1;

    // Events with 2 taus with min pt
    if ((Muon1TLV.Pt() > tau_pt_cut) && (Muon2TLV.Pt() > muon_pt_cut)){
      pass_cuts[1] = 1;
    }
    // Events with 2 taus with max eta
    if((pass_cuts[1] == 1) && (abs(Muon1TLV.Eta()) < muon_eta_cut) && (abs(Muon2TLV.Eta()) < muon_eta_cut)){
      pass_cuts[2] = 1;
    }
    //Min impact para meter tau1 cut
    if((pass_cuts[2] == 1)){
      pass_cuts[3] = 1;
    }
    // Min tau system mass
    if ((pass_cuts[3] == 1)/* && (tauMass > tauMass_cut)*/){
      pass_cuts[4] = 1;
    }
    //Min MET cut
    if ((pass_cuts[4] == 1) && (MET > MET_cut)){
      pass_cuts[5] = 1;
    }
    // Number of bjets cut
    if ((pass_cuts[5] == 1) && (nBJets == 0)){
      pass_cuts[6] = 1;
    }
    // Jets with min pt
    if ((pass_cuts[6] == 1) && (jet_pt_condition > 1)){
      pass_cuts[7] = 1;
    }
    //Transverse mass cut
    if((pass_cuts[7] == 1) && (transmass > transmass_cut)){
      pass_cuts[8] = 1;
    }
    // Opposite hemisfere in dijet cut
    if((pass_cuts[8] == 1) && ((jetLeadingVec.Eta()*jetSleadingVec.Eta()) < 0) ){
      pass_cuts[9] = 1;
    }
    // Delta eta in dijet pair cut
    if ((pass_cuts[9] == 1) && (delta_eta_diJet > deltaEta_diJet_cut)){
      pass_cuts[10] = 1;
    }
    //Min DiJetMass cut
    if ((pass_cuts[10] == 1) && (DiJetMass_final > diJetmass_cut)){
      pass_cuts[11] = 1;
    }

    //Fill histograms
    for (int i = 0; i < nDir; i++){
      _hmap_Nevents[i]->Fill(0.0);
      _hmap_n_jets[i]->Fill(n_jets);
      _hmap_n_muon[i]->Fill(nmuon_counter);

      if (pass_cuts[i] == 1){
        _hmap_Nevents[i]->Fill(1.0);
        if(Jet_leading_vec.Pt() > 1.0){
	  _hmap_lead_jet_pT[i]->Fill(Jet_leading_vec.Pt());
	  _hmap_lead_jet_eta[i]->Fill(Jet_leading_vec.Eta());
	  _hmap_lead_jet_phi[i]->Fill(Jet_leading_vec.Phi());
          _hmap_slead_jet_pT[i]->Fill(Jet_sleading_vec.Pt());
          _hmap_slead_jet_eta[i]->Fill(Jet_sleading_vec.Eta());
          _hmap_slead_jet_phi[i]->Fill(Jet_sleading_vec.Phi());
        }
	_hmap_muon1_pT[i]->Fill(Muon1TLV.Pt());
	_hmap_muon1_eta[i]->Fill(Muon1TLV.Eta());
	_hmap_muon1_phi[i]->Fill(Muon1TLV.Phi());
	_hmap_muon2_pT[i]->Fill(Muon2TLV.Pt());
	_hmap_muon2_eta[i]->Fill(Muon2TLV.Eta());
	_hmap_muon2_phi[i]->Fill(Muon2TLV.Phi());
        if(ht > 0.0){_hmap_ht[i]->Fill(ht);}
        if(DiJetMass_final > 0.0){_hmap_dijet_mass[i]->Fill(DiJetMass_final);}
        if(delta_eta_diJet > 0.){_hmap_dijet_deltaEta[i]->Fill(delta_eta_diJet);}
      }
    }
  } //end entry loop for muon channel

  theFile->cd();
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      _hmap_Nevents[d]->Write();
      _hmap_n_jets[d]->Write();
      _hmap_n_muon[d]->Write();
      _hmap_lead_jet_pT[d]->Write();
      _hmap_lead_jet_eta[d]->Write();
      _hmap_lead_jet_phi[d]->Write();
      _hmap_slead_jet_pT[d]->Write();
      _hmap_slead_jet_eta[d]->Write();
      _hmap_slead_jet_phi[d]->Write();
      _hmap_muon1_pT[d]->Write();
      _hmap_muon1_eta[d]->Write();
      _hmap_muon1_phi[d]->Write();
      _hmap_muon2_pT[d]->Write();
      _hmap_muon2_eta[d]->Write();
      _hmap_muon2_phi[d]->Write();
      _hmap_ht[d]->Write();
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
      _hmap_Nevents[i]       = new TH1F("Nevents", "Nevents", 3,0.,3);
      _hmap_lead_jet_pT[i]   = new TH1F("jet_lead_pT",    "Pt leading jet", 100, 0., 1000.);
      _hmap_lead_jet_eta[i]  = new TH1F("jet_lead_eta",   "#eta jet", 50, -5.0, 5.0);
      _hmap_lead_jet_phi[i]  = new TH1F("jet_lead_phi",   "#phi jet", 70, -3.6, 3.6);
      _hmap_slead_jet_pT[i]   = new TH1F("jet_slead_pT",    "Pt sleading jet", 100, 0., 1000.);
      _hmap_slead_jet_eta[i]  = new TH1F("jet_slead_eta",   "#eta jet", 50, -5.0, 5.0);
      _hmap_slead_jet_phi[i]  = new TH1F("jet_slead_phi",   "#phi jet", 70, -3.6, 3.6);
      _hmap_n_jets[i]        = new TH1F("N_jets",         "N(jet)", 4, 0., 4);
      _hmap_n_muon[i]         = new TH1F("N_taus",          "", 4, 0., 4);
      _hmap_muon1_pT[i]       = new TH1F("muon1_pT",        "p_{T}(#muon_{1})", 200, 0., 600.);
      _hmap_muon1_eta[i]      = new TH1F("muon1_eta",       "#eta(#muon_{1})", 50, -3.5, 3.5);
      _hmap_muon1_phi[i]      = new TH1F("muon1_phi",       "#phi(#muon_{1})", 70, -3.6, 3.6);
      _hmap_muon2_pT[i]       = new TH1F("muon2_pT",        "p_{T}(#muon_{2})", 150, 0., 300.);
      _hmap_muon2_eta[i]      = new TH1F("muon2_eta",       "#eta(#muon_{2})", 50, -3.5, 3.5);
      _hmap_muon2_phi[i]      = new TH1F("muon2_phi",       "#phi(#muon_{2})", 70, -3.6, 3.6);
      _hmap_ht[i]            = new TH1F("HT", "HT", 200, 0, 1000);
      _hmap_dijet_mass[i]     = new TH1F("diJet_Mass", "diJet_Mass", 1000, 0, 5000);
      _hmap_dijet_deltaEta[i] = new TH1F("diJet_deltaEta", "diJet_deltaEta", 160, 0, 8);
    }
}
