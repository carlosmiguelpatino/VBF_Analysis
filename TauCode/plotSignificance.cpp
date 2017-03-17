#include <string>
#include <math.h>

using namespace std;

void plotSignificance(){


  string delta_names[]= {"38","39", "40", "41", "42", "43", "44"};

  int n_points = 8;
  double eta_values [] = { 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5};
  double significance_values [7];

  for (int i = 0; i <=6; i++){

    string fileDY="../dEtadiJetOpti/normalizedEta" +  delta_names[i] + "/VBF_diJetMass/normalizedHistos_DYToLL.root";
    string fileW="../dEtadiJetOpti/normalizedEta" +  delta_names[i] + "/VBF_diJetMass/normalizedHistos_wjets.root";
    string fileTtbar="../dEtadiJetOpti/normalizedEta" +  delta_names[i] + "/VBF_diJetMass/normalizedHistos_ttbar_semi.root";
    string fileSignal="../dEtadiJetOpti/normalizedEta" +  delta_names[i] + "/VBF_diJetMass/normalizedHistos_signal.root";

    TFile *DY_file = new TFile(fileDY.c_str());
    TFile *Wjets_file = new TFile(fileW.c_str());
    TFile *ttbar_file = new TFile(fileTtbar.c_str());
    TFile *signal_file = new TFile(fileSignal.c_str());

    TH1F * DY_histo= (TH1F*)DY_file->Get("HT_lumi");
    TH1F * W_histo= (TH1F*)Wjets_file->Get("HT_lumi");
    TH1F * ttbar_histo= (TH1F*)ttbar_file->Get("HT_lumi");
    TH1F * signal_histo= (TH1F*)signal_file->Get("HT_lumi");

    double DY_events = DY_histo->Integral();
    double W_events = W_histo->Integral();
    double ttbar_events = ttbar_histo->Integral();
    double signal_events = signal_histo->Integral();


    int nBackgrounds = DY_events + W_events + ttbar_events;

    double significance = signal_events/sqrt(nBackgrounds + signal_events);

    significance_values[i] = significance;

  }

  TCanvas *c = new TCanvas;

  TGraph *g = new TGraph(n_points,eta_values,significance_values);
  g->Draw("AC*");

  TImage *img = TImage::Create();

  img->FromPad(c);
  img->WriteImage("significance.png");

}
