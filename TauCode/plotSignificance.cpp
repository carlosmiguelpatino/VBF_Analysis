#include <string>
#include <math.h>

using namespace std;

void plotSignificance(string normalizedDirectories, const unsigned n_points, const unsigned initial_value, const unsigned delta){

  double variable_values [n_points];
  double significance_values [n_points];

  std::cout << "Value,Significance" << '\n';

  for (int i = 0; i < n_points; i++){

    variable_values[i] = (i*delta) + initial_value;
    string variable_path = Form("%d", (i*delta) + initial_value);

    string fileDY= normalizedDirectories +  variable_path + "/VBF_diJetMass/normalizedHistos_DYToLL.root";
    string fileW= normalizedDirectories +  variable_path + "/VBF_diJetMass/normalizedHistos_wjets.root";
    string fileTtbar= normalizedDirectories +  variable_path + "/VBF_diJetMass/normalizedHistos_ttbar_semi.root";
    string fileSignal= normalizedDirectories +  variable_path + "/VBF_diJetMass/normalizedHistos_signal.root";

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

    cout << (i*delta) + initial_value << "," <<significance << '\n';

  }

  TCanvas *c = new TCanvas;

  TGraph *g = new TGraph(n_points,variable_values,significance_values);
  g->SetTitle("Significance vs diJet Pair Mass");
  g->GetXaxis()->SetTitle("m(jj)");
  //g->GetYaxis()->SetTitle("Significance");
  g->Draw("AC*");

  TImage *img = TImage::Create();

  img->FromPad(c);
  img->WriteImage("../TauPlots/significance.png");

}
