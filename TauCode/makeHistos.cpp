#include <string>

using namespace std;

void makeHistos (){

  const char *filePath = "/home/carlosmiguelpatino/Documents/VBF_Analysis/TauSamples/signal/output_signal.root";

  TFile *f = new TFile(filePath);

  const char *variables[] = { "jet_lead_pT", "jet_lead_eta", "jet_slead_pt", "jet_slead_eta", "tau1_pt", "tau1_eta", "tau2_pt", "tau2_eta", "HT", "ST" };

  const char *directories[] = {"No_cuts", "N_bjets", "VBF_diJetMass"};

  for (size_t i = 0; i < sizeof(directories); i++) {

    for (size_t j = 0; j < sizeof(variables); j++) {

      string directory = str(directories[i]);
      string variable = str(variables[j]);
      string histogramPath = directory + "/" + variable;
      string outputPath = "/home/carlosmiguelpatino/Documents/Plots_VBFAnalysis/QCD0" + directories[i] + variables[j] + ".png";

      TCanvas *c = new TCanvas;

      TH1F *histogram = (TH1F*)f->Get(histogramPath);

      histogram->GetXaxis()->SetRangeUser(0.,900.);
      histogram->Draw();

      TImage *img = TImage::Create();

      img->FromPad(c);
      img->WriteImage(outputPath);

      delete histogram;
      delete c;
      delete img;
    }
  }
}
