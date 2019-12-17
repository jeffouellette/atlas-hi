

void MakeNchWeights (const bool do2018 = true) {

  TFile* f_pythia = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/eventWeightsFile.root", "read");
  TH1D* h_pythia = (TH1D*) f_pythia->Get ("h_PbPbNchDist_mc");

  TFile* f_hijing = new TFile (Form ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MinbiasAnalysis/PbPb_Hijing_%i_eventWeights.root", do2018 ? 18:15), "read");
  TH1D* h_hijing = (TH1D*) f_hijing->Get ("h_PbPbNchDist_hijing");

  TH1D* pythia_weights = (TH1D*) h_pythia->Clone ("h_PythiaNchWeights");
  pythia_weights->Divide (h_hijing, h_pythia);

  TFile* f_out = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/pythiaNchWeights.root", "recreate");
  pythia_weights->Write ();

  f_pythia->Close ();
  f_hijing->Close ();
  f_out->Close ();
}
