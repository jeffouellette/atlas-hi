

void CompareData () {

  TFile* f_hitight = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/DataAnalysis/Variations/TrackHITightWPVariation/savedHists.root", "read");
  TH1D* h_hitight = (TH1D*) f_hitight->Get ("h_z_trk_raw_pt_mumu_iPtZ3_iPhi1_iCent3_data_trackHITightVar");
  TH1D* h_n_hitight = (TH1D*) f_hitight->Get ("h_z_counts_mumu_iPtZ3_iCent3_data_trackHITightVar");

  TFile* f_hiloose = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/savedHists.root", "read");
  TH1D* h_hiloose = (TH1D*) f_hiloose->Get ("h_z_trk_raw_pt_mumu_iPtZ3_iPhi1_iCent3_data");
  TH1D* h_n_hiloose = (TH1D*) f_hiloose->Get ("h_z_counts_mumu_iPtZ3_iCent3_data");

  h_hitight->Scale (1./h_n_hitight->GetBinContent (1));
  h_hiloose->Scale (1./h_n_hiloose->GetBinContent (1));

  h_hitight->SetLineColor (kBlack);
  h_hiloose->SetLineColor (kAzure+2);
  h_hitight->SetMarkerColor (kBlack);
  h_hiloose->SetMarkerColor (kAzure+2);

  TCanvas* c  = new TCanvas ("c", "", 800, 800);
  c->Divide (1,2);
  c->Draw ();

  c->cd (1);
  h_hitight->Draw ("e1");
  h_hiloose->Draw ("same e1");

  TH1D* ratio = (TH1D*) h_hitight->Clone ("ratio");
  ratio->Divide (h_hitight, h_hiloose);
  c->cd (2);
  ratio->Draw ("e1");

  
  

}
