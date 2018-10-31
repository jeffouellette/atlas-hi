const double etabins[5] = {-2.5, -1.8, 0, 1.8, 2.5};


void MakePlots () {

  TFile* dataFile = new TFile ("myOut_pPb_data_perf_0.root", "read");
  TFile* overlayFile = new TFile ("myOut_pPb_mc_pythia8_perf_0.overlay.root", "read");
  TFile* signalFile = new TFile ("myOut_pPb_mc_pythia8_perf_0.signal.root", "read");

  TH2D* h_data_rtrk1 = (TH2D*)dataFile->Get("h_rtrk1");
  TH2D* h_overlay_rtrk1 = (TH2D*)overlayFile->Get("h_rtrk1");
  TH2D* h_signal_rtrk1 = (TH2D*)signalFile->Get("h_rtrk1");

  TH1D* h_data_rtrk1_pt[4];
  TH1D* h_overlay_rtrk1_pt[4];
  TH1D* h_signal_rtrk1_pt[4];
 
  for (int iEta = 0; iEta < 4; iEta++) {
    h_data_rtrk1_pt[iEta] = (TH1D*)h_data_rtrk1->ProjectionY (Form ("data_eta%i", iEta), iEta, iEta);
    h_overlay_rtrk1_pt[iEta] = (TH1D*)h_overlay_rtrk1->ProjectionY (Form ("overlay_eta%i", iEta), iEta, iEta);
    h_signal_rtrk1_pt[iEta] = (TH1D*)h_signal_rtrk1->ProjectionY (Form ("signal_eta%i", iEta), iEta, iEta);
  }

  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);
  canvas->Draw();

  for (int iEta = 0; iEta < 4; iEta++) {
    h_data_rtrk1_pt[iEta]->SetMarkerColor (kBlack);
    h_data_rtrk1_pt[iEta]->SetLineColor (kBlack);
    h_data_rtrk1_pt[iEta]->Draw ("e1");
    myMarkerText (0.68, 0.79, kBlack, kFullCircle, "Data", 1.25, 0.04);

    h_overlay_rtrk1_pt[iEta]->SetMarkerColor(kBlue);
    h_overlay_rtrk1_pt[iEta]->SetLineColor(kBlue);
    h_overlay_rtrk1_pt[iEta]->Draw("same e1");
    myMarkerText (0.68, 0.73, kBlue, kFullCircle, "Pythia8 with Overlay", 1.25, 0.04);

    h_signal_rtrk1_pt[iEta]->SetMarkerColor (kRed);
    h_signal_rtrk1_pt[iEta]->SetLineColor (kRed);
    h_signal_rtrk1_pt[iEta]->Draw("same e1");
    myMarkerText (0.68, 0.67, kRed, kFullCircle, "Pythia8 Signal Only", 1.25, 0.04);

    myText (0.68, 0.85, kBlack, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 0.04);

    canvas->SaveAs (Form ("data_rtrk1_eta%i.pdf", iEta));
  }
}
