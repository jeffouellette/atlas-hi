

void PlotTrackMatchingProb () {

  TFile* inFile = new TFile ("../rootFiles/TrackingEfficiencies/Nominal/trackingEfficiencies_18.root", "read");

  TH1D* h_pp = (TH1D*) inFile->Get ("h_truth_matching_prob_pp");
  TH1D* h_PbPb = (TH1D*) inFile->Get ("h_truth_matching_prob_PbPb");

  h_pp->SetMarkerStyle (kFullCircle);
  h_PbPb->SetMarkerStyle (kOpenSquare);

  h_pp->SetMarkerColor (kRed+1);
  h_pp->SetLineColor (kRed+1);
  h_PbPb->SetMarkerColor (kAzure+2);
  h_PbPb->SetLineColor (kAzure+2);

  h_pp->GetYaxis ()->SetRangeUser (0.1, 1e10);
  h_PbPb->GetYaxis ()->SetRangeUser (0.1, 1e10);

  TCanvas* c1 = new TCanvas ("c1", "", 800, 800);
  gPad->SetLogy ();

  h_pp->Draw ("e1");
  h_PbPb->Draw ("e1 same");

  TLine* l0p3 = new TLine (0.3, 0.1, 0.3, 1e5);
  TLine* l0p5 = new TLine (0.5, 0.1, 0.5, 1e5);
  TLine* l0p6 = new TLine (0.6, 0.1, 0.6, 1e5);

  l0p3->SetLineWidth (1);
  l0p3->SetLineStyle (2);
  l0p5->SetLineWidth (1);
  l0p5->SetLineStyle (2);
  l0p6->SetLineWidth (1);
  l0p6->SetLineStyle (2);

  l0p3->Draw ("same");
  l0p5->Draw ("same");
  l0p6->Draw ("same");

  myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Simulation", 0.044);
  myMarkerTextNoLine (0.25, 0.83, kRed+1, kFullCircle, "Powheg+Pythia Z#rightarrow ll", 1.4, 0.032);
  myMarkerTextNoLine (0.25, 0.78, kAzure+2, kOpenSquare, "Powheg+Pythia Z#rightarrow ll with overlay", 1.4, 0.032);

  const float i_pp_0p3 = h_pp->Integral (h_pp->FindFixBin (0.3), h_pp->FindFixBin (1)) / h_pp->Integral (0, 1);
  const float i_pp_0p5 = h_pp->Integral (h_pp->FindFixBin (0.5), h_pp->FindFixBin (1)) / h_pp->Integral (0, 1);
  const float i_pp_0p6 = h_pp->Integral (h_pp->FindFixBin (0.6), h_pp->FindFixBin (1)) / h_pp->Integral (0, 1);
  const float i_PbPb_0p3 = h_PbPb->Integral (h_PbPb->FindFixBin (0.3), h_PbPb->FindFixBin (1)) / h_PbPb->Integral (0, 1);
  const float i_PbPb_0p5 = h_PbPb->Integral (h_PbPb->FindFixBin (0.5), h_PbPb->FindFixBin (1)) / h_PbPb->Integral (0, 1);
  const float i_PbPb_0p6 = h_PbPb->Integral (h_PbPb->FindFixBin (0.6), h_PbPb->FindFixBin (1)) / h_PbPb->Integral (0, 1);

  myText (0.35, 0.270+0.40, kBlack, "f_{ch}(> 0.3)", 0.026);
  myText (0.50, 0.270+0.40, kBlack, "f_{ch}(> 0.5)", 0.026);
  myText (0.62, 0.270+0.40, kBlack, "f_{ch}(> 0.6)", 0.026);

  myText (0.35, 0.235+0.40, kRed+1, Form ("%.3f%%", i_pp_0p3), 0.026);
  myText (0.50, 0.235+0.40, kRed+1, Form ("%.3f%%", i_pp_0p5), 0.026);
  myText (0.62, 0.235+0.40, kRed+1, Form ("%.3f%%", i_pp_0p6), 0.026);
  myText (0.35, 0.200+0.40, kAzure+2, Form ("%.4f%%", i_PbPb_0p3), 0.026);
  myText (0.50, 0.200+0.40, kAzure+2, Form ("%.4f%%", i_PbPb_0p5), 0.026);
  myText (0.62, 0.200+0.40, kAzure+2, Form ("%.4f%%", i_PbPb_0p6), 0.026);
}
