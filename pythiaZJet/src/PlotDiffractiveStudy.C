#ifndef __PlotDiffractiveStudy_C__
#define __PlotDiffractiveStudy_C__

void PlotDiffractiveStudy () {
  TFile* ndFile = new TFile ("../nondiffractive_only.root", "read");
  TTree* ndTree = (TTree*) ndFile->Get ("tree");

  TFile* allFile = new TFile ("../nondiffractive_and_diffractive.root", "read");
  TTree* allTree = (TTree*) allFile->Get ("tree");

  int code = 0;
  int part_n = 0;
  vector<float>* part_pt = nullptr, *part_eta = nullptr, *part_phi = nullptr;

  ndTree->SetBranchAddress ("code",     &code);
  ndTree->SetBranchAddress ("part_n",   &part_n);
  ndTree->SetBranchAddress ("part_pt",  &part_pt);
  ndTree->SetBranchAddress ("part_eta", &part_eta);
  ndTree->SetBranchAddress ("part_phi", &part_phi);

  allTree->SetBranchAddress ("code",     &code);
  allTree->SetBranchAddress ("part_n",   &part_n);
  allTree->SetBranchAddress ("part_pt",  &part_pt);
  allTree->SetBranchAddress ("part_eta", &part_eta);
  allTree->SetBranchAddress ("part_phi", &part_phi);

  const double ptBins[7] = {1., 2., 4., 8., 15., 30., 60.};

  TH1D* h_nd = new TH1D ("h_trk_pt_nd", ";#it{p}_{T}^{ ch} [GeV]", 6, ptBins);
  TH1D* h_all = new TH1D ("h_trk_pt_all", ";#it{p}_{T}^{ ch} [GeV]", 6, ptBins);
  int n_nd = 0;
  int n_all = 0;

  TH1D* h_nch_nondiff = new TH1D ("h_nch_nondiff", "N_{ch}", 10, -0.5, 9.5);
  TH1D* h_nch_diff = new TH1D ("h_nch_diff", "N_{ch}", 10, -0.5, 9.5);

  for (int iEvt = 0; iEvt < ndTree->GetEntries (); iEvt++) {
    ndTree->GetEntry (iEvt);

    if (part_n < 2)
      continue;

    n_nd++;

    for (int iPart = 0; iPart < part_n; iPart++) {
      if (fabs (part_eta->at (iPart)) > 2.5)
        continue;
      h_nd->Fill (part_pt->at (iPart));
    }
  }

  for (int iEvt = 0; iEvt < allTree->GetEntries (); iEvt++) {
    allTree->GetEntry (iEvt);

    if (code != 101)
      h_nch_diff->Fill (part_n);
    else
      h_nch_nondiff->Fill (part_n);

    if (part_n < 2)
      continue;

    n_all++;

    for (int iPart = 0; iPart < part_n; iPart++) {
      if (fabs (part_eta->at (iPart)) > 2.5)
        continue;
      h_all->Fill (part_pt->at (iPart));
    }
  }

  cout << ndTree->GetEntries () << " events in non-diffractive-only tree" << endl;
  cout << allTree->GetEntries () << " events in inclusive events tree" << endl;

  h_nd->Scale (1. / (n_nd * TMath::Pi ()), "width");
  h_all->Scale (1. / (n_all * TMath::Pi ()), "width");

  TFile* oFile = new TFile ("diff_hists.root", "recreate");
  h_nd->Write ("", TObject::kOverwrite);
  h_all->Write ("", TObject::kOverwrite);

  oFile->Close ();


  TCanvas* c1 = new TCanvas ("c1", "", 800, 800);
  {
    c1->cd ();

    TPad* uPad = new TPad ("uPad", "", 0.0, 0.4, 1.0, 1.0);
    TPad* dPad = new TPad ("dPad", "", 0.0, 0.0, 1.0, 0.4);
    uPad->SetBottomMargin (0);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    uPad->Draw ();
    dPad->Draw ();


    uPad->cd ();
    uPad->SetLogx ();
    uPad->SetLogy ();

    h_nd->SetMarkerColor (kRed);
    h_nd->SetLineColor (kRed);
    h_nd->SetMarkerStyle (kOpenCircle);
    h_all->SetMarkerColor (kBlack);
    h_all->SetLineColor (kBlack);
    h_all->SetMarkerStyle (kOpenCircle);

    h_nd->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T}d#Delta#phi [GeV^{-1}]");
    h_nd->GetXaxis ()->SetMoreLogLabels ();

    h_nd->GetXaxis ()->SetTitleFont (43);
    h_nd->GetXaxis ()->SetTitleSize (24);
    h_nd->GetXaxis ()->SetLabelFont (43);
    h_nd->GetXaxis ()->SetLabelSize (24);

    h_nd->GetYaxis ()->SetTitleFont (43);
    h_nd->GetYaxis ()->SetTitleSize (24);
    h_nd->GetYaxis ()->SetLabelFont (43);
    h_nd->GetYaxis ()->SetLabelSize (24);

    h_nd->GetXaxis ()->SetTitleOffset (3.0);
    h_nd->GetYaxis ()->SetTitleOffset (2.0);

    h_nd->Draw ("e1");
    h_all->Draw ("e1 same");


    dPad->cd ();
    dPad->SetLogx ();
    TH1D* h_rat = (TH1D*) h_nd->Clone ("h_rat");
    h_rat->Divide (h_all);

    h_rat->SetMarkerColor (kBlack);
    h_rat->SetLineColor (kBlack);
    h_rat->SetMarkerStyle (kOpenCircle);
    h_rat->GetYaxis ()->SetTitle ("Non-diff. / All processes");
    h_rat->GetYaxis ()->CenterTitle ();
    h_rat->GetYaxis ()->SetRangeUser (0.7, 2);
    h_rat->GetXaxis ()->SetMoreLogLabels ();

    h_rat->GetXaxis ()->SetTitleFont (43);
    h_rat->GetXaxis ()->SetTitleSize (24);
    h_rat->GetXaxis ()->SetLabelFont (43);
    h_rat->GetXaxis ()->SetLabelSize (24);

    h_rat->GetYaxis ()->SetTitleFont (43);
    h_rat->GetYaxis ()->SetTitleSize (24);
    h_rat->GetYaxis ()->SetLabelFont (43);
    h_rat->GetYaxis ()->SetLabelSize (24);

    h_rat->GetXaxis ()->SetTitleOffset (3.0);
    h_rat->GetYaxis ()->SetTitleOffset (2.0);

    h_rat->Draw ("e1");

    TLine* line = new TLine (0,0,0,0);
    line->SetLineColor (kBlack);
    line->SetLineStyle (2);
    line->DrawLine (1, 1, 60, 1);

    line->SetLineColor (kPink-8);
    line->SetLineWidth (1);
    line->DrawLine (1, 1.05, 60, 1.05);
    line->DrawLine (1, 0.95, 60, 0.95);
  }



  TCanvas* c2 = new TCanvas ("c2", "", 800, 800);
  {
    h_nch_diff->Draw ("hist");

    TPad* uPad = new TPad ("uPad", "", 0.0, 0.4, 1.0, 1.0);
    TPad* dPad = new TPad ("dPad", "", 0.0, 0.0, 1.0, 0.4);
    uPad->SetBottomMargin (0);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    uPad->Draw ();
    dPad->Draw ();


    uPad->cd ();
    uPad->SetLogy ();

    h_nch_diff->SetMarkerColor (kRed);
    h_nch_diff->SetLineColor (kRed);
    h_nch_diff->SetMarkerStyle (kOpenCircle);
    h_nch_nondiff->SetMarkerColor (kBlack);
    h_nch_nondiff->SetLineColor (kBlack);
    h_nch_nondiff->SetMarkerStyle (kOpenCircle);

    h_nch_diff->GetXaxis ()->SetTitle ("N_{ch}");
    h_nch_diff->GetYaxis ()->SetTitle ("N_{event}");
    h_nch_diff->GetXaxis ()->SetMoreLogLabels ();

    h_nch_diff->GetXaxis ()->SetTitleFont (43);
    h_nch_diff->GetXaxis ()->SetTitleSize (24);
    h_nch_diff->GetXaxis ()->SetLabelFont (43);
    h_nch_diff->GetXaxis ()->SetLabelSize (24);

    h_nch_diff->GetYaxis ()->SetTitleFont (43);
    h_nch_diff->GetYaxis ()->SetTitleSize (24);
    h_nch_diff->GetYaxis ()->SetLabelFont (43);
    h_nch_diff->GetYaxis ()->SetLabelSize (24);

    h_nch_diff->GetXaxis ()->SetTitleOffset (3.0);
    h_nch_diff->GetYaxis ()->SetTitleOffset (2.0);

    h_nch_diff->Draw ("e1");
    h_nch_nondiff->Draw ("e1 same");

    dPad->cd ();
    TH1D* h_rat = (TH1D*) h_nch_diff->Clone ("h_rat");
    h_rat->Divide (h_nch_nondiff);

    h_rat->SetMarkerColor (kBlack);
    h_rat->SetLineColor (kBlack);
    h_rat->SetMarkerStyle (kOpenCircle);
    h_rat->GetYaxis ()->SetTitle ("Non-diff. / All processes");
    h_rat->GetYaxis ()->CenterTitle ();
    h_rat->GetYaxis ()->SetRangeUser (0.7, 2);
    h_rat->GetXaxis ()->SetMoreLogLabels ();

    h_rat->GetXaxis ()->SetTitleFont (43);
    h_rat->GetXaxis ()->SetTitleSize (24);
    h_rat->GetXaxis ()->SetLabelFont (43);
    h_rat->GetXaxis ()->SetLabelSize (24);

    h_rat->GetYaxis ()->SetTitleFont (43);
    h_rat->GetYaxis ()->SetTitleSize (24);
    h_rat->GetYaxis ()->SetLabelFont (43);
    h_rat->GetYaxis ()->SetLabelSize (24);

    h_rat->GetXaxis ()->SetTitleOffset (3.0);
    h_rat->GetYaxis ()->SetTitleOffset (2.0);

    h_rat->Draw ("e1");

    TLine* line = new TLine (0,0,0,0);
    line->SetLineColor (kBlack);
    line->SetLineStyle (2);
    line->DrawLine (1, 1, 60, 1);

    line->SetLineColor (kPink-8);
    line->SetLineWidth (1);
    line->DrawLine (1, 1.05, 60, 1.05);
    line->DrawLine (1, 0.95, 60, 0.95);
  }

}

#endif
