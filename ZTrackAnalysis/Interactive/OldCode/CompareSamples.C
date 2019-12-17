const double pi = TMath::Pi ();

TH1D* oldEpt = nullptr;
TH1D* newEpt = nullptr;
TH1D* oldMpt = nullptr;
TH1D* newMpt = nullptr;
TH1D* oldTpt = nullptr;
TH1D* newTpt = nullptr;

TH2D* oldEetaphi = nullptr;
TH2D* newEetaphi = nullptr;
TH2D* oldMetaphi = nullptr;
TH2D* newMetaphi = nullptr;
TH2D* oldTetaphi = nullptr;
TH2D* newTetaphi = nullptr;


void CompareSamples () {

  TFile* oldFile = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/DataAnalysis/340718.root", "read");
  TFile* newFile = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/AlternateProcessing/340718.root", "read");

  TTree* oldTree = (TTree*)oldFile->Get ("ppZTrackTree");
  TTree* newTree = (TTree*)newFile->Get ("ppZTrackTree");

  bool isEE = false;
  int ntrk = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0, z_pt = 0, z_eta = 0, z_phi = 0, z_m = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;

  oldTree->SetBranchAddress ("isEE",    &isEE);
  oldTree->SetBranchAddress ("l1_pt",   &l1_pt);
  oldTree->SetBranchAddress ("l1_eta",  &l1_eta);
  oldTree->SetBranchAddress ("l1_phi",  &l1_phi);
  oldTree->SetBranchAddress ("l2_pt",   &l2_pt);
  oldTree->SetBranchAddress ("l2_eta",  &l2_eta);
  oldTree->SetBranchAddress ("l2_phi",  &l2_phi);
  oldTree->SetBranchAddress ("z_pt",    &z_pt);
  oldTree->SetBranchAddress ("z_eta",   &z_eta);
  oldTree->SetBranchAddress ("z_phi",   &z_phi);
  oldTree->SetBranchAddress ("z_m",     &z_m);
  oldTree->SetBranchAddress ("ntrk",    &ntrk);
  oldTree->SetBranchAddress ("trk_pt",  &trk_pt);
  oldTree->SetBranchAddress ("trk_eta", &trk_eta);
  oldTree->SetBranchAddress ("trk_phi", &trk_phi);

  newTree->SetBranchAddress ("isEE",    &isEE);
  newTree->SetBranchAddress ("l1_pt",   &l1_pt);
  newTree->SetBranchAddress ("l1_eta",  &l1_eta);
  newTree->SetBranchAddress ("l1_phi",  &l1_phi);
  newTree->SetBranchAddress ("l2_pt",   &l2_pt);
  newTree->SetBranchAddress ("l2_eta",  &l2_eta);
  newTree->SetBranchAddress ("l2_phi",  &l2_phi);
  newTree->SetBranchAddress ("z_pt",    &z_pt);
  newTree->SetBranchAddress ("z_eta",   &z_eta);
  newTree->SetBranchAddress ("z_phi",   &z_phi);
  newTree->SetBranchAddress ("z_m",     &z_m);
  newTree->SetBranchAddress ("ntrk",    &ntrk);
  newTree->SetBranchAddress ("trk_pt",  &trk_pt);
  newTree->SetBranchAddress ("trk_eta", &trk_eta);
  newTree->SetBranchAddress ("trk_phi", &trk_phi);

  oldEpt = new TH1D ("old_electron_pt", "", 30, 0, 250);
  newEpt = new TH1D ("new_electron_pt", "", 30, 0, 250);
  oldMpt = new TH1D ("old_muon_pt", "", 30, 0, 250);
  newMpt = new TH1D ("new_muon_pt", "", 30, 0, 250);
  oldTpt = new TH1D ("old_trk_pt", "", 30, 0, 250);
  newTpt = new TH1D ("new_trk_pt", "", 30, 0, 250);

  oldEetaphi = new TH2D ("old_electron_eta_phi", "", 50, -2.5, 2.5, 50, -pi, pi);
  newEetaphi = new TH2D ("new_electron_eta_phi", "", 50, -2.5, 2.5, 50, -pi, pi);
  oldMetaphi = new TH2D ("old_muon_eta_phi", "", 50, -2.5, 2.5, 50, -pi, pi);
  newMetaphi = new TH2D ("new_muon_eta_phi", "", 50, -2.5, 2.5, 50, -pi, pi);
  oldTetaphi = new TH2D ("old_trk_eta_phi", "", 50, -2.5, 2.5, 50, -pi, pi);
  newTetaphi = new TH2D ("new_trk_eta_phi", "", 50, -2.5, 2.5, 50, -pi, pi);

  oldEpt->Sumw2 ();
  newEpt->Sumw2 ();
  oldMpt->Sumw2 ();
  newMpt->Sumw2 ();
  oldTpt->Sumw2 ();
  newTpt->Sumw2 ();

  for (int i = 0; i < oldTree->GetEntries (); i++) {
    oldTree->GetEntry (i);
    if (isEE) {
      oldEpt->Fill (l1_pt);
      oldEpt->Fill (l2_pt);
      oldEetaphi->Fill (l1_eta, l1_phi);
      oldEetaphi->Fill (l2_eta, l2_phi);
    } else {
      oldMpt->Fill (l1_pt);
      oldMpt->Fill (l2_pt);
      oldMetaphi->Fill (l1_eta, l1_phi);
      oldMetaphi->Fill (l2_eta, l2_phi);
    }
    for (int t = 0; t < ntrk; t++) {
      oldTpt->Fill (trk_pt->at (t));
      oldTetaphi->Fill (trk_eta->at (t), trk_phi->at (t));
    }
  }
  for (int i = 0; i < newTree->GetEntries (); i++) {
    newTree->GetEntry (i);
    if (isEE) {
      newEpt->Fill (l1_pt);
      newEpt->Fill (l2_pt);
      newEetaphi->Fill (l1_eta, l1_phi);
      newEetaphi->Fill (l2_eta, l2_phi);
    } else {
      newMpt->Fill (l1_pt);
      newMpt->Fill (l2_pt);
      newMetaphi->Fill (l1_eta, l1_phi);
      newMetaphi->Fill (l2_eta, l2_phi);
    }
    for (int t = 0; t < ntrk; t++) {
      newTpt->Fill (trk_pt->at (t));
      newTetaphi->Fill (trk_eta->at (t), trk_phi->at (t));
    }
  }

  TCanvas* EPtCanvas = new TCanvas ("EPtCanvas", "", 800, 1200);
  EPtCanvas->Divide (1, 2);
  EPtCanvas->cd (1); 
  gPad->SetLogy ();
  oldEpt->SetLineColor (kBlack);
  newEpt->SetLineColor (kBlue);
  oldEpt->SetMarkerColor (kBlack);
  newEpt->SetMarkerColor (kBlue);
  oldEpt->GetXaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  oldEpt->Draw ("e1");
  newEpt->Draw ("e1 same");
  myText (0.65, 0.8, kBlue, "New", 0.06);
  myText (0.65, 0.65, kBlack, "Old", 0.06);
  EPtCanvas->cd (2);
  TH1D* ratioEpt = (TH1D*)newEpt->Clone ("ratioEpt");
  ratioEpt->Divide (oldEpt);
  ratioEpt->GetXaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  ratioEpt->GetYaxis ()->SetTitle ("New / Old");
  ratioEpt->Draw ("e1");
  EPtCanvas->SaveAs ("Plots/newPPrecon/electron_pt.pdf");
 
 
  TCanvas* MPtCanvas = new TCanvas ("MPtCanvas", "", 800, 1200);
  MPtCanvas->Divide (1, 2);
  MPtCanvas->cd (1); 
  gPad->SetLogy ();
  oldMpt->SetLineColor (kBlack);
  newMpt->SetLineColor (kBlue);
  oldMpt->SetMarkerColor (kBlack);
  newMpt->SetMarkerColor (kBlue);
  oldMpt->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  oldMpt->Draw ("e1");
  newMpt->Draw ("e1 same");
  myText (0.65, 0.8, kBlue, "New", 0.06);
  myText (0.65, 0.65, kBlack, "Old", 0.06);
  MPtCanvas->cd (2);
  TH1D* ratioMpt = (TH1D*)newMpt->Clone ("ratioMpt");
  ratioMpt->Divide (oldMpt);
  ratioMpt->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  ratioMpt->GetYaxis ()->SetTitle ("New / Old");
  ratioMpt->Draw ("e1");
  MPtCanvas->SaveAs ("Plots/newPPrecon/muon_pt.pdf");


  TCanvas* EetaphiCanvas = new TCanvas ("EetaphiCanvas", "", 800, 600);
  EetaphiCanvas->SetRightMargin (0.18);

  oldEetaphi->GetXaxis ()->SetTitle ("#eta");
  oldEetaphi->GetYaxis ()->SetTitle ("#phi");
  oldEetaphi->Draw ("colz");
  myText (0.6, 0.8, kBlack, "Old e's", 0.06);
  EetaphiCanvas->SaveAs ("Plots/newPPrecon/old_electron_etaphi.pdf");

  newEetaphi->GetXaxis ()->SetTitle ("#eta");
  newEetaphi->GetYaxis ()->SetTitle ("#phi");
  newEetaphi->Draw ("colz");
  myText (0.6, 0.8, kBlack, "New e's", 0.06);
  EetaphiCanvas->SaveAs ("Plots/newPPrecon/new_electron_etaphi.pdf");

  TH2D* ratioEetaphi = (TH2D*)newEetaphi->Clone ("ratioEetaphi");
  ratioEetaphi->Divide (oldEetaphi);
  ratioEetaphi->GetXaxis ()->SetTitle ("#eta");
  ratioEetaphi->GetYaxis ()->SetTitle ("#phi");
  ratioEetaphi->Draw ("colz");
  myText (0.6, 0.8, kBlack, "New / Old e's", 0.06);
  EetaphiCanvas->SaveAs ("Plots/newPPrecon/ratio_electron_etaphi.pdf");
  

  TCanvas* MetaphiCanvas = new TCanvas ("MetaphiCanvas", "", 800, 600);
  MetaphiCanvas->SetRightMargin (0.18);

  oldMetaphi->GetXaxis ()->SetTitle ("#eta");
  oldMetaphi->GetYaxis ()->SetTitle ("#phi");
  oldMetaphi->Draw ("colz");
  myText (0.6, 0.8, kBlack, "Old #mu's", 0.06);
  MetaphiCanvas->SaveAs ("Plots/newPPrecon/old_muon_etaphi.pdf");

  newMetaphi->GetXaxis ()->SetTitle ("#eta");
  newMetaphi->GetYaxis ()->SetTitle ("#phi");
  newMetaphi->Draw ("colz");
  myText (0.6, 0.8, kBlack, "New #mu's", 0.06);
  MetaphiCanvas->SaveAs ("Plots/newPPrecon/new_muon_etaphi.pdf");

  TH2D* ratioMetaphi = (TH2D*)newMetaphi->Clone ("ratioMetaphi");
  ratioMetaphi->Divide (oldMetaphi);
  ratioMetaphi->GetXaxis ()->SetTitle ("#eta");
  ratioMetaphi->GetYaxis ()->SetTitle ("#phi");
  ratioMetaphi->Draw ("colz");
  myText (0.6, 0.8, kBlack, "New / Old #mu's", 0.06);
  MetaphiCanvas->SaveAs ("Plots/newPPrecon/ratio_muon_etaphi.pdf");


  TCanvas* TetaphiCanvas = new TCanvas ("TetaphiCanvas", "", 800, 600);
  TetaphiCanvas->SetRightMargin (0.18);

  oldTetaphi->GetXaxis ()->SetTitle ("#eta");
  oldTetaphi->GetYaxis ()->SetTitle ("#phi");
  oldTetaphi->Draw ("colz");
  myText (0.6, 0.8, kBlack, "Old tracks", 0.06);
  TetaphiCanvas->SaveAs ("Plots/newPPrecon/old_track_etaphi.pdf");

  newTetaphi->GetXaxis ()->SetTitle ("#eta");
  newTetaphi->GetYaxis ()->SetTitle ("#phi");
  newTetaphi->Draw ("colz");
  myText (0.6, 0.8, kBlack, "New tracks", 0.06);
  TetaphiCanvas->SaveAs ("Plots/newPPrecon/new_track_etaphi.pdf");

  TH2D* ratioTetaphi = (TH2D*)newTetaphi->Clone ("ratioTetaphi");
  ratioTetaphi->Divide (oldTetaphi);
  ratioTetaphi->GetXaxis ()->SetTitle ("#eta");
  ratioTetaphi->GetYaxis ()->SetTitle ("#phi");
  ratioTetaphi->Draw ("colz");
  myText (0.6, 0.8, kBlack, "New / Old tracks", 0.06);
  TetaphiCanvas->SaveAs ("Plots/newPPrecon/ratio_track_etaphi.pdf");
}
