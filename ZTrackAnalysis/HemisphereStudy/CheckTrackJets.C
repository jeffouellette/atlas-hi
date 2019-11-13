#include <Utilities.h>
#include <GlobalParams.h>

using namespace atlashi;

void CheckTrackJets (TString _inFileName) {

  TFile* inFile = new TFile (_inFileName, "read");
  TTree* tree = (TTree*) inFile->Get ("tree");

  int njet = 0;
  float jetpt[20] = {};
  float jetphi[20] = {};
  float jeteta[20] = {};

  float zpt=0, zphi=0, zy=0, zm=0;

  tree->SetBranchAddress ("z_pt", &zpt);
  tree->SetBranchAddress ("z_phi",&zphi);
  tree->SetBranchAddress ("z_y",  &zy);
  tree->SetBranchAddress ("z_m",  &zm);
  tree->SetBranchAddress ("njt",  &njet);
  tree->SetBranchAddress ("jtpt", jetpt);
  tree->SetBranchAddress ("jtphi",jetphi);
  tree->SetBranchAddress ("jteta",jeteta);

  TH1D* h_zjetdphi = new TH1D ("h_zjetdphi", ";#Delta#phi;Y (#Delta#phi)", 12, -pi/2, 3*pi/2);
  h_zjetdphi->Sumw2 ();

  TH1D* h_zleadjetdphi = new TH1D ("h_zleadjetdphi", ";#Delta#phi;Y (#Delta#phi)", 12, -pi/2, 3*pi/2);
  h_zleadjetdphi->Sumw2 ();

  TH1D* h_jetjetdphi = new TH1D ("h_jetjetdphi", ";#Delta#phi;Y (#Delta#phi)", 12, -pi/2, 3*pi/2);
  h_jetjetdphi->Sumw2 ();

  TH1D* h_xzj = new TH1D ("h_xzj", ";#it{x}_{J}^{Z};Counts", 20, 0, 2);
  h_xzj->Sumw2 ();

  TH1D* h_xj = new TH1D ("h_xj", ";#it{x}_{J};Counts", 20, 0, 2);
  h_xj->Sumw2 ();
  

  for (int iEvt = 0; iEvt < tree->GetEntries (); iEvt++) {
    tree->GetEntry (iEvt);

    for (int iJ = 0; iJ < njet; iJ++) {
      //if (jetpt[iJ] < 20)
      //  continue;

      float dphi = DeltaPhi (zphi, jetphi[iJ], true);
      if (dphi < -pi/2) dphi = dphi + 2*pi;

      h_zjetdphi->Fill (dphi);

      dphi = DeltaPhi (zphi, jetphi[iJ]);
      float xzj = jetpt[iJ] * cos (pi-dphi) / zpt;
      h_xzj->Fill (xzj);
      //if (1 < xzj && xzj < 1.1)
      //  cout << iEvt << ", ";

      for (int iJ2 = 0; iJ2 < iJ; iJ2++) {
        //if (jetpt[iJ2] < 20)
        //  continue;
        dphi = DeltaPhi (jetphi[iJ], jetphi[iJ2], true);
        if (dphi < -pi/2) dphi = dphi + 2*pi;

        h_jetjetdphi->Fill (dphi);
      }
    }

    int lJ = -1, slJ = -1;
    for (int iJ = 0; iJ < njet; iJ++) {
      //if (jetpt[iJ] < 20)
      //  continue;

      if (lJ == -1 || jetpt[lJ] < jetpt[iJ]) {
        slJ = lJ;
        lJ = iJ;
      }
      else if (slJ == -1 || jetpt[slJ] < jetpt[iJ])
        slJ = iJ;
    }

    if (slJ != -1) {
      float dphi = DeltaPhi (jetphi[lJ], jetphi[slJ]);
      h_xj->Fill (jetpt[slJ] * cos (pi-dphi) / jetpt[lJ]);
    }

    if (lJ != -1) {
      float dphi = DeltaPhi (zphi, jetphi[lJ], true);
      if (dphi < -pi/2) dphi = dphi + 2*pi;

      h_zleadjetdphi->Fill (dphi);
    }
  }
  //cout << endl;

  TCanvas* c1 = new TCanvas ("c1", "", 800, 800);

  h_zjetdphi->Scale (1./tree->GetEntries ());
  h_zjetdphi->GetYaxis ()->SetRangeUser (0, 0.6);
  h_zjetdphi->SetLineColor (kBlue+1);
  h_zjetdphi->SetMarkerColor (kBlue+1);
  h_zjetdphi->Draw ("e1");

  h_zleadjetdphi->Scale (1./tree->GetEntries ());
  h_zleadjetdphi->GetYaxis ()->SetRangeUser (0, 0.6);
  h_zleadjetdphi->SetLineColor (kRed+1);
  h_zleadjetdphi->SetMarkerColor (kRed+1);
  h_zleadjetdphi->Draw ("same e1");
  
  h_jetjetdphi->Scale (1./tree->GetEntries ());
  h_jetjetdphi->GetYaxis ()->SetRangeUser (0, 0.6);
  h_jetjetdphi->SetLineColor (kGreen+2);
  h_jetjetdphi->SetMarkerColor (kGreen+2);
  h_jetjetdphi->Draw ("same e1");

  myText (0.22, 0.9, kBlue+1, "Z-jet correlations", 0.04);
  myText (0.22, 0.84, kRed+1, "Z-leading jet correlations", 0.04);
  myText (0.22, 0.78, kGreen+2, "Jet-jet correlations", 0.04);


  TCanvas* c2 = new TCanvas ("c2", "", 800, 800);
  h_xzj->SetLineColor (kBlack);
  h_xzj->SetMarkerColor (kBlack);
  h_xzj->SetMarkerStyle (kOpenCircle);
  h_xzj->Draw ("e1");

  h_xj->SetLineColor (kRed+1);
  h_xj->SetMarkerColor (kRed+1);
  h_xj->Draw ("same e1");

  myText (0.55, 0.9, kRed+1, "Lead. jet / sublead. jet (#it{x}_{J})", 0.04);
  myText (0.55, 0.84, kBlack, "Inclusive Z / jet (#it{x}_{J}^{Z})", 0.04);
}
