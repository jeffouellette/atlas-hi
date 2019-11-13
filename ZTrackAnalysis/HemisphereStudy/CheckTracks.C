#include <Utilities.h> 
#include <GlobalParams.h>

using namespace atlashi;

void CheckTracks () {
  TFile* f = new TFile ("outFile.root", "read");
  TTree* t = (TTree*) f->Get ("ppZTrackTree");

  bool isEE;
  float zpt = 0, zphi = 0, zy = 0, zm=0;
  float l1pt = 0, l1phi= 0, l1eta=0, l2pt = 0, l2phi = 0, l2eta = 0;
  vector<float>* tpt = nullptr, *teta = nullptr, *tphi = nullptr;

  t->SetBranchAddress ("isEE", &isEE);
  t->SetBranchAddress ("trk_pt", &tpt);
  t->SetBranchAddress ("trk_eta", &teta);
  t->SetBranchAddress ("trk_phi", &tphi);
  t->SetBranchAddress ("z_pt",&zpt);
  t->SetBranchAddress ("z_phi",&zphi);
  t->SetBranchAddress ("z_y",&zy);
  t->SetBranchAddress ("z_m",&zm);
  t->SetBranchAddress ("l1_pt", &l1pt);
  t->SetBranchAddress ("l1_phi",&l1phi);
  t->SetBranchAddress ("l1_eta",&l1eta);
  t->SetBranchAddress ("l2_pt",&l2pt);
  t->SetBranchAddress ("l2_eta",&l2eta);
  t->SetBranchAddress ("l2_phi",&l2phi);

  TH1D* h_z_trk_dphi_trig = new TH1D ("h_z_trk_dphi_trig", ";#Delta#phi;Y (#Delta#phi)", 12, -pi/2, 3*pi/2);
  h_z_trk_dphi_trig->Sumw2 ();
  int n_trig = 0;

  TH1D* h_z_trk_dphi_incl = new TH1D ("h_z_trk_dphi_incl", ";#Delta#phi;Y (#Delta#phi)", 12, -pi/2, 3*pi/2);
  h_z_trk_dphi_incl->Sumw2 ();
  int n_incl = 0;

  //vector<int> events = {8971, 10881, 23213, 28369, 30878, 35675, 37190, 41221, 42104, 43394, 49626, 51854, 51962, 52482, 52828, 57058, 65594, 67511, 78125, 79179, 81683, 91797, 100453, 102026, 103341, 108935, 114154, 115263, 116915, 117069, 121697, 127167, 130492, 132766, 138399, 138571};
  //for (int iEvt : events) {
  const int nEvts = t->GetEntries ();
  cout << "Looping over " << nEvts << " pp events..." << endl;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    t->GetEntry (iEvt);
    if (zpt < 25)
      continue;
    //if (!isEE)
    //  continue;

    bool triggerEvent = false;
    //float triggerPhi = 0;
    int triggerPart = 0;
    //cout << endl << "iEvt " << iEvt << endl;

    //TLorentzVector ztlv;
    //ztlv.SetPxPyPzE (zpt*TMath::Cos(zphi), zpt*TMath::Sin(zphi), TMath::Sqrt(zpt*zpt+zm*zm)*TMath::SinH(zy), TMath::Sqrt(zpt*zpt+zm*zm)*TMath::CosH(zy));

    //cout << "zpt = " << zpt << ",\t zphi/2pi = " << 0.1591549 * (InTwoPi ((double)zphi)) << ",\t zy = " << zy << ",\t zeta = " << ztlv.Eta () << ",\t zm = " << zm << endl << endl;
    //cout << "l1pt = " << l1pt << ",\t l1phi/2pi = " << 0.1591549 * (InTwoPi(l1phi)) << ",\t l1eta = " << l1eta << endl;
    //cout << "l2pt = " << l2pt << ",\t l2phi/2pi = " << 0.1591549 * InTwoPi(l2phi) << ",\t l2eta = " << l2eta << endl;

    //float sumptx = 0, sumpty = 0;
    for (int it = 0; it <tpt->size (); it++) {
      if (tpt->at (it) < 2)
        continue;

      //cout << "trkpt = " << tpt->at (it) << ",\t trkphi/2pi = " << 0.1591549 * InTwoPi (tphi->at (it)) << ",\t trketa = " << teta->at (it) << endl;
      //sumptx += tpt->at (it) * cos (tphi->at (it));
      //sumpty += tpt->at (it) * sin (tphi->at (it));
      float dphi = DeltaPhi (zphi, tphi->at (it));
 
      if (tpt->at (it) > 25 && dphi < pi/2) {
        triggerEvent = true;
        //triggerPhi = tphi->at (it);
        triggerPart = it;
      }

      dphi = DeltaPhi (zphi, tphi->at (it), true);
      if (dphi < -pi/2) dphi = dphi + 2*pi;

      h_z_trk_dphi_incl->Fill (dphi);
    }
    n_incl++;

    if (triggerEvent) {
      for (int it = 0; it <tpt->size (); it++) {
        if (it == triggerPart)
          continue;

        if (tpt->at (it) < 2)
          continue;

        float dphi = DeltaPhi (zphi, tphi->at (it));

        if (tpt->at (it) > 25 && dphi < pi/2)
          continue; // skip the trigger particle

        dphi = DeltaPhi (tphi->at (triggerPart), tphi->at (it), true);
        if (dphi < -pi/2) dphi = dphi + 2*pi;

        h_z_trk_dphi_trig->Fill (dphi);
      }
      n_trig++;
    }
    //cout << "sumptx - zpx= " << sumptx - ztlv.Px () << endl;
    //cout << "sumpty - zpy= " << sumpty - ztlv.Py () << endl;
    //cout << "Diff sq. = " << TMath::Sqrt (pow (sumptx-ztlv.Px (), 2) + pow (sumpty-ztlv.Py (), 2)) << endl;

    //h_mpt->Fill (TMath::Sqrt (pow (sumpty-ztlv.Py (), 2) + pow (sumptx-ztlv.Px (), 2))); }
  }
  cout << "Looping complete!" << endl;
  cout << "# inclusive events: " << n_incl << endl;
  cout << "# triggered events: " << n_trig << endl;


  h_z_trk_dphi_incl->Scale (1./n_incl);
  h_z_trk_dphi_trig->Scale (1./n_trig);

  TCanvas* c1 = new TCanvas ("c1", "", 800, 600);

  TGraphAsymmErrors* g = make_graph (h_z_trk_dphi_incl);
  ResetXErrors (g);

  g->GetXaxis ()->SetTitle ("#Delta#phi");
  g->GetYaxis ()->SetTitle ("Y (#Delta#phi)");

  g->GetYaxis ()->SetRangeUser (0, 4);

  g->SetMarkerStyle (kFullCircle);
  g->SetLineColor (kBlue+1);
  g->SetMarkerColor (kBlue+1);

  g->GetXaxis ()->SetTitleOffset (0.6);
  g->GetYaxis ()->SetTitleOffset (0.8);
  g->GetXaxis ()->SetTitleSize (0.08);
  g->GetYaxis ()->SetTitleSize (0.06);
  g->GetXaxis ()->SetLabelSize (0.06);
  g->GetYaxis ()->SetLabelSize (0.06);

  g->Draw ("AP");


  g = make_graph (h_z_trk_dphi_trig);
  ResetXErrors (g);

  g->GetXaxis ()->SetTitle ("#Delta#phi");
  g->GetYaxis ()->SetTitle ("Y (#Delta#phi)");

  g->GetYaxis ()->SetRangeUser (0, 4);

  g->SetMarkerStyle (kFullCircle);
  g->SetLineColor (kOrange+8);
  g->SetMarkerColor (kOrange+8);

  g->GetXaxis ()->SetTitleOffset (0.6);
  g->GetYaxis ()->SetTitleOffset (0.8);
  g->GetXaxis ()->SetTitleSize (0.08);
  g->GetYaxis ()->SetTitleSize (0.06);
  g->GetXaxis ()->SetLabelSize (0.06);
  g->GetYaxis ()->SetLabelSize (0.06);

  g->Draw ("P");


  TH1D* h_z_trk_dphi_diff = new TH1D ("h_z_trk_dphi_diff", ";#Delta#phi;Y (#Delta#phi)", 12, -pi/2, 3*pi/2);
  h_z_trk_dphi_diff->Sumw2 ();

  h_z_trk_dphi_diff->Add (h_z_trk_dphi_trig);
  h_z_trk_dphi_diff->Add (h_z_trk_dphi_incl, -1);

  g = make_graph (h_z_trk_dphi_diff);
  ResetXErrors (g);

  g->GetXaxis ()->SetTitle ("#Delta#phi");
  g->GetYaxis ()->SetTitle ("Y (#Delta#phi)");

  g->GetYaxis ()->SetRangeUser (0, 4);

  g->SetMarkerStyle (kOpenCircle);
  g->SetLineColor (kBlack);
  g->SetMarkerColor (kBlack);

  g->GetXaxis ()->SetTitleOffset (0.6);
  g->GetYaxis ()->SetTitleOffset (0.8);
  g->GetXaxis ()->SetTitleSize (0.08);
  g->GetYaxis ()->SetTitleSize (0.06);
  g->GetXaxis ()->SetLabelSize (0.06);
  g->GetYaxis ()->SetLabelSize (0.06);

  g->Draw ("P");

  myText (0.22, 0.90, kBlue+1, "Inclusive events", 0.04);
  myText (0.22, 0.84, kOrange+8, "Events w/ 25 GeV trk < #pi/2 from Z", 0.04);
  myText (0.22, 0.78, kBlack, "Difference", 0.04);

  myText (0.72, 0.90, kBlack, "#it{pp}, 5.02 TeV", 0.04);
  myText (0.72, 0.84, kBlack, "#it{p}_{T}^{Z} > 25 GeV", 0.04);
  myText (0.72, 0.78, kBlack, "#it{p}_{T}^{h} > 2 GeV", 0.04);
  
}
