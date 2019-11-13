#include <Utilities.h>
#include <GlobalParams.h>

using namespace atlashi;

void TreeSlimmer (TString _inFileName, TString _outFileName) {

  TFile* inFile = new TFile (_inFileName, "read");
  TTree* inTree = (TTree*) inFile->Get ("ppZTrackTree");

  TFile* outFile = new TFile (_outFileName, "recreate");
  TTree* outTree = new TTree ("tree", "Zs + Jets in pp data");
  outTree->SetDirectory (outFile);

  int rn, en;
  float zpt, zphi, zy, zeta, zm, missingPx, missingPy;
  int njet = 0;
  vector<float>* jpt = nullptr, *jeta = nullptr, *jphi = nullptr, *je = nullptr, *trkpt = nullptr, *trketa = nullptr, *trkphi = nullptr;
  float jetpt[10] = {};
  float jeteta[10] = {};
  float jetphi[10] = {};
  float jete[10] = {};

  inTree->SetBranchAddress ("run_number", &rn);
  inTree->SetBranchAddress ("event_number", &en);
  inTree->SetBranchAddress ("z_pt", &zpt);
  inTree->SetBranchAddress ("z_y", &zy);
  inTree->SetBranchAddress ("z_m", &zm);
  inTree->SetBranchAddress ("z_phi", &zphi);
  inTree->SetBranchAddress ("njet", &njet);
  inTree->SetBranchAddress ("jet_pt", &jpt);
  inTree->SetBranchAddress ("jet_eta", &jeta);
  inTree->SetBranchAddress ("jet_phi", &jphi);
  inTree->SetBranchAddress ("jet_e", &je);
  inTree->SetBranchAddress ("trk_pt", &trkpt);
  inTree->SetBranchAddress ("trk_eta", &trketa);
  inTree->SetBranchAddress ("trk_phi", &trkphi);

  outTree->Branch ("run_number", &rn);
  outTree->Branch ("event_number", &en);
  outTree->Branch ("z_pt", &zpt);
  outTree->Branch ("z_eta", &zeta);
  outTree->Branch ("z_y", &zy);
  outTree->Branch ("z_m", &zm);
  outTree->Branch ("z_phi", &zphi);
  outTree->Branch ("njt", &njet);
  outTree->Branch ("jtpt", jetpt, "jtpt[njt]/F");
  outTree->Branch ("jteta", jeteta, "jteta[njt]/F");
  outTree->Branch ("jtphi", jetphi, "jtphi[njt]/F");
  outTree->Branch ("jte", jete, "jte[njt]/F");
  outTree->Branch ("missingPx", &missingPx);
  outTree->Branch ("missingPy", &missingPy);


  const int nEvt = inTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvt; iEvt++) {
    inTree->GetEntry (iEvt);

    if (zpt < 25)
      continue;

    TLorentzVector tlv;
    tlv.SetPtEtaPhiM (zpt*cos(zphi), zpt*sin(zphi), sqrt(zpt*zpt+zm*zm)*sinh(zy), sqrt(zpt*zpt+zm*zm)*cosh(zy));
    zeta = tlv.Eta ();

    missingPx = zpt*cos (zphi);
    missingPy = zpt*sin (zphi);

    bool triggerEvent = false;
    for (int iP = 0; iP < trkpt->size (); iP++) {
      if (trkpt->at (iP) > 20 && DeltaPhi (trkphi->at (iP), zphi) < pi/2)
        triggerEvent = true;
    }
    if (!triggerEvent)
      continue;

    //njet = jpt->size ();
    for (int iJ = 0; iJ < njet; iJ++) {
      jetpt[iJ] = jpt->at (iJ);
      jeteta[iJ] = jeta->at (iJ);
      jetphi[iJ] = jphi->at (iJ);
      jete[iJ] = je->at (iJ);

      missingPx += jetpt[iJ]*cos (jetphi[iJ]);
      missingPy += jetpt[iJ]*sin (jetphi[iJ]);
    }

    outTree->Fill ();
  }

  inFile->Close ();

  outTree->Write ();
  outFile->Close ();

}
