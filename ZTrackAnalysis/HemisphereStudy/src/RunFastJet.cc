#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH2.h>

#include <vector>

#include <Utilities.h>

using namespace fastjet;
using namespace std;
using namespace atlashi;


int main (int argc, char* argv[]) {

  if (argc-1 != 2) {
    cout << "Error! Need two arguments: inFile outFile" << endl;
    return 1;
  }

  TString _inFileName = argv[1];
  TString _outFileName = argv[2];

  TFile* inFile = new TFile (_inFileName, "read");
  TTree* inTree = (TTree*) inFile->Get ("tree");

  int code = 0, z_n = 0, part_n = 0, en = 0, rn = 0;
  vector<float>* in_z_pt = nullptr, *in_z_eta = nullptr, *in_z_phi = nullptr, *in_z_m = nullptr, *part_pt = nullptr, *part_eta = nullptr, *part_phi = nullptr;
  float z_pt = 0, z_eta = 0, z_y = 0, z_phi = 0, z_m = 0, missingPx = 0, missingPy = 0;
  int jet_n = 0;
  float jet_pt[10] = {};
  float jet_eta[10] = {};
  float jet_phi[10] = {};

  inTree->SetBranchAddress ("code", &code);
  inTree->SetBranchAddress ("z_n", &z_n);
  inTree->SetBranchAddress ("z_pt", &in_z_pt);
  inTree->SetBranchAddress ("z_eta", &in_z_eta);
  inTree->SetBranchAddress ("z_phi", &in_z_phi);
  inTree->SetBranchAddress ("z_m", &in_z_m);
  inTree->SetBranchAddress ("part_n", &part_n);
  inTree->SetBranchAddress ("part_pt", &part_pt);
  inTree->SetBranchAddress ("part_eta", &part_eta);
  inTree->SetBranchAddress ("part_phi", &part_phi);

  TFile* outFile = new TFile (_outFileName, "recreate");
  TTree* outTree = new TTree ("tree", "Zs and jets in pp Pythia8 MC");

  outTree->Branch ("code", &code);
  outTree->Branch ("run_number", &rn);
  outTree->Branch ("event_number", &en);
  outTree->Branch ("z_pt", &z_pt);
  outTree->Branch ("z_eta", &z_eta);
  outTree->Branch ("z_y", &z_y);
  outTree->Branch ("z_phi", &z_phi);
  outTree->Branch ("z_m", &z_m);
  outTree->Branch ("njt", &jet_n);
  outTree->Branch ("jtpt", jet_pt, "jtpt[njt]/F");
  outTree->Branch ("jteta", jet_eta, "jteta[njt]/F");
  outTree->Branch ("jtphi", jet_phi, "jtphi[njt]/F");
  outTree->Branch ("missingPx", &missingPx);
  outTree->Branch ("missingPy", &missingPy);

  outTree->SetDirectory (outFile); 

  vector <PseudoJet> particles;

  JetDefinition jet_defE (antikt_algorithm, 0.4, E_scheme);

  TH1D* h_njet_incl = new TH1D ("h_njet_incl", ";N_{jet};Events", 7, -0.5, 6.5);
  h_njet_incl->Sumw2 ();
  TH1D* h_njet_trig = new TH1D ("h_njet_trig", ";N_{jet};Events", 7, -0.5, 6.5);
  h_njet_trig->Sumw2 ();

  const int nEvts = inTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    inTree->GetEntry (iEvt);

    assert (z_n == 1);
    z_pt = in_z_pt->at (0);
    z_eta = in_z_eta->at (0);
    z_phi = in_z_phi->at (0);
    z_m = in_z_m->at (0);

    TLorentzVector tlv;
    tlv.SetPtEtaPhiM (z_pt, z_eta, z_phi, z_m);
    z_y = tlv.Rapidity ();

    missingPx = z_pt*cos (z_phi);
    missingPy = z_pt*sin (z_phi);

    particles.clear ();
    for (int iP = 0; iP < part_pt->size (); iP++) {
      tlv.SetPtEtaPhiM (part_pt->at (iP), part_eta->at (iP), part_phi->at (iP), 0);
      particles.push_back (PseudoJet (tlv.Px (), tlv.Py (), tlv.Pz (), tlv.E ()));
    }

    ClusterSequence csE (particles, jet_defE);

    vector <PseudoJet> jetsE = sorted_by_pt (csE.inclusive_jets ());

    jet_n = 0;
    for (PseudoJet jet : jetsE) {
      if (jet.pt () < 15)
        continue;

      jet_pt[jet_n] = jet.pt ();
      jet_eta[jet_n] = jet.eta ();
      jet_phi[jet_n] = jet.phi ();
      jet_n++;

      missingPx += jet.pt ()*cos (jet.phi ());
      missingPy += jet.pt ()*sin (jet.phi ());
    }

    h_njet_incl->Fill (jet_n);

    bool keepEvent = false;
    for (int iP = 0; iP < part_pt->size (); iP++) {
      if (part_pt->at (iP) > 20 && DeltaPhi (z_phi, part_phi->at (iP)) < pi/2)
        keepEvent = true;
    }
    if (!keepEvent)
      continue;

    h_njet_trig->Fill (jet_n);

    outTree->Fill ();
  }

  h_njet_incl->Scale (1./h_njet_incl->Integral ());
  h_njet_trig->Scale (1./h_njet_trig->Integral ());

  outTree->Write ();
  h_njet_incl->Write ();
  h_njet_trig->Write ();
  outFile->Close ();

  return 0; 
}
