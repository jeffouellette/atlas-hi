/**
 * Takes in a file path and file pattern, and "adds" all the TTrees together in a TChain and slims it down for generating event-level weights.
 */
void MakeZMixFile (const char* path, const char* filePattern) {

  TChain* inTree = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  inTree->Add (Form ("%s/%s", path, filePattern));

  inTree->SetBranchStatus ("event_number", 0);
  inTree->SetBranchStatus ("run_number", 0);
  inTree->SetBranchStatus ("event_weight", 0);
  //inTree->SetBranchStatus ("isEE", 0);
  //inTree->SetBranchStatus ("fcal_et", 0);
  inTree->SetBranchStatus ("q2", 0);
  inTree->SetBranchStatus ("psi2", 0);
  inTree->SetBranchStatus ("vz", 0);
  inTree->SetBranchStatus ("l1_pt", 0);
  inTree->SetBranchStatus ("l1_eta", 0);
  inTree->SetBranchStatus ("l1_phi", 0);
  inTree->SetBranchStatus ("l1_charge", 0);
  inTree->SetBranchStatus ("l1_trk_pt", 0);
  inTree->SetBranchStatus ("l1_trk_eta", 0);
  inTree->SetBranchStatus ("l1_trk_phi", 0);
  inTree->SetBranchStatus ("l2_pt", 0);
  inTree->SetBranchStatus ("l2_eta", 0);
  inTree->SetBranchStatus ("l2_phi", 0);
  inTree->SetBranchStatus ("l2_charge", 0);
  inTree->SetBranchStatus ("l2_trk_pt", 0);
  inTree->SetBranchStatus ("l2_trk_eta", 0);
  inTree->SetBranchStatus ("l2_trk_phi", 0);
  inTree->SetBranchStatus ("ntrk", 0);
  inTree->SetBranchStatus ("trk_pt", 0);
  inTree->SetBranchStatus ("trk_eta", 0);
  inTree->SetBranchStatus ("trk_phi", 0);
  inTree->SetBranchStatus ("trk_charge", 0);

  TFile* outFile = new TFile (Form ("%s/zMixFile.root", path), "recreate");
  outFile->Delete ("PbPbZTrackTree;*");
  outFile->Delete ("ppZTrackTree;*");

  outFile->cd ();
  TTree* outTree = inTree->CloneTree ();
  outTree->SetName ("PbPbZTree");
  outTree->SetDirectory (outFile);

  inTree = new TChain ("ppZTrackTree", "ppZTrackTree");
  inTree->Add (Form ("%s/%s", path, filePattern));

  inTree->SetBranchStatus ("event_number", 0);
  inTree->SetBranchStatus ("run_number", 0);
  inTree->SetBranchStatus ("event_weight", 0);
  //inTree->SetBranchStatus ("isEE", 0);
  inTree->SetBranchStatus ("fcal_et", 0);
  inTree->SetBranchStatus ("q2", 0);
  inTree->SetBranchStatus ("psi2", 0);
  inTree->SetBranchStatus ("vz", 0);
  inTree->SetBranchStatus ("l1_pt", 0);
  inTree->SetBranchStatus ("l1_eta", 0);
  inTree->SetBranchStatus ("l1_phi", 0);
  inTree->SetBranchStatus ("l1_charge", 0);
  inTree->SetBranchStatus ("l1_trk_pt", 0);
  inTree->SetBranchStatus ("l1_trk_eta", 0);
  inTree->SetBranchStatus ("l1_trk_phi", 0);
  inTree->SetBranchStatus ("l2_pt", 0);
  inTree->SetBranchStatus ("l2_eta", 0);
  inTree->SetBranchStatus ("l2_phi", 0);
  inTree->SetBranchStatus ("l2_charge", 0);
  inTree->SetBranchStatus ("l2_trk_pt", 0);
  inTree->SetBranchStatus ("l2_trk_eta", 0);
  inTree->SetBranchStatus ("l2_trk_phi", 0);
  inTree->SetBranchStatus ("ntrk", 0);
  inTree->SetBranchStatus ("trk_pt", 0);
  inTree->SetBranchStatus ("trk_eta", 0);
  inTree->SetBranchStatus ("trk_phi", 0);
  inTree->SetBranchStatus ("trk_charge", 0);

  outTree = inTree->CloneTree ();
  outTree->SetName ("ppZTree");

  outTree->SetDirectory (outFile);
  outFile->Write ();
  outFile->Close ();

}

void MakeDataZMixFile () {
  MakeZMixFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal", "outFile.root");
}

void MakeElectronLHMediumZMixFile () {
  MakeZMixFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Variations/ElectronLHMediumWPVariation", "outFile.root");
}

void MakeMuonTightZMixFile () {
  MakeZMixFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Variations/MuonTightWPVariation", "outFile.root");
}
