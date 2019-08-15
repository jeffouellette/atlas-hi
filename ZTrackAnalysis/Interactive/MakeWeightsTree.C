/**
 * Takes in a file path and file pattern, and "adds" all the TTrees together in a TChain and slims it down for generating event-level weights.
 */
void MakeWeightsTree (const char* path, const char* filePattern, const char* outFileName = "eventWeightsTree.root") {

  TChain* inTree = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  inTree->Add (Form ("%s/%s", path, filePattern));

  inTree->SetBranchStatus ("event_number", 0);
  inTree->SetBranchStatus ("run_number", 0);
  inTree->SetBranchStatus ("isEE", 0);
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
  inTree->SetBranchStatus ("z_y", 0);
  inTree->SetBranchStatus ("z_phi", 0);
  inTree->SetBranchStatus ("z_m", 0);
  inTree->SetBranchStatus ("trk_pt", 0);
  inTree->SetBranchStatus ("trk_eta", 0);
  inTree->SetBranchStatus ("trk_phi", 0);
  inTree->SetBranchStatus ("trk_charge", 0);

  TFile* outFile = new TFile (Form ("%s/%s", path, outFileName), "recreate");
  outFile->Delete ("PbPbZTrackTree;*");
  outFile->Delete ("ppZTrackTree;*");

  outFile->cd ();
  TTree* outTree = inTree->CloneTree ();


  inTree = new TChain ("ppZTrackTree", "ppZTrackTree");
  inTree->Add (Form ("%s/%s", path, filePattern));

  inTree->SetBranchStatus ("event_number", 0);
  inTree->SetBranchStatus ("run_number", 0);
  inTree->SetBranchStatus ("isEE", 0);
  inTree->SetBranchStatus ("fcal_et", 0);
  inTree->SetBranchStatus ("q2", 0);
  inTree->SetBranchStatus ("psi2", 0);
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
  inTree->SetBranchStatus ("z_y", 0);
  inTree->SetBranchStatus ("z_phi", 0);
  inTree->SetBranchStatus ("z_m", 0);
  inTree->SetBranchStatus ("trk_pt", 0);
  inTree->SetBranchStatus ("trk_eta", 0);
  inTree->SetBranchStatus ("trk_phi", 0);
  inTree->SetBranchStatus ("trk_charge", 0);

  outTree = inTree->CloneTree ();

  outTree->SetDirectory (outFile);
  outFile->Close ();
}


void MakeDataWeightsTree () {
  MakeWeightsTree ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal", "outFile.root");
}


void MakeMinbiasWeightsTree () {
  MakeWeightsTree ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MinbiasAnalysis/Nominal", "3*.root");
}


void MakeMCWeightsTree () {
  MakeWeightsTree ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal", "*Z*.root");
}


//void Make2015HijingWeightsTree () {
//  MakeWeightsTree ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal", "PbPb_Hijing_15.root", "hijing2015_eventWeightsTree.root");
//}
//
//void Make2018HijingWeightsTree () {
//  MakeWeightsTree ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal", "PbPb_Hijing_18.root", "hijing2018_eventWeightsTree.root");
//}
