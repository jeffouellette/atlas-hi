#ifndef __CalcLeptonBkgs_C__
#define __CalcLeptonBkgs_C__

#include "Params.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

void CalcLeptonBkgs () {

  TChain* chain = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  double ttbarCounts = 0, ZtautauCounts = 0, allCounts = 0;

  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_pp_ttbar/*.root");
  cout << "Running over ttbar events..." << endl;
  ttbarCounts = chain->GetEntries ();
  chain->Reset ();


  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_pp_Ztautau/*.root");
  cout << "Running over ttbar events..." << endl;
  ZtautauCounts = chain->GetEntries ();
  chain->Reset ();


  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_nn_Zee/Group*.root");
  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_nn_Zmumu/Group*.root");
  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_np_Zee/Group*.root");
  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_np_Zmumu/Group*.root");
  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_pn_Zee/Group*.root");
  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_pn_Zmumu/Group*.root");
  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_pp_Zee/Group*.root");
  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_pp_Zmumu/Group*.root");
  cout << "Running over all events..." << endl;
  allCounts = chain->GetEntries ();
  chain->Reset ();


  const double allEffXS = allCounts * 642.820 / (415887 + 638346 + 638321 + 976912 + 357861 + 541661 + 561006 + 840470);
  const double ttbarEffXS = ttbarCounts * 56.695 * 0.54439 / 67734;
  const double ZtautauEffXS = ZtautauCounts * 642.820 / 193395;

  cout << "effective cross sections are" << endl;
  cout << "all Zs         = " << allEffXS << endl;
  cout << "ttbar bkg      = " << ttbarEffXS << endl;
  cout << "Z->tautau bkg  = " << ZtautauEffXS << endl;

  cout << "# ttbar events     = " << ttbarCounts << endl;
  cout << "# Z->tautau events = " << ZtautauCounts << endl;
  cout << "# events           = " << allCounts << endl;

  cout << "ttbar bkg rate     = " << (ttbarEffXS / allEffXS) * 100 << "%" << endl;
  cout << "Z->tautau bkg rate = " << (ZtautauEffXS / allEffXS) * 100 << "%" << endl;
}

#endif
