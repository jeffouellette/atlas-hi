#ifndef __CalcLeptonBkgs_C__
#define __CalcLeptonBkgs_C__

#include "Params.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <fstream>

void CalcLeptonBkgs () {

  TChain* chain = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  double ttbarEECounts = 0, ttbarMuMuCounts = 0, ZtautauEECounts = 0, ZtautauMuMuCounts = 0, eeAllCounts = 0, mumuAllCounts = 0;

  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_pp_ttbar/*.root");
  cout << "Running over ttbar events..." << endl;
  ttbarEECounts = chain->GetEntries ("isEE");
  ttbarMuMuCounts = chain->GetEntries ("!isEE");
  chain->Reset ();


  chain->Add ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18_pp_Ztautau/*.root");
  cout << "Running over ttbar events..." << endl;
  ZtautauEECounts = chain->GetEntries ("isEE");
  ZtautauMuMuCounts = chain->GetEntries ("!isEE");
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
  eeAllCounts = chain->GetEntries ("isEE");
  mumuAllCounts = chain->GetEntries ("!isEE");
  chain->Reset ();

  ofstream of;
  of.open ("DileptonBkgs.out");


  const double eeEffXS = eeAllCounts * 642.820 / (415887 + 638346 + 638321 + 976912);
  const double ttbarEEEffXS = ttbarEECounts * 56.695 * 0.54439 / 67734;
  const double ZtautauEEEffXS = ZtautauEECounts * 642.820 / 193395;

  const double mumuEffXS = mumuAllCounts * 642.820 / (357861 + 541661 + 561006 + 840470);
  const double ttbarMuMuEffXS = ttbarMuMuCounts * 56.695 * 0.54439 / 67734;
  const double ZtautauMuMuEffXS = ZtautauMuMuCounts * 642.820 / 193395;

  of << "electron effective cross sections are" << endl;
  of << "all Z->ee         = " << eeEffXS << endl;
  of << "ttbar->ee bkg      = " << ttbarEEEffXS << endl;
  of << "Z->tautau->ee bkg  = " << ZtautauEEEffXS << endl;
  of << endl;
  of << "# ttbar->ee events     = " << ttbarEECounts << endl;
  of << "# Z->tautau->ee events = " << ZtautauEECounts << endl;
  of << "# Z->ee events         = " << eeAllCounts << endl;
  of << endl;
  of << "ttbar->ee bkg rate     = " << (ttbarEEEffXS / eeEffXS) * 100 << "%" << endl;
  of << "Z->tautau->ee bkg rate = " << (ZtautauEEEffXS / eeEffXS) * 100 << "%" << endl;

  of << endl << endl;

  of << "muon effective cross sections are" << endl;
  of << "all Z->mumu         = " << mumuEffXS << endl;
  of << "ttbar->mumu bkg      = " << ttbarMuMuEffXS << endl;
  of << "Z->tautau->mumu bkg  = " << ZtautauMuMuEffXS << endl;
  of << endl;
  of << "# ttbar->mumu events     = " << ttbarMuMuCounts << endl;
  of << "# Z->tautau->mumu events = " << ZtautauMuMuCounts << endl;
  of << "# Z->mumu events         = " << mumuAllCounts << endl;
  of << endl;
  of << "ttbar->mumu bkg rate     = " << (ttbarMuMuEffXS / mumuEffXS) * 100 << "%" << endl;
  of << "Z->tautau->mumu bkg rate = " << (ZtautauMuMuEffXS / mumuEffXS) * 100 << "%" << endl;
}

#endif
