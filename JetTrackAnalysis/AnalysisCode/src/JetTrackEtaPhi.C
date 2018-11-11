#include "JetTrackEtaPhi.h"
#include "Params.h"

#include <TreeVariables.h>
#include <GlobalParams.h>
#include <Utils.h>
#include <Trigger.h>

#include <TH2D.h>
#include <TFile.h>

#include <iostream>

using namespace std;

using namespace atlashi;

namespace JetTrackAnalysis {

void JetTrackEtaPhi (const char* directory,
                     const int dataSet,
                     const bool isMC,
                     const bool isPeriodA,
                     const char* inFileName,
                     const double crossSection_microbarns,
                     const double filterEfficiency,
                     const int numberEvents)
{

  SetupDirectories ("", "JetTrackAnalysis/");

  const bool isSignalOnlySample = isMC && TString (inFileName).Contains ("signalonly");
  const TString identifier = GetIdentifier (dataSet, inFileName, isMC, isSignalOnlySample, isPeriodA);
  cout << "File Identifier: " << identifier << endl;

  /**** Find the relevant TTree for this run ****/
  TFile* file = GetFile (directory, dataSet, isMC, inFileName);
  TTree* tree = NULL;
  if (file) tree = (TTree*)file->Get ("bush");
  if (tree == NULL || file == NULL) {
   cout << "Error: In JetTrackEtaPhi.C: TTree not obtained for given data set. Quitting." << endl;
   return;
  }
  /**** End find TTree ****/

  TreeVariables* t = new TreeVariables (tree, false);
  if (crossSection_microbarns != 0)
   t->SetGetMCInfo (false, crossSection_microbarns, filterEfficiency, numberEvents);
  t->SetGetTracks ();
  t->SetGetHIJets ();
  t->SetGetSimpleJets ();
  t->SetGetFCals ();
  t->SetBranchAddresses ();
  //tree->SetBranchAddress ("akt4hi_njet", &t->jet_n);
  //tree->SetBranchAddress ("akt4hi_xcalib_jet_pt", &t->jet_pt);
  //tree->SetBranchAddress ("akt4hi_xcalib_jet_eta", &t->jet_eta);
  //tree->SetBranchAddress ("akt4hi_xcalib_jet_phi", &t->jet_phi);
  //tree->SetBranchAddress ("akt4hi_xcalib_jet_e", &t->jet_e);
  if (isMC) {
   tree->SetBranchAddress ("akt4_truth_njet", &t->truth_jet_n);
   tree->SetBranchAddress ("akt4_truth_jet_pt", &t->truth_jet_pt);
   tree->SetBranchAddress ("akt4_truth_jet_eta", &t->truth_jet_eta);
   tree->SetBranchAddress ("akt4_truth_jet_phi", &t->truth_jet_phi);
   tree->SetBranchAddress ("akt4_truth_jet_e", &t->truth_jet_e);
  }

  vector <Trigger*> triggers = {};
  if (!isMC) {
   for (int iTrig = 0; iTrig < nJetTrig; iTrig++) {
    Trigger* trig = new Trigger (jetTrigNames[iTrig].Data (), jetTrigMinPtCuts[iTrig], jetTrigMinEtaCuts[iTrig], jetTrigMaxEtaCuts[iTrig]);
    trig->minPt = jetTrigMinPtCuts[iTrig];
    trig->maxPt = jetTrigMaxPtCuts[iTrig];
    tree->SetBranchAddress (jetTrigNames[iTrig], &(trig->trigBool));
    tree->SetBranchAddress (jetTrigNames[iTrig]+"_prescale", &(trig->trigPrescale));
    triggers.push_back (trig);
   }
   //if (dataSet < 313629) {
   // for (int iTrig = 0; iTrig < nJetTrigIon; iTrig++) {
   //  Trigger* trig = new Trigger (jetTrigNamesIon[iTrig].Data (), jetTrigMinPtCutsIon[iTrig], jetTrigMinEtaCutsIon[iTrig], jetTrigMaxEtaCutsIon[iTrig]);
   //  trig->minPt = jetTrigMinPtCutsIon[iTrig];
   //  trig->maxPt = jetTrigMaxPtCutsIon[iTrig];
   //  tree->SetBranchAddress (jetTrigNamesIon[iTrig], &(trig->trigBool));
   //  tree->SetBranchAddress (jetTrigNamesIon[iTrig]+"_prescale", &(trig->trigPrescale));
   //  triggers.push_back (trig);
   // }
   //}
   //else {
   // for (int iTrig = 0; iTrig < nJetTrigPP; iTrig++) {
   //  Trigger* trig = new Trigger (jetTrigNamesPP[iTrig].Data (), jetTrigMinPtCutsPP[iTrig], jetTrigMinEtaCutsPP[iTrig], jetTrigMaxEtaCutsPP[iTrig]);
   //  trig->minPt = jetTrigMinPtCutsPP[iTrig];
   //  trig->maxPt = jetTrigMaxPtCutsPP[iTrig];
   //  tree->SetBranchAddress (jetTrigNamesPP[iTrig], &(trig->trigBool));
   //  tree->SetBranchAddress (jetTrigNamesPP[iTrig]+"_prescale", &(trig->trigPrescale));
   //  triggers.push_back (trig);
   // }
   //}
  }

  TH2D* JetTrackDeltaEtaDeltaPhi[numCentBins][2] = {};
  TH2D* JetPtSpectrum[numCentBins] = {};
  TH2D* JetCounts[numCentBins] = {};
  TH2D* TrackPtSpectrum[numCentBins] = {};
  TH2D* TrackCounts[numCentBins] = {};

  TH1D* TrackSelectionHist[numCentBins] = {};
  TH1D* JetSelectionHist[numCentBins] = {};

  TH2D* DijetDeltaEtaDeltaPhi[numCentBins] = {};
  TH2D* DijetEta1Eta2[numCentBins] = {};
  TH1D* DijetAj[numCentBins] = {};

  for (short iCent = 0; iCent < numCentBins; iCent++) {
   for (short iCut = 0; iCut < 2; iCut++) {
    JetTrackDeltaEtaDeltaPhi[iCent][iCut] = new TH2D (Form ("JetTrackDeltaEtaDeltaPhi_cent%i_cut%i", iCent, iCut), "", 80, 0, 8, numphibins, phibins);
    JetTrackDeltaEtaDeltaPhi[iCent][iCut]->Sumw2 ();
   }

   JetPtSpectrum[iCent] = new TH2D (Form ("JetPtSpectrum_cent%i", iCent), "", 28, 20, 140, numetabins, etabins);
   JetPtSpectrum[iCent]->Sumw2 ();

   JetCounts[iCent] = new TH2D (Form ("JetCounts_cent%i", iCent), "", 28, 20, 140, numetabins, etabins);
   JetCounts[iCent]->Sumw2 ();

   TrackPtSpectrum[iCent] = new TH2D (Form ("TrackPtSpectrum_cent%i", iCent), "", 99, 0.5, 50, 10, -2.5, 2.5);
   TrackPtSpectrum[iCent]->Sumw2 ();

   TrackCounts[iCent] = new TH2D (Form ("TrackCounts_cent%i", iCent), "", 99, 0.5, 50, 10, -2.5, 2.5);
   TrackCounts[iCent]->Sumw2 ();

   TrackSelectionHist[iCent] = new TH1D (Form ("TrackSelectionHist_cent%i", iCent), "", 10, -0.5, 9.5);
   TrackSelectionHist[iCent]->Sumw2 ();

   JetSelectionHist[iCent] = new TH1D (Form ("JetSelectionHist_cent%i", iCent), "", 11, -0.5, 10.5);
   JetSelectionHist[iCent]->Sumw2 ();

   DijetDeltaEtaDeltaPhi[iCent] = new TH2D (Form ("DijetDeltaEtaDeltaPhi_cent%i", iCent), "", 160, -8, 8, numphibins, phibins);
   DijetDeltaEtaDeltaPhi[iCent]->Sumw2 ();

   DijetEta1Eta2[iCent] = new TH2D (Form ("DijetEta1Eta2_cent%i", iCent), "", 98, -4.9, 4.9, 98, -4.9, 4.9);
   DijetEta1Eta2[iCent]->Sumw2 ();

   DijetAj[iCent] = new TH1D (Form ("DijetAj_cent%i", iCent), "", 100, 0, 1);
   DijetAj[iCent]->Sumw2 ();

  }

  const int nentries = tree->GetEntries ();

  for (int entry = 0; entry < nentries; entry++) {
   if (entry % 10000 == 0) cout << entry << endl;

   tree->GetEntry (entry);

   double fcal_et = (isPeriodA ? t->fcalA_et : t->fcalC_et);
   // find the centrality class
   int iCent = 0;
   while (iCent < numCentBins && fcal_et < centCuts[iCent]) iCent++;
   iCent--;
   if (iCent < 0 || iCent > numCentBins-1)
    continue;

   // find leading jet
   int lj = -1;
   int slj = -1;
   //int sslj = -1;
   for (int iJet = 0; iJet < t->jet_n; iJet++) {
    if (lj == -1 || t->jet_pt->at (lj) < t->jet_pt->at (iJet)) {
     lj = iJet;
    }
   }
   if (lj == -1)
    continue;
   JetSelectionHist[iCent]->Fill (0); // found leading jet

   // find sub-leading jet
   for (int iJet = 0; iJet < t->jet_n; iJet++) {
    if (iJet == lj || DeltaR (t->jet_eta->at (lj), t->jet_eta->at (iJet), t->jet_phi->at (lj), t->jet_phi->at (iJet)) < 0.8)
     continue;
    if (slj == -1 || t->jet_pt->at (slj) < t->jet_pt->at (iJet)) {
     slj = iJet;
    }
   }
   if (slj == -1)
    continue;
   JetSelectionHist[iCent]->Fill (1); // found subleading jet

   //// find sub-sub-leading jet
   //for (int iJet = 0; iJet < t->jet_n; iJet++) {
   // if (iJet == lj || iJet == slj)
   //  continue;
   // if (sslj == -1 || t->jet_pt->at (sslj) < t->jet_pt->at (iJet)) {
   //  sslj = iJet;
   // }
   //}

   // find whether triggered on leading jet
   float prescale = -1;
   if (!isMC) {
    for (Trigger* trig : triggers) {
     if (t->jet_pt->at (lj) < trig->minPt || trig->maxPt < t->jet_pt->at (lj) || t->jet_eta->at (lj) < trig->lowerEta || trig->upperEta < t->jet_eta->at (lj)|| !trig->trigBool)
      continue;
     else if (trig->trigPrescale > 0 && (prescale < 0 || trig->trigPrescale < prescale))
      prescale = trig->trigPrescale;
    }
    if (prescale < 0)
     continue; // skip events which the leading jet did not trigger
   }
   else {
    prescale = t->crossSection_microbarns / t->filterEfficiency / t->numberEvents;
   }
   JetSelectionHist[iCent]->Fill (2); // trigger fired

   for (int iJet = 0; iJet < t->jet_n; iJet++) {
    if (InDisabledHEC (t->jet_eta->at (iJet), t->jet_phi->at (iJet)))
     continue;
    JetPtSpectrum[iCent]->Fill (t->jet_pt->at (iJet), t->jet_eta->at (iJet), prescale);
    JetCounts[iCent]->Fill (t->jet_pt->at (iJet), t->jet_eta->at (iJet));
   }
   //JetPtSpectrum[iCent]->Fill (t->jet_pt->at (lj), t->jet_eta->at (lj), prescale);
   //JetCounts[iCent]->Fill (t->jet_pt->at (lj), t->jet_eta->at (lj));

   {
    float dphi = DeltaPhi (t->jet_phi->at (lj), t->jet_phi->at (slj), true);
    if (dphi < -pi/2.)
     dphi += 2*pi;
    DijetDeltaEtaDeltaPhi[iCent]->Fill (t->jet_eta->at (lj) - t->jet_eta->at (slj), dphi, prescale);
   }

   if (InDisabledHEC (t->jet_eta->at (lj), t->jet_phi->at (lj)))
    continue;
   JetSelectionHist[iCent]->Fill (3); // leading jet not in disabled HEC
   if (InDisabledHEC (t->jet_eta->at (slj), t->jet_phi->at (slj)))
    continue;
   JetSelectionHist[iCent]->Fill (4); // subleading jet not in disabled HEC

   if (DeltaPhi (t->jet_phi->at (lj), t->jet_phi->at (slj)) < 3.*pi / 4.)
    continue;
   JetSelectionHist[iCent]->Fill (5); // leading and subleading at least delta phi apart

   DijetEta1Eta2[iCent]->Fill (t->jet_eta->at (lj), t->jet_eta->at (slj), prescale);
   DijetAj[iCent]->Fill ((t->jet_pt->at (lj) - t->jet_pt->at (slj)) / (t->jet_pt->at (lj) + t->jet_pt->at (slj)), prescale);

   //if (sslj != -1 && trijetMaxPtRatio < t->jet_pt->at (sslj) / t->jet_pt->at (lj))
   // continue;
   JetSelectionHist[iCent]->Fill (6); // no significant 3rd jet

   // fill jet spectra & track correlations
   if (isMC) {
    double minDeltaR = 1000;
    for (int iTJet = 0; iTJet < t->truth_jet_n; iTJet++) {
     const double dR = DeltaR (t->jet_eta->at (lj), t->truth_jet_eta->at (iTJet), t->jet_phi->at (lj), t->truth_jet_phi->at (iTJet));
     if (dR < minDeltaR) {
      minDeltaR = dR;
     } 
    }
    if (minDeltaR > 0.2)
     continue;
   }

   if (t->jet_pt->at (lj) < jet_pt_cut) 
    continue; // jet pT cut
   JetSelectionHist[iCent]->Fill (7);

   if (isPeriodA && t->jet_eta->at (lj) > 3.2) 
    continue; // reject leading jets outside barrel to avoid a centrality bias
   else if (!isPeriodA && t->jet_eta->at (lj) < -3.2)
    continue;
   JetSelectionHist[iCent]->Fill (8);

   if (isPeriodA && t->jet_eta->at (slj) > 3.2) 
    continue; // reject subleading jets outside barrel to avoid a centrality bias
   else if (!isPeriodA && t->jet_eta->at (slj) < -3.2)
    continue;
   JetSelectionHist[iCent]->Fill (9);

   for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
    if (1e-3 * (t->trk_pt->at (iTrk)) < trk_pt_cut)
     continue; // track pT cut
    if (!t->trk_quality_4->at (iTrk))
     continue; // cut on track quality

    float dphi = DeltaPhi (t->jet_phi->at (lj), t->trk_phi->at (iTrk), true);
    if (dphi < -pi/2.)
     dphi += 2*pi;

    JetTrackDeltaEtaDeltaPhi[iCent][0]->Fill (abs (t->jet_eta->at (lj) - t->trk_eta->at (iTrk)), dphi, prescale);

    if (abs (t->trk_eta->at (iTrk) - t->jet_eta->at (slj)) < 2)
     continue; // reject tracks within delta eta < 2 of the opposing dijet

    JetTrackDeltaEtaDeltaPhi[iCent][1]->Fill (abs (t->jet_eta->at (lj) - t->trk_eta->at (iTrk)), dphi, prescale);
   }
   JetSelectionHist[iCent]->Fill (10, prescale);
   

   // fill track spectra
   for (int iTrk = 0; iTrk < t->ntrk; iTrk++) {
    if (!t->trk_quality_4->at (iTrk))
     continue; // cut on track quality
    TrackPtSpectrum[iCent]->Fill (1e-3 * (t->trk_pt->at (iTrk)), t->trk_eta->at (iTrk), prescale);
    TrackCounts[iCent]->Fill (1e-3 * (t->trk_pt->at (iTrk)), t->trk_eta->at (iTrk)); // need to convert MeV to GeV
   }
  }


  for (short iCent = 0; iCent < numCentBins; iCent++) {
   JetPtSpectrum[iCent]->Scale (1., "width");
   TrackPtSpectrum[iCent]->Scale (1., "width");
   //JetTrackDeltaEtaDeltaPhi[iCent][0]->Scale (1., "width");
   //JetTrackDeltaEtaDeltaPhi[iCent][1]->Scale (1., "width");
   DijetDeltaEtaDeltaPhi[iCent]->Scale (1., "width");
   DijetEta1Eta2[iCent]->Scale (1., "width");
   DijetAj[iCent]->Scale (1., "width");

   JetPtSpectrum[iCent]->GetXaxis ()->SetTitle ("#it{p}_{T}^{Jet} #left[GeV#right]");
   JetPtSpectrum[iCent]->GetYaxis ()->SetTitle ("#eta^{Jet}_{det}");
   JetPtSpectrum[iCent]->GetZaxis ()->SetTitle ("d^{3}N / d#it{p}_{T} d#eta d#left(Cent#right)");

   JetCounts[iCent]->GetXaxis ()->SetTitle ("#it{p}_{T}^{Jet} #left[GeV#right]");
   JetCounts[iCent]->GetYaxis ()->SetTitle ("#eta^{Jet}_{det}");
   JetCounts[iCent]->GetZaxis ()->SetTitle ("Counts");

   TrackPtSpectrum[iCent]->GetXaxis ()->SetTitle ("#it{p}_{T}^{Trk} #left[GeV#right]");
   TrackPtSpectrum[iCent]->GetYaxis ()->SetTitle ("#eta^{Trk}_{det}");
   TrackPtSpectrum[iCent]->GetZaxis ()->SetTitle ("d^{3}N / d#it{p}_{T} d#eta d#left(Cent#right)");

   TrackCounts[iCent]->GetXaxis ()->SetTitle ("#it{p}_{T}^{Trk} #left[GeV#right]");
   TrackCounts[iCent]->GetYaxis ()->SetTitle ("#eta^{Trk}_{det}");
   TrackCounts[iCent]->GetZaxis ()->SetTitle ("Counts");

   for (short iCut = 0; iCut < 2; iCut++) {
    JetTrackDeltaEtaDeltaPhi[iCent][iCut]->GetXaxis ()->SetTitle ("#left|#eta^{Jet}_{det} - #eta^{Trk}_{det}#right|");
    JetTrackDeltaEtaDeltaPhi[iCent][iCut]->GetXaxis ()->SetTitleOffset (1);
    JetTrackDeltaEtaDeltaPhi[iCent][iCut]->GetYaxis ()->SetTitle ("#phi^{Jet}_{det} - #phi^{Trk}_{det}");
    JetTrackDeltaEtaDeltaPhi[iCent][iCut]->GetYaxis ()->SetTitleOffset (1);
    JetTrackDeltaEtaDeltaPhi[iCent][iCut]->GetZaxis ()->SetTitle ("d^{3}N / d#Delta#eta d#Delta#phi d#left(Cent#right)");
    JetTrackDeltaEtaDeltaPhi[iCent][iCut]->GetZaxis ()->SetTitleOffset (1.4);
   }

   DijetDeltaEtaDeltaPhi[iCent]->GetXaxis ()->SetTitle ("#eta^{Lead Jet}_{det} - #eta^{Sublead. Jet}_{det}");
   DijetDeltaEtaDeltaPhi[iCent]->GetXaxis ()->SetTitleOffset (1);
   DijetDeltaEtaDeltaPhi[iCent]->GetYaxis ()->SetTitle ("#phi^{Lead Jet}_{det} - #phi^{Sublead. Jet}_{det}");
   DijetDeltaEtaDeltaPhi[iCent]->GetYaxis ()->SetTitleOffset (1); 
   DijetDeltaEtaDeltaPhi[iCent]->GetZaxis ()->SetTitle ("d^{3}N / d#Delta#eta d#Delta#phi d#left(Cent#right)");
   DijetDeltaEtaDeltaPhi[iCent]->GetZaxis ()->SetTitleOffset (1.4);

   DijetEta1Eta2[iCent]->GetXaxis ()->SetTitle ("Leading - Subleading #Delta#eta");
   DijetEta1Eta2[iCent]->GetXaxis ()->SetTitleOffset (1);
   DijetEta1Eta2[iCent]->GetYaxis ()->SetTitle ("Leading - Subleading #Delta#phi");
   DijetEta1Eta2[iCent]->GetYaxis ()->SetTitleOffset (1);
   DijetEta1Eta2[iCent]->GetZaxis ()->SetTitle ("d^{3}N / d#Delta#eta d#Delta#phi d#left(Cent#right)");
   DijetEta1Eta2[iCent]->GetZaxis ()->SetTitleOffset (1.4);

   DijetAj[iCent]->GetXaxis ()->SetTitle ("A_{J}");
   DijetAj[iCent]->GetXaxis ()->SetTitleOffset (1);
   DijetAj[iCent]->GetYaxis ()->SetTitle ("d^{2}N / dA_{J} d#left(Cent#right)");
   DijetAj[iCent]->GetYaxis ()->SetTitleOffset (1);
  }

  TFile* outFile = new TFile (Form ("%s/JetTrackEtaPhi/dataSet_%s.root", rootPath.Data (), identifier.Data ()), "recreate");

  for (short iCent = 0; iCent < numCentBins; iCent++) {
   JetPtSpectrum[iCent]->Write ();
   JetCounts[iCent]->Write ();
   TrackPtSpectrum[iCent]->Write ();
   TrackCounts[iCent]->Write ();

   JetSelectionHist[iCent]->Write ();
   TrackSelectionHist[iCent]->Write ();

   for (short iCut = 0; iCut < 2; iCut++) {
    JetTrackDeltaEtaDeltaPhi[iCent][iCut]->Write ();
   }

   DijetDeltaEtaDeltaPhi[iCent]->Write ();
   DijetEta1Eta2[iCent]->Write ();
   DijetAj[iCent]->Write ();
  }

  outFile->Close ();
}

} // end namespace
