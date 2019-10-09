#ifndef __Systematic_h__
#define __Systematic_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>
#include <algorithm>
#include <map>

using namespace std;
using namespace atlashi;

class Systematic : public PhysicsAnalysis {

  protected:
  vector<PhysicsAnalysis*> variations;
  map <PhysicsAnalysis*, bool> variationDirs; // variation direction values are true for both, false for asymmetric
  vector<Systematic*> systematics;

  // map taking TH1Ds to TGAEs for systematics
  map <TH1D*, TGAE*> graphMap;

  // TFile with systematics graphs
  TFile* graphFile = nullptr;
  bool graphsLoaded = false;
  

  void NullifyErrors ();
  void AddGraphPair (TH1D* h);
  void CreateSysGraphs ();

  void WriteIAAs () override;

  public:

  string description;
  bool cancelIAA = true;

  Systematic (PhysicsAnalysis* nom, const char* _name, const char* _desc) : PhysicsAnalysis (){
    name = _name;
    description = _desc;
    CopyAnalysis (nom, true);
    CreateSysGraphs ();
    NullifyErrors ();
  }

  Systematic (PhysicsAnalysis* nom, const char* _name, const char* _desc, const char* _graphFileName) : PhysicsAnalysis () {
    name = _name;
    description = _desc;
    CopyAnalysis (nom, true);
    LoadGraphs (_graphFileName);
    NullifyErrors ();
  }

  vector<PhysicsAnalysis*>& GetVariations ()  { return variations;  }
  vector<Systematic*>&      GetSystematics () { return systematics; }

  virtual TGAE* GetTGAE (TH1D* h) override;

  //virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;
  virtual void LoadGraphs (const char* graphFileName);
  virtual void SaveGraphs (const char* graphFileName);

  virtual void AddVariation (PhysicsAnalysis* a, const bool applyBothWays = true);
  virtual void AddSystematic (Systematic* a);

  virtual void AddVariations (); // variations add linearly
  virtual void AddSystematics (); // systematics add in quadrature

  virtual void PrintRelSystematics (); // print relative systematics

  void PlotTrkYieldSystematics (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  void PlotTrkYieldSystematicsPtZ (const bool useTrkPt = true, const short pSpc = 2);
  void PlotSignalTrkYieldSystematics (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  void PlotSignalTrkYieldSystematicsPtZ (const bool useTrkPt = true, const short pSpc = 2);
  void PlotIAASystematics (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  void PlotIAASystematicsPtZ (const bool useTrkPt = true, const short pSpc = 2);

  virtual void PrintIAA (const bool printErrs, const bool useTrkPt = true, const short iCent = numCentBins-1, const short iPtZ = nPtZBins-1, const short iSpc = 2) override;

};




void Systematic :: AddGraphPair (TH1D* h) {
  if (!h)
    return;
  TGAE* g = make_graph (h);
  g->SetName (string (h->GetName ()).replace (0, 1, "g").c_str ());;
  graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
}




void Systematic :: CreateSysGraphs () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 2; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {

        //  AddGraphPair (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);

        //  AddGraphPair (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
        //  AddGraphPair (h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]);
        //  
        //  AddGraphPair (h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
        //  AddGraphPair (h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent]);
        //  //AddGraphPair (h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent]);

        //  AddGraphPair (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
        //  AddGraphPair (h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]);
        //  AddGraphPair (h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
        //  AddGraphPair (h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]);
        //  //AddGraphPair (h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]);
        //} // end loop over iPhi

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          AddGraphPair (h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]);
          AddGraphPair (h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent]);
        } // end loop over iPtTrk

        AddGraphPair (h_z_trk_zpt[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zpt_sub[iSpc][iPtZ][iCent]);
        //AddGraphPair (h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]);
        //AddGraphPair (h_z_trk_zpt_icp[iSpc][iPtZ][iCent]);

        AddGraphPair (h_z_trk_zxzh[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
        //AddGraphPair (h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
        //AddGraphPair (h_z_trk_zxzh_icp[iSpc][iPtZ][iCent]);
        
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Returns a TGraphAsymmErrors corresponding to this systematic
////////////////////////////////////////////////////////////////////////////////////////////////
TGAE* Systematic :: GetTGAE (TH1D* h) {
  if (graphMap[h])
    return graphMap[h];
  return nullptr;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
//void Systematic :: LoadHists (const char* histFileName, const bool _finishHists) {
//  PhysicsAnalysis :: LoadHists (histFileName, _finishHists);
//  CreateSysGraphs ();
//}
//void Systematic :: LoadHists (const char* histFileName) {
//  PhysicsAnalysis :: LoadHists (histFileName, false);
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
//        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
//          h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_phi_sub_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()));
//        }
//        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//        //  h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_pt_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//        //  h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_pt_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//        //  h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_xzh_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//        //  h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_xzh_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//        //}
//        h_z_trk_zpt_sub[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_zpt_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
//        h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_zpt_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
//        h_z_trk_zxzh_sub[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_zxzh_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
//        h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_zxzh_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
//      }
//    }
//  }
//}




//////////////////////////////////////////////////////////////////////////////////////////////////
//// Save histograms
//////////////////////////////////////////////////////////////////////////////////////////////////
//void Systematic :: SaveHists (const char* histFileName) {
//  SetupDirectories ("", "ZTrackAnalysis/");
//  if (!histsLoaded)
//    return;
//
//  PhysicsAnalysis :: SaveHists (histFileName);
//
//  TDirectory* _gDirectory = gDirectory;
//  if (!histFile) {
//    histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "update");
//    histFile->cd ();
//  }
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//
//      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
//        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
//          SafeWrite (h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]);
//        }
//        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//        //  SafeWrite (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
//        //  SafeWrite (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
//        //  SafeWrite (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
//        //}
//        SafeWrite (h_z_counts[iSpc][iPtZ][iCent]);
//      }
//    }
//  }
//  
//  histFile->Close ();
//  histFile = nullptr;
//  histsLoaded = false;
//
//  _gDirectory->cd ();
//  return;
//}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: LoadGraphs (const char* graphFileName) {
  SetupDirectories ("", "ZTrackAnalysis/");

  TDirectory* _gDirectory = gDirectory;
  graphFile = new TFile (Form ("%s/%s", rootPath.Data (), graphFileName), "read");

  TH1D* h = nullptr;
  TGAE* g = nullptr;

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          h = h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent];
          g = (TGAE*) graphFile->Get (Form ("g_z_trk_phi_sub_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()));
          graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        }
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  h = h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
        //  g = (TGAE*) graphFile->Get (Form ("g_z_trk_pt_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        //  graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        //  //h = (h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
        //  //g = (TGAE*) graphFile->Get (Form ("g_z_trk_pt_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        //  //graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        //  h = h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];
        //  g = (TGAE*) graphFile->Get (Form ("g_z_trk_xzh_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        //  graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        //  //h = (h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
        //  //g = (TGAE*) graphFile->Get (Form ("g_z_trk_xzh_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        //  //graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        //}
        h = h_z_trk_zpt_sub[iSpc][iPtZ][iCent];
        g = (TGAE*) graphFile->Get (Form ("g_z_trk_zpt_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        //h = h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent];
        //g = (TGAE*) graphFile->Get (Form ("g_z_trk_zpt_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        //graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        h = h_z_trk_zxzh_sub[iSpc][iPtZ][iCent];
        g = (TGAE*) graphFile->Get (Form ("g_z_trk_zxzh_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        //h = h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent];
        //g = (TGAE*) graphFile->Get (Form ("g_z_trk_zxzh_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        //graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
      }
    }
  }

  _gDirectory->cd ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: SaveGraphs (const char* graphFileName) {
  SetupDirectories ("", "ZTrackAnalysis/");

  TDirectory* _gDirectory = gDirectory;
  graphFile = new TFile (Form ("%s/%s", rootPath.Data (), graphFileName), "recreate");

  for (map <TH1D*, TGAE*> :: iterator element = graphMap.begin (); element != graphMap.end (); ++element) {
    TGAE* g = element->second;
    //TH1D* h = element->first;

    //h->Write ();
    g->Write ();
  }

  graphFile->Close ();

  _gDirectory->cd ();

  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Adds a variation to consider
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: AddVariation (PhysicsAnalysis* a, const bool applyBothWays) {
  if (find (variations.begin (), variations.end (), a) == variations.end ()) {
    variations.push_back (a);
    variationDirs.insert (pair <PhysicsAnalysis*, bool> (a, applyBothWays));
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Adds a systematic error set to consider
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: AddSystematic (Systematic* a) {
  if (find (systematics.begin (), systematics.end (), a) == systematics.end ())
    systematics.push_back (a);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Sets the errors in this systematic as a combination of all added variations.
// Takes the maximum error for each point.
// Intended for combining up & down variations, but expandable for additional categories
// (e.g. track quality criteria)
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: AddVariations () {

  //NullifyErrors ();

  for (PhysicsAnalysis* a : variations) {

    cout << "Adding variation " << a->Name () << " to systematic " << name << endl;

    a->SubtractBackground ();
    a->CalculateIAA ();
    //a->CalculateICP ();

    TGAE* sys = nullptr;
    TH1D* var = nullptr;

    for (short iSpc = 2; iSpc < 3; iSpc++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

        //// Hadron yield systematics, signal & signal+bkg levels
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  for (short iCent = 0; iCent < numCentBins; iCent++) {
        //    sys = GetTGAE (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
        //    var = a->h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent];
        //    if (sys && var) CalcSystematics (sys, var, variationDirs[a]);

        //    sys = GetTGAE (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
        //    var = a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
        //    if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        //    sys = GetTGAE (h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]);
        //    var = a->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
        //    if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        //    sys = GetTGAE (h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
        //    var = a->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
        //    if (sys && var) CalcSystematics (sys, var, variationDirs[a]);

        //    sys = GetTGAE (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
        //    var = a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
        //    if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        //    sys = GetTGAE (h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]);
        //    var = a->h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];
        //    if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        //    sys = GetTGAE (h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
        //    var = a->h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
        //    if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        //  } // end loop over cents
        //} // end loop over phi

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            sys = GetTGAE (h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]);
            var = a->h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
            sys = GetTGAE (h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent]);
            var = a->h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          } // end loop over cents
        } // end loop over iPtTrk

        for (short iCent = 0; iCent < numCentBins; iCent++) {
          sys = GetTGAE (h_z_trk_zpt[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zpt[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          sys = GetTGAE (h_z_trk_zpt_sub[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zpt_sub[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          //sys = GetTGAE (h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]);
          //var = a->h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent];
          //if (sys && var) CalcSystematics (sys, var, variationDirs[a]);

          sys = GetTGAE (h_z_trk_zxzh[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zxzh[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          sys = GetTGAE (h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zxzh_sub[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          //sys = GetTGAE (h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);
          //var = a->h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent];
          //if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        } // end loop over cents


        //// IAA, ICP systematics
        //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        //  for (short iCent = 1; iCent < numCentBins; iCent++) {
        //    sys = GetTGAE (h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent]);
        //    var = a->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];
        //    if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        //    sys = GetTGAE (h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]);
        //    var = a->h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent];
        //    if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        //  } // end loop over cents

        //  //for (short iCent = 2; iCent < numCentBins; iCent++) {
        //  //  sys = GetTGAE (h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent]);
        //  //  var = a->h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent];
        //  //  if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        //  //  sys = GetTGAE (h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]);
        //  //  var = a->h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent];
        //  //  if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        //  //} // end loop over cents
        //} // end loop over phi

        if (cancelIAA) {
          for (short iCent = 1; iCent < numCentBins; iCent++) {
            sys = GetTGAE (h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]);
            var = a->h_z_trk_zpt_iaa[iSpc][iPtZ][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);

            sys = GetTGAE (h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
            var = a->h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          }
        }

        //for (short iCent = 1; iCent < numCentBins; iCent++) {
        //  sys = GetTGAE (h_z_trk_zpt_icp[iSpc][iPtZ][iCent]);
        //  var = a->h_z_trk_zpt_icp[iSpc][iPtZ][iCent];
        //  if (sys && var) CalcSystematics (sys, var, variationDirs[a]);

        //  sys = GetTGAE (h_z_trk_zxzh_icp[iSpc][iPtZ][iCent]);
        //  var = a->h_z_trk_zxzh_icp[iSpc][iPtZ][iCent];
        //  if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        //} // end loop over cents

      } // end loop over pT^Z bins
    } // end loop over species

  } // end loop over variations


  if (!cancelIAA) {
    TGAE* pp_sys, *PbPb_sys, *iaa_sys;
    for (short iSpc = 2; iSpc < 3; iSpc++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
        pp_sys = GetTGAE (h_z_trk_zpt_sub[iSpc][iPtZ][0]);
        for (short iCent = 1; iCent < numCentBins; iCent++) {
          PbPb_sys = GetTGAE (h_z_trk_zpt_sub[iSpc][iPtZ][iCent]);

          iaa_sys = GetTGAE (h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]);

          for (int ix = 0; ix < iaa_sys->GetN (); ix++) {
            double x, y, yhierr, yloerr;

            pp_sys->GetPoint (ix, x, y);
            yhierr = pp_sys->GetErrorYhigh (ix);
            yloerr = pp_sys->GetErrorYlow (ix);
            const double pp_relsys_hi = yhierr / y;
            const double pp_relsys_lo = yloerr / y;

            //if (name == "trkEffSys")
            //  cout << pp_relsys_hi << ", " << pp_relsys_lo << endl;

            PbPb_sys->GetPoint (ix, x, y);
            yhierr = PbPb_sys->GetErrorYhigh (ix);
            yloerr = PbPb_sys->GetErrorYlow (ix);
            const double PbPb_relsys_hi = yhierr / y;
            const double PbPb_relsys_lo = yloerr / y;

            //if (name == "trkEffSys")
            //  cout << PbPb_relsys_hi << ", " << PbPb_relsys_lo << endl;

            const double iaa_relsys_hi = TMath::Sqrt (pow (PbPb_relsys_hi, 2) + pow (pp_relsys_hi, 2));
            const double iaa_relsys_lo = TMath::Sqrt (pow (PbPb_relsys_lo, 2) + pow (pp_relsys_lo, 2));

            iaa_sys->GetPoint (ix, x, y);
            iaa_sys->SetPointEYhigh (ix, iaa_relsys_hi * y);
            iaa_sys->SetPointEYlow (ix, iaa_relsys_lo * y);
          }
        } // end loop over cents

        pp_sys = GetTGAE (h_z_trk_zxzh_sub[iSpc][iPtZ][0]);
        for (short iCent = 1; iCent < numCentBins; iCent++) {
          PbPb_sys = GetTGAE (h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);

          iaa_sys = GetTGAE (h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);

          for (int ix = 0; ix < iaa_sys->GetN (); ix++) {
            double x, y, yhierr, yloerr;

            pp_sys->GetPoint (ix, x, y);
            yhierr = pp_sys->GetErrorYhigh (ix);
            yloerr = pp_sys->GetErrorYlow (ix);
            const double pp_relsys_hi = yhierr / y;
            const double pp_relsys_lo = yloerr / y;

            PbPb_sys->GetPoint (ix, x, y);
            yhierr = PbPb_sys->GetErrorYhigh (ix);
            yloerr = PbPb_sys->GetErrorYlow (ix);
            const double PbPb_relsys_hi = yhierr / y;
            const double PbPb_relsys_lo = yloerr / y;

            const double iaa_relsys_hi = TMath::Sqrt (pow (PbPb_relsys_hi, 2) + pow (pp_relsys_hi, 2));
            const double iaa_relsys_lo = TMath::Sqrt (pow (PbPb_relsys_lo, 2) + pow (pp_relsys_lo, 2));

            iaa_sys->GetPoint (ix, x, y);
            iaa_sys->SetPointEYhigh (ix, iaa_relsys_hi * y);
            iaa_sys->SetPointEYlow (ix, iaa_relsys_lo * y);
          }

        } // end loop over cents
      } // end loop over pT^Z bins
    } // end loop over species
  }

  variations.clear ();
  variationDirs.clear ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Addition of independent systematics in quadrature; adds the errors of s to this systematic.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: NullifyErrors () {
  for (map <TH1D*, TGAE*> :: iterator element = graphMap.begin (); element != graphMap.end (); ++element) {
    if (element->second) ResetTGAEErrors (element->second);
    if (element->first) ResetHistErrors (element->first);
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Addition of independent systematics in quadrature; adds the errors of s to this systematic.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: AddSystematics () {

  NullifyErrors ();

  for (Systematic* s : systematics) {

    TGAE* master = nullptr, *sys = nullptr;

    for (short iSpc = 2; iSpc < 3; iSpc++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

        //// Hadron yield systematics, signal & signal+bkg levels
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  for (short iCent = 0; iCent < numCentBins; iCent++) {
        //    master = GetTGAE (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
        //    sys = s->GetTGAE (s->h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
        //    if (master && sys) AddErrorsInQuadrature (master, sys);

        //    master = GetTGAE (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
        //    sys = s->GetTGAE (s->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
        //    if (master && sys) AddErrorsInQuadrature (master, sys);
        //    master = GetTGAE (h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]);
        //    sys = s->GetTGAE (s->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]);
        //    if (master && sys) AddErrorsInQuadrature (master, sys);
        //    master = GetTGAE (h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
        //    sys = s->GetTGAE (s->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
        //    if (master && sys) AddErrorsInQuadrature (master, sys);

        //    master = GetTGAE (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
        //    sys = s->GetTGAE (s->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
        //    if (master && sys) AddErrorsInQuadrature (master, sys);
        //    master = GetTGAE (h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]);
        //    sys = s->GetTGAE (s->h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]);
        //    if (master && sys) AddErrorsInQuadrature (master, sys);
        //    master = GetTGAE (h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
        //    sys = s->GetTGAE (s->h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
        //    if (master && sys) AddErrorsInQuadrature (master, sys);
        //  } // end loop over cents
        //} // end loop over phi

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            master = GetTGAE (h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]);
            sys = s->GetTGAE (s->h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);

            master = GetTGAE (h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent]);
            sys = s->GetTGAE (s->h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);
          } // end loop over cents
        } // end loop over iPtTrk

        for (short iCent = 0; iCent < numCentBins; iCent++) {
          master = GetTGAE (h_z_trk_zpt[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zpt[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);
          master = GetTGAE (h_z_trk_zpt_sub[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zpt_sub[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);
          //master = GetTGAE (h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]);
          //sys = s->GetTGAE (s->h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]);
          //if (master && sys) AddErrorsInQuadrature (master, sys);

          master = GetTGAE (h_z_trk_zxzh[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zxzh[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);
          master = GetTGAE (h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);
          //master = GetTGAE (h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);
          //sys = s->GetTGAE (s->h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);
          //if (master && sys) AddErrorsInQuadrature (master, sys);
        } // end loop over cents


        //// IAA, ICP systematics
        //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        //  for (short iCent = 1; iCent < numCentBins; iCent++) {
        //    master = GetTGAE (h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent]);
        //    sys = s->GetTGAE (s->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent]);
        //    if (master && sys) AddErrorsInQuadrature (master, sys);
        //    master = GetTGAE (h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]);
        //    sys = s->GetTGAE (s->h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]);
        //    if (master && sys) AddErrorsInQuadrature (master, sys);
        //  } // end loop over cents
        //  //for (short iCent = 2; iCent < numCentBins; iCent++) {
        //  //  master = GetTGAE (h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent]);
        //  //  sys = s->GetTGAE (s->h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent]);
        //  //  if (master && sys) AddErrorsInQuadrature (master, sys);
        //  //  master = GetTGAE (h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]);
        //  //  sys = s->GetTGAE (s->h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]);
        //  //  if (master && sys) AddErrorsInQuadrature (master, sys);
        //  //} // end loop over cents
        //} // end loop over phi

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          master = GetTGAE (h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);

          master = GetTGAE (h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);
        }

        //for (short iCent = 1; iCent < numCentBins; iCent++) {
        //  master = GetTGAE (h_z_trk_zpt_icp[iSpc][iPtZ][iCent]);
        //  sys = s->GetTGAE (s->h_z_trk_zpt_icp[iSpc][iPtZ][iCent]);
        //  if (master && sys) AddErrorsInQuadrature (master, sys);

        //  master = GetTGAE (h_z_trk_zxzh_icp[iSpc][iPtZ][iCent]);
        //  sys = s->GetTGAE (s->h_z_trk_zxzh_icp[iSpc][iPtZ][iCent]);
        //  if (master && sys) AddErrorsInQuadrature (master, sys);
        //} // end loop over cents
      } // end loop over pT^Z bins
    } // end loop over species

  } // end loop over systematics
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Prints the maximum relative systematics in this bin.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: PrintRelSystematics () {
  TH1D* centralVals = nullptr, *h = nullptr, *highs = nullptr, *lows = nullptr;

  cout << endl;
  for (Systematic* s : systematics) {
    cout << s->description;
    for (short iCent : {0, numCentBins-1}) {
      for (short iPtZ : {nPtZBins-2, nPtZBins-1}) {
        cout << " & ";
        float _maxRelSys = 0;
        for (bool useTrkPt : {true}) {
          centralVals = (useTrkPt ? h_z_trk_zpt_sub[2][iPtZ][iCent] : h_z_trk_zxzh_sub[2][iPtZ][iCent]);
          h = (useTrkPt ? s->h_z_trk_zpt_sub[2][iPtZ][iCent] : s->h_z_trk_zxzh_sub[2][iPtZ][iCent]);

          highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
          if (!centralVals) continue;

          SaveRelativeErrors (s->GetTGAE (h), GetTGAE (centralVals), highs, lows);

          highs->Scale (100.);
          lows->Scale (100.);

          _maxRelSys = fmax (fabs (_maxRelSys), fabs (highs->GetMaximum ()));
          _maxRelSys = fmax (fabs (_maxRelSys), fabs (lows->GetMaximum ()));

          SaferDelete (highs);
          SaferDelete (lows);
        } // end loop over useTrkPt

        int i = 0;
        while (_maxRelSys < 1 && i < 3) {
          _maxRelSys *= 10;
          i++;
        }
        _maxRelSys = floor (_maxRelSys);
        _maxRelSys = _maxRelSys / pow (10, i);

        //if (_maxRelSys < 0.1)
        //  cout << "< 0.1\\\%";
        //else
          cout << _maxRelSys << "\\\%";
      } // end loop over iPtZ
    } // end loop iCent
    cout << "\\\\ \\hline" << endl;
  } // end loop over systematics
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot this set of systematics
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: PlotTrkYieldSystematics (const short pSpc, const short pPtZ) {
  const char* canvasName = Form ("c_TrkYieldSys");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
    gPad->SetLogx ();
  }
  c->cd ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {

          TH1D* centralVals = nullptr, *highs = nullptr, *lows = nullptr;
          centralVals = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];

          if (!centralVals) continue;

          highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
          SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

          highs->Scale (100.);
          lows->Scale (100.);

          highs->GetXaxis ()->SetMoreLogLabels ();
          highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
          highs->GetYaxis ()->SetTitle ("#deltaY / Y (#it{p}_{T}^{ch}) [%]");

          highs->SetLineColor (kGray+1);
          highs->SetLineStyle (1);
          highs->SetLineWidth (3);

          highs->DrawCopy ("][ hist");

          if (systematics.size () == 0)
            myLineText (0.65, 0.88, kGray+1, 1, description.c_str (), 1, 0.04);
          else
            myLineText (0.65, 0.92, kGray+1, 1, "Total", 1, 0.026);

          lows->GetXaxis ()->SetMoreLogLabels ();
          lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
          lows->GetYaxis ()->SetTitle ("#deltaY / Y (#it{p}_{T}^{ch}) [%]");

          lows->SetLineColor (kGray+1);
          lows->SetLineStyle (1);
          lows->SetLineWidth (3);

          lows->DrawCopy ("][ same hist");

          delete highs, lows;

          short iSys = 0;
          for (Systematic* sys : systematics) {

            if (sys->Name () == string ("bkgStatSys") || sys->Name () == string ("bkgMixSys"))
              continue;

            TH1D* h = sys->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
            highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
            lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
            
            SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

            highs->Scale (100.);
            lows->Scale (100.);

            highs->GetXaxis ()->SetMoreLogLabels ();
            highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            highs->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
            highs->GetYaxis ()->SetTitle ("#deltaY / Y (#it{p}_{T}^{ch}) [%]");

            highs->SetLineColor (colors[iSys+1]);
            highs->SetLineStyle (iSys+2);

            highs->DrawCopy ("][ hist same");

            lows->GetXaxis ()->SetMoreLogLabels ();
            lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            lows->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
            lows->GetYaxis ()->SetTitle ("#deltaY / Y (#it{p}_{T}^{ch}) [%]");

            lows->SetLineColor (colors[iSys+1]);
            lows->SetLineStyle (iSys+2);

            lows->DrawCopy ("][ same hist");

            myLineText (0.65, 0.89-0.026*iSys, colors[iSys+1], iSys+2, sys->description.c_str (), 1, 0.026);

            delete highs, lows;
            iSys++;
          }

          myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
          if (iCent == 0)
            myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
          else {
            myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
          }

          if (iPtZ == nPtZBins-1)
            myText (0.24, 0.86, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04);
          else
            myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

          const char* lo = GetPiString (phiLowBins[iPhi]);
          const char* hi = GetPiString (phiHighBins[iPhi]);
          myText (0.24, 0.80, kBlack, Form ("%s < #left|#Delta#phi#right| < %s", lo, hi), 0.04);

          c->SaveAs (Form ("%s/TrkYieldSystematics/%s_iPtZ%i_iPhi%i_iCent%i.pdf", plotPath.Data (), spc, iPtZ, iPhi, iCent));

        } // end loop over centralities

      } // end loop over Phi bins

    } // end loop over PtZ
  } // end loop over species
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot this set of systematics
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: PlotTrkYieldSystematicsPtZ (const bool useTrkPt, const short pSpc) {
  const char* canvasName = Form ("c_TrkYieldSys");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
    gPad->SetLogx ();
  }
  c->cd ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {

        TH1D* centralVals = nullptr, *highs = nullptr, *lows = nullptr;
        centralVals = (useTrkPt ? h_z_trk_zpt[iSpc][iPtZ][iCent] : h_z_trk_zxzh[iSpc][iPtZ][iCent]);

        if (!centralVals) continue;

        highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
        lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
        SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

        highs->Scale (100.);
        lows->Scale (100.);

        highs->GetXaxis ()->SetMoreLogLabels ();
        highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        highs->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
        highs->GetYaxis ()->SetTitle (useTrkPt ? "#deltaY / Y (#it{p}_{T}^{ch}) [%]" : "#deltaY / Y (#it{x}_{hZ}) [%]");

        highs->SetLineColor (kGray+1);
        highs->SetLineStyle (1);
        highs->SetLineWidth (3);

        highs->DrawCopy ("][ hist");

        if (systematics.size () == 0)
          myLineText (0.65, 0.88, kGray+1, 1, description.c_str (), 1, 0.04);
        else
          myLineText (0.65, 0.92, kGray+1, 1, "Total", 1, 0.026);

        lows->GetXaxis ()->SetMoreLogLabels ();
        lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        lows->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
        lows->GetYaxis ()->SetTitle (useTrkPt ? "#deltaY / Y (#it{p}_{T}^{ch}) [%]" : "#deltaY / Y (#it{x}_{hZ}) [%]");

        lows->SetLineColor (kGray+1);
        lows->SetLineStyle (1);
        lows->SetLineWidth (3);

        lows->DrawCopy ("][ same hist");

        delete highs, lows;

        short iSys = 0;
        for (Systematic* sys : systematics) {

          if (sys->Name () == string ("bkgStatSys") || sys->Name () == string ("bkgMixSys"))
            continue;

          TH1D* h = (useTrkPt ? sys->h_z_trk_zpt[iSpc][iPtZ][iCent] : sys->h_z_trk_zxzh[iSpc][iPtZ][iCent]);
          highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
          
          SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

          highs->Scale (100.);
          lows->Scale (100.);

          highs->GetXaxis ()->SetMoreLogLabels ();
          highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
          highs->GetYaxis ()->SetTitle (useTrkPt ? "#deltaY / Y (#it{p}_{T}^{ch}) [%]" : "#deltaY / Y (#it{x}_{hZ}) [%]");

          highs->SetLineColor (colors[iSys+1]);
          highs->SetLineStyle (iSys+2);

          highs->DrawCopy ("][ hist same");

          lows->GetXaxis ()->SetMoreLogLabels ();
          lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
          lows->GetYaxis ()->SetTitle (useTrkPt ? "#deltaY / Y (#it{p}_{T}^{ch}) [%]" : "#deltaY / Y (#it{x}_{hZ}) [%]");

          lows->SetLineColor (colors[iSys+1]);
          lows->SetLineStyle (iSys+2);

          lows->DrawCopy ("][ same hist");

          myLineText (0.65, 0.89-0.026*iSys, colors[iSys+1], iSys+2, sys->description.c_str (), 1, 0.026);

          delete highs, lows;
          iSys++;
        }

        myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
        if (iCent == 0)
          myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
        else {
          myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
        }

        if (iPtZ == nPtZBins-1)
          myText (0.24, 0.86, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04);
        else
          myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

        const char* lo = GetPiString (phiLowBins[1]);
        const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
        myText (0.24, 0.80, kBlack, Form ("%s < #left|#Delta#phi#right| < %s", lo, hi), 0.04);

        c->SaveAs (Form ("%s/TrkYieldSystematics/%s_%s_iPtZ%i_iCent%i.pdf", plotPath.Data (), spc, useTrkPt ? "pTTrk":"xHZ", iPtZ, iCent));

      } // end loop over centralities

    } // end loop over PtZ
  } // end loop over species
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot this set of systematics on the signal track yield
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: PlotSignalTrkYieldSystematics (const short pSpc, const short pPtZ) {
  const char* canvasName = Form ("c_SignalTrkYieldSysPtZ");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
    gPad->SetLogx ();
  }
  c->cd ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {

          gPad->Clear ();

          TH1D* centralVals = nullptr, *highs = nullptr, *lows = nullptr;
          centralVals = h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];

          if (!centralVals) continue;

          highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
          SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

          highs->Scale (100.);
          lows->Scale (100.);

          highs->GetXaxis ()->SetMoreLogLabels ();
          highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
          highs->GetYaxis ()->SetTitle ("#deltaY_{sub} / Y_{sub} (#it{p}_{T}^{ch}) [%]");

          highs->SetLineColor (kGray+1);
          highs->SetLineStyle (1);
          highs->SetLineWidth (3);

          highs->DrawCopy ("][ hist");

          if (systematics.size () == 0)
            myLineText (0.65, 0.88, kGray+1, 1, description.c_str (), 1, 0.04);
          else
            myLineText (0.65, 0.92, kGray+1, 1, "Total", 1, 0.026);

          lows->GetXaxis ()->SetMoreLogLabels ();
          lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
          lows->GetYaxis ()->SetTitle ("Y_{sub} (#it{p}_{T}^{ch}) Relative error");
          lows->GetYaxis ()->SetTitle ("#deltaY_{sub} / Y_{sub} (#it{p}_{T}^{ch}) [%]");

          lows->SetLineColor (kGray+1);
          lows->SetLineStyle (1);
          lows->SetLineWidth (3);

          lows->DrawCopy ("][ same hist");

          delete highs, lows;

          short iSys = 0;
          for (Systematic* sys : systematics) {

            TH1D* h = sys->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
            highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
            lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
            
            SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

            highs->Scale (100.);
            lows->Scale (100.);

            highs->GetXaxis ()->SetMoreLogLabels ();
            highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            highs->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
            highs->GetYaxis ()->SetTitle ("#deltaY_{sub} / Y_{sub} (#it{p}_{T}^{ch}) [%]");

            highs->SetLineColor (colors[iSys+1]);
            highs->SetLineStyle (iSys+2);

            highs->DrawCopy ("][ hist same");

            lows->GetXaxis ()->SetMoreLogLabels ();
            lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            lows->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
            lows->GetYaxis ()->SetTitle ("#deltaY_{sub} / Y_{sub} (#it{p}_{T}^{ch}) [%]");

            lows->SetLineColor (colors[iSys+1]);
            lows->SetLineStyle (iSys+2);

            lows->DrawCopy ("][ same hist");

            myLineText (0.65, 0.89-0.026*iSys, colors[iSys+1], iSys+2, sys->description.c_str (), 1, 0.026);

            delete highs, lows;
            iSys++;
          }

          myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
          if (iCent == 0)
            myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
          else {
            myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
          }

          if (iPtZ == nPtZBins-1)
            myText (0.24, 0.86, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04);
          else
            myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

          const char* lo = GetPiString (phiLowBins[iPhi]);
          const char* hi = GetPiString (phiHighBins[iPhi]);
          myText (0.24, 0.80, kBlack, Form ("%s < #left|#Delta#phi#right| < %s", lo, hi), 0.04);

          c->SaveAs (Form ("%s/TrkYieldSignalSystematics/%s_iPtZ%i_iPhi%i_iCent%i.pdf", plotPath.Data (), spc, iPtZ, iPhi, iCent));

        } // end loop over centralities

      } // end loop over Phi bins

    } // end loop over PtZ
  } // end loop over species
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot this set of systematics on the signal track yield
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: PlotSignalTrkYieldSystematicsPtZ (const bool useTrkPt, const short pSpc) {
  const char* canvasName = Form ("c_SignalTrkYieldSys");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
    gPad->SetLogx ();
  }
  c->cd ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {

        gPad->Clear ();

        TH1D* centralVals = nullptr, *highs = nullptr, *lows = nullptr;
        centralVals = (useTrkPt ? h_z_trk_zpt_sub[iSpc][iPtZ][iCent] : h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);

        if (!centralVals) continue;

        highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
        lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
        SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

        highs->Scale (100.);
        lows->Scale (100.);

        highs->GetXaxis ()->SetMoreLogLabels ();
        highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        highs->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
        highs->GetYaxis ()->SetTitle (useTrkPt ? "#deltaY_{sub} / Y_{sub} (#it{p}_{T}^{ch}) [%]" : "#deltaY_{sub} / Y_{sub} (#it{x}_{hZ}) [%]");

        highs->SetLineColor (kGray+1);
        highs->SetLineStyle (1);
        highs->SetLineWidth (3);

        highs->DrawCopy ("][ hist");

        if (systematics.size () == 0)
          myLineText (0.65, 0.88, kGray+1, 1, description.c_str (), 1, 0.04);
        else
          myLineText (0.65, 0.92, kGray+1, 1, "Total", 1, 0.026);

        lows->GetXaxis ()->SetMoreLogLabels ();
        lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        lows->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
        lows->GetYaxis ()->SetTitle (useTrkPt ? "#deltaY_{sub} / Y_{sub} (#it{p}_{T}^{ch}) [%]" : "#deltaY_{sub} / Y_{sub} (#it{x}_{hZ}) [%]");

        lows->SetLineColor (kGray+1);
        lows->SetLineStyle (1);
        lows->SetLineWidth (3);

        lows->DrawCopy ("][ same hist");

        delete highs, lows;

        short iSys = 0;
        for (Systematic* sys : systematics) {

          TH1D* h = (useTrkPt ? sys->h_z_trk_zpt_sub[iSpc][iPtZ][iCent] : sys->h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
          highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
          
          SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

          highs->Scale (100.);
          lows->Scale (100.);

          highs->GetXaxis ()->SetMoreLogLabels ();
          highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
          highs->GetYaxis ()->SetTitle (useTrkPt ? "#deltaY_{sub} / Y_{sub} (#it{p}_{T}^{ch}) [%]" : "#deltaY_{sub} / Y_{sub} (#it{x}_{hZ}) [%]");

          highs->SetLineColor (colors[iSys+1]);
          highs->SetLineStyle (iSys+2);

          highs->DrawCopy ("][ hist same");

          lows->GetXaxis ()->SetMoreLogLabels ();
          lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
          lows->GetYaxis ()->SetTitle (useTrkPt ? "#deltaY_{sub} / Y_{sub} (#it{p}_{T}^{ch}) [%]" : "#deltaY_{sub} / Y_{sub} (#it{x}_{hZ}) [%]");

          lows->SetLineColor (colors[iSys+1]);
          lows->SetLineStyle (iSys+2);

          lows->DrawCopy ("][ same hist");

          myLineText (0.65, 0.89-0.026*iSys, colors[iSys+1], iSys+2, sys->description.c_str (), 1, 0.026);

          delete highs, lows;
          iSys++;
        }

        myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
        if (iCent == 0)
          myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
        else {
          myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
        }

        if (iPtZ == nPtZBins-1)
          myText (0.24, 0.86, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04);
        else
          myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

        const char* lo = GetPiString (phiLowBins[1]);
        const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
        myText (0.24, 0.80, kBlack, Form ("%s < #left|#Delta#phi#right| < %s", lo, hi), 0.04);

        c->SaveAs (Form ("%s/TrkYieldSignalSystematics/%s_%s_iPtZ%i_iCent%i.pdf", plotPath.Data (), spc, useTrkPt ? "pTTrk":"xHZ", iPtZ, iCent));

      } // end loop over centralities

    } // end loop over PtZ
  } // end loop over species
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot this set of systematics on the signal track yield
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: PlotIAASystematics (const short pSpc, const short pPtZ) {
  const char* canvasName = Form ("c_iaasys");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
    gPad->SetLogx ();
  }
  c->cd ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {

          TH1D* centralVals = nullptr, *highs = nullptr, *lows = nullptr;
          centralVals = h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];

          if (!centralVals) continue;

          highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
          SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

          highs->Scale (100.);
          lows->Scale (100.);

          highs->GetXaxis ()->SetMoreLogLabels ();
          highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
          highs->GetYaxis ()->SetTitle ("#deltaI_{AA} / I_{AA} (#it{p}_{T}^{ch}) [%]");

          highs->SetLineColor (kGray+1);
          highs->SetLineStyle (1);
          highs->SetLineWidth (3);

          highs->DrawCopy ("][ hist");

          if (systematics.size () == 0)
            myLineText (0.65, 0.88, kGray+1, 1, description.c_str (), 1, 0.04);
          else
            myLineText (0.65, 0.92, kGray+1, 1, "Total", 1, 0.026);

          lows->GetXaxis ()->SetMoreLogLabels ();
          lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
          lows->GetYaxis ()->SetTitle ("#deltaI_{AA} / I_{AA} (#it{p}_{T}^{ch}) [%]");

          lows->SetLineColor (kGray+1);
          lows->SetLineStyle (1);
          lows->SetLineWidth (3);

          lows->DrawCopy ("][ same hist");

          delete highs, lows;

          bool drawn = false;
          short iSys = 0;
          for (Systematic* sys : systematics) {

            TH1D* h = sys->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];
            highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
            lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
            
            SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

            highs->Scale (100.);
            lows->Scale (100.);

            highs->GetXaxis ()->SetMoreLogLabels ();
            highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            highs->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
            highs->GetYaxis ()->SetTitle ("I_{AA} (#it{p}_{T}^{ch}) Relative error");
            highs->GetYaxis ()->SetTitle ("#deltaI_{AA} / I_{AA} (#it{p}_{T}^{ch}) [%]");

            highs->SetLineColor (colors[iSys+1]);
            highs->SetLineStyle (iSys+2);

            if (!drawn)
              highs->DrawCopy ("][ hist");
            else
              highs->DrawCopy ("][ hist same");
            drawn = true;

            lows->GetXaxis ()->SetMoreLogLabels ();
            lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            lows->GetXaxis ()->SetTitle ("#it{p}_{T}^{ch} [GeV]");
            lows->GetYaxis ()->SetTitle ("#deltaI_{AA} / I_{AA} (#it{p}_{T}^{ch}) [%]");

            lows->SetLineColor (colors[iSys+1]);
            lows->SetLineStyle (iSys+2);

            lows->DrawCopy ("][ same hist");

            myLineText (0.65, 0.89-0.026*iSys, colors[iSys+1], iSys+2, sys->description.c_str (), 1, 0.026);

            delete highs, lows;
            iSys++;
          }

          myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
          if (iCent == 0)
            myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
          else {
            myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
          }

          if (iPtZ == nPtZBins-1)
            myText (0.24, 0.86, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04);
          else
            myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

          const char* lo = GetPiString (phiLowBins[iPhi]);
          const char* hi = GetPiString (phiHighBins[iPhi]);
          myText (0.24, 0.80, kBlack, Form ("%s < #left|#Delta#phi#right| < %s", lo, hi), 0.04);

          c->SaveAs (Form ("%s/IAASystematics/%s_iPtZ%i_iPhi%i_iCent%i.pdf", plotPath.Data (), spc, iPtZ, iPhi, iCent));

        } // end loop over centralities

      } // end loop over Phi bins

    } // end loop over PtZ
  } // end loop over species
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot this set of systematics on the signal track yield
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: PlotIAASystematicsPtZ (const bool useTrkPt, const short pSpc) {
  const char* canvasName = Form ("c_iaasys");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
    gDirectory->Add (c);
    gPad->SetLogx ();
  }
  c->cd ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 1; iCent < numCentBins; iCent++) {

        TH1D* centralVals = nullptr, *highs = nullptr, *lows = nullptr;
        centralVals = (useTrkPt ? h_z_trk_zpt_iaa[iSpc][iPtZ][iCent] : h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);

        if (!centralVals) continue;

        highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
        lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
        SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

        highs->Scale (100.);
        lows->Scale (100.);

        highs->GetXaxis ()->SetMoreLogLabels ();
        highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        highs->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
        highs->GetYaxis ()->SetTitle (useTrkPt ? "#deltaI_{AA} / I_{AA} (#it{p}_{T}^{ch}) [%]" : "#deltaI_{AA} / I_{AA} (#it{x}_{hZ}) [%]");

        highs->SetLineColor (kGray+1);
        highs->SetLineStyle (1);
        highs->SetLineWidth (3);

        highs->DrawCopy ("][ hist");

        if (systematics.size () == 0)
          myLineText (0.65, 0.88, kGray+1, 1, description.c_str (), 1, 0.04);
        else
          myLineText (0.65, 0.92, kGray+1, 1, "Total", 1, 0.026);

        lows->GetXaxis ()->SetMoreLogLabels ();
        lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        lows->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
        lows->GetYaxis ()->SetTitle (useTrkPt ? "#deltaI_{AA} / I_{AA} (#it{p}_{T}^{ch}) [%]" : "#deltaI_{AA} / I_{AA} (#it{x}_{hZ}) [%]");

        lows->SetLineColor (kGray+1);
        lows->SetLineStyle (1);
        lows->SetLineWidth (3);

        lows->DrawCopy ("][ same hist");

        delete highs, lows;

        short iSys = 0;
        for (Systematic* sys : systematics) {

          TH1D* h = (useTrkPt ? sys->h_z_trk_zpt_iaa[iSpc][iPtZ][iCent] : sys->h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
          highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
          
          SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

          highs->Scale (100.);
          lows->Scale (100.);

          highs->GetXaxis ()->SetMoreLogLabels ();
          highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
          highs->GetYaxis ()->SetTitle (useTrkPt ? "#deltaI_{AA} / I_{AA} (#it{p}_{T}^{ch}) [%]" : "#deltaI_{AA} / I_{AA} (#it{x}_{hZ}) [%]");

          highs->SetLineColor (colors[iSys+1]);
          highs->SetLineStyle (iSys+2);

          highs->DrawCopy ("][ hist same");

          lows->GetXaxis ()->SetMoreLogLabels ();
          lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ch} [GeV]" : "#it{x}_{hZ}");
          lows->GetYaxis ()->SetTitle (useTrkPt ? "#deltaI_{AA} / I_{AA} (#it{p}_{T}^{ch}) [%]" : "#deltaI_{AA} / I_{AA} (#it{x}_{hZ}) [%]");

          lows->SetLineColor (colors[iSys+1]);
          lows->SetLineStyle (iSys+2);

          lows->DrawCopy ("][ same hist");

          myLineText (0.65, 0.89-0.026*iSys, colors[iSys+1], iSys+2, sys->description.c_str (), 1, 0.026);

          delete highs, lows;
          iSys++;
        }

        myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
        if (iCent == 0)
          myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
        else {
          myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
        }

        if (iPtZ == nPtZBins-1)
          myText (0.24, 0.86, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.04);
        else
          myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

        const char* lo = GetPiString (phiLowBins[1]);
        const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
        myText (0.24, 0.80, kBlack, Form ("%s < #left|#Delta#phi#right| < %s", lo, hi), 0.04);

        c->SaveAs (Form ("%s/IAASystematics/%s_%s_iPtZ%i_iCent%i.pdf", plotPath.Data (), spc, useTrkPt ? "pTTrk":"xHZ", iPtZ, iCent));

      } // end loop over centralities

    } // end loop over PtZ
  } // end loop over species
  
}  




void Systematic :: PrintIAA (const bool printErrs, const bool useTrkPt, const short iCent, const short iPtZ, const short iSpc) {
  TGAE* g = nullptr;
  if (useTrkPt)
    g = GetTGAE (h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]);
  else
    g = GetTGAE (h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
  cout << name << "\t";
  for (int ix = 0; ix < g->GetN (); ix++) {
    double x, y;
    g->GetPoint (ix, x, y);
    cout << y << "\t";
    if (printErrs)
      cout << "+" << g->GetErrorYhigh (ix) << ",\t-" << g->GetErrorYlow (ix) << "\t";
    cout << endl;
  }
}




void Systematic :: WriteIAAs () {
  SetupDirectories ("", "ZTrackAnalysis/");

  TDirectory* _gDirectory = gDirectory;

  const char* outFileName = "DataAnalysis/Nominal/data18hi_iaa_fits.root"; 
  TFile* outFile = new TFile (Form ("%s/%s", rootPath.Data (), outFileName), "recreate");

  TF1**** f_z_trk_zpt_iaa = Get3DArray <TF1*> (3, nPtZBins, numCentBins);
  TF1**** f_z_trk_zxzh_iaa = Get3DArray <TF1*> (3, nPtZBins, numCentBins);

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent] = new TF1 (Form ("f_z_trk_zpt_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "[0]+[1]*log(x)+[2]*(log(x))^2+[3]*(log(x))^3", ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]);
        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->SetParameter (0, 1);
        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->SetParameter (1, 0);
        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->SetParameter (2, 0);
        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->SetParameter (3, 0);
        h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->Fit (f_z_trk_zpt_iaa[iSpc][iPtZ][iCent], "RN0Q");

        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent] = new TF1 (Form ("f_z_trk_zxzh_iaa_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()), "[0]+[1]*log(x)+[2]*(log(x))^2+[3]*(log(x))^3", xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[iPtZ]]);
        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->SetParameter (0, 1);
        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->SetParameter (1, 0);
        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->SetParameter (2, 0);
        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->SetParameter (3, 0);
        h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->Fit (f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent], "RN0Q");

        f_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->Write ();
        h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]->Write ();
        f_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->Write ();
        h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]->Write ();
      }
    }
  }

  outFile->Close ();

  Delete3DArray (f_z_trk_zpt_iaa, 3, nPtZBins, numCentBins);
  Delete3DArray (f_z_trk_zxzh_iaa, 3, nPtZBins, numCentBins);

  _gDirectory->cd ();
}

   
#endif
