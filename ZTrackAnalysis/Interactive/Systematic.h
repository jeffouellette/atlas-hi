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
  map <PhysicsAnalysis*, short> variationDirs; // variation direction values are -1 for down, 0 for both, 1 for up
  vector<Systematic*> systematics;

  // map taking TH1Ds to TGAEs for systematics
  map <TH1D*, TGAE*> graphMap;

  // TFile with systematics graphs
  TFile* graphFile = nullptr;
  bool graphsLoaded = false;
  

  void NullifyErrors ();
  void AddGraphPair (TH1D* h);
  void CreateSysGraphs ();

  public:

  string description;

  Systematic (PhysicsAnalysis* nom, const char* _name = "systematics", const char* _desc = "systematic") : PhysicsAnalysis (){
    name = _name;
    directory = "Systematics/";
    description = _desc;
    SetupDirectories (directory, "ZTrackAnalysis/");
    CopyAnalysis (nom, true);
    CreateSysGraphs ();
    NullifyErrors ();
  }

  vector<PhysicsAnalysis*>& GetVariations ()  { return variations;  }
  vector<Systematic*>&      GetSystematics () { return systematics; }

  virtual TGAE* GetTGAE (TH1D* h) override;

  virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;

  void AddVariation (PhysicsAnalysis* a, const short dir = 0);
  void AddSystematic (Systematic* a);

  void AddVariations (); // variations add linearly
  void AddSystematics (); // systematics add in quadrature

  void PlotTrkYieldSystematics (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  void PlotTrkYieldSystematicsPtZ (const bool useTrkPt = true, const short pSpc = 2);
  void PlotSignalTrkYieldSystematics (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  void PlotSignalTrkYieldSystematicsPtZ (const bool useTrkPt = true, const short pSpc = 2);
  void PlotIAASystematics (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  void PlotIAASystematicsPtZ (const bool useTrkPt = true, const short pSpc = 2);

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
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));

      for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {

          AddGraphPair (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);

          AddGraphPair (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
          AddGraphPair (h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]);
          
          AddGraphPair (h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
          AddGraphPair (h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent]);
          AddGraphPair (h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent]);

          AddGraphPair (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
          AddGraphPair (h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]);
          AddGraphPair (h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
          AddGraphPair (h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]);
          AddGraphPair (h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over iPhi

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          AddGraphPair (h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]);
          AddGraphPair (h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent]);
        } // end loop over iPtTrk

        AddGraphPair (h_z_trk_zpt[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zpt_sub[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zpt_icp[iSpc][iPtZ][iCent]);

        AddGraphPair (h_z_trk_zxzh[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
        AddGraphPair (h_z_trk_zxzh_icp[iSpc][iPtZ][iCent]);
        
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




void Systematic :: LoadHists (const char* histFileName, const bool _finishHists) {
  PhysicsAnalysis :: LoadHists (histFileName, _finishHists);
  CreateSysGraphs ();
}
//void Systematic :: LoadHists (const char* histFileName) {
//  PhysicsAnalysis :: LoadHists (histFileName, false);
//
//  for (short iCent = 0; iCent < numCentBins; iCent++) {
//    for (short iSpc = 0; iSpc < 3; iSpc++) {
//      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
//          h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_phi_sub_%s_iPtZ%i_iPtTrk%i_iCent%i_%s", spc, iPtZ, iPtTrk, iCent, name.c_str ()));
//        }
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_pt_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//          h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_pt_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//          h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_xzh_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//          h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_xzh_sig_to_bkg_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
//        }
//        h_z_trk_zpt_sub[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_zpt_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
//        h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_zpt_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
//        h_z_trk_zxzh_sub[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_zxzh_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
//        h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent] = (TH1D*) histFile->Get (Form ("h_z_trk_zxzh_sig_to_bkg_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
//      }
//    }
//  }
//}




//void Systematic :: SaveHists (const char* histFileName) {
//  SetupDirectories (directory.c_str (), "ZTrackAnalysis/");
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
//      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
//          SafeWrite (h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]);
//        }
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          SafeWrite (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
//          SafeWrite (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
//          SafeWrite (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
//        }
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




//void Systematic :: LoadGraphs (const char* graphFileName) {
//  if (!histsLoaded)
//    loadHists (
//
//  SetupDirectories (directory.c_str (), "ZTrackAnalysis/");
//
//  TDirectory* _gDirectory = gDirectory;
//  graphFile = new TFile (Form ("%s/%s", rootPath.Data (), graphFileName), "read");
//
//  _gDirectory->cd ();
//}




////////////////////////////////////////////////////////////////////////////////////////////////
// Adds a variation to consider
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: AddVariation (PhysicsAnalysis* a, const short dir) {
  if (find (variations.begin (), variations.end (), a) == variations.end ()) {
    variations.push_back (a);
    variationDirs.insert (pair <PhysicsAnalysis*, short> (a, dir));
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

  NullifyErrors ();

  for (PhysicsAnalysis* a : variations) {

    cout << "Adding variation " << a->Name () << " to systematic " << name << endl;

    a->SubtractBackground ();
    a->CalculateIAA ();
    a->CalculateICP ();

    TGAE* sys = nullptr;
    TH1D* var = nullptr;

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) { 

        // Hadron yield systematics, signal & signal+bkg levels
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            sys = GetTGAE (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);

            sys = GetTGAE (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
            sys = GetTGAE (h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
            sys = GetTGAE (h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);

            sys = GetTGAE (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
            sys = GetTGAE (h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
            sys = GetTGAE (h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          } // end loop over cents
        } // end loop over phi

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
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
          sys = GetTGAE (h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);

          sys = GetTGAE (h_z_trk_zxzh[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zxzh[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          sys = GetTGAE (h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zxzh_sub[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          sys = GetTGAE (h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        } // end loop over cents


        // IAA, ICP systematics
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 1; iCent < numCentBins; iCent++) {
            sys = GetTGAE (h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
            sys = GetTGAE (h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          } // end loop over cents

          for (short iCent = 2; iCent < numCentBins; iCent++) {
            sys = GetTGAE (h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
            sys = GetTGAE (h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
          } // end loop over cents
        } // end loop over phi

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          sys = GetTGAE (h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zpt_iaa[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);

          sys = GetTGAE (h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        }

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          sys = GetTGAE (h_z_trk_zpt_icp[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zpt_icp[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);

          sys = GetTGAE (h_z_trk_zxzh_icp[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zxzh_icp[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, variationDirs[a]);
        } // end loop over cents

      } // end loop over pT^Z bins
    } // end loop over species

  } // end loop over variations
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

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) { 

        // Hadron yield systematics, signal & signal+bkg levels
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            master = GetTGAE (h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);

            master = GetTGAE (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = GetTGAE (h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = GetTGAE (h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);

            master = GetTGAE (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = GetTGAE (h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = GetTGAE (h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);
          } // end loop over cents
        } // end loop over phi

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
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
          master = GetTGAE (h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);

          master = GetTGAE (h_z_trk_zxzh[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zxzh[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);
          master = GetTGAE (h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);
          master = GetTGAE (h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);
        } // end loop over cents


        // IAA, ICP systematics
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 1; iCent < numCentBins; iCent++) {
            master = GetTGAE (h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = GetTGAE (h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);
          } // end loop over cents
          for (short iCent = 2; iCent < numCentBins; iCent++) {
            master = GetTGAE (h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = GetTGAE (h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]);
            sys = s->GetTGAE (s->h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent]);
            if (master && sys) AddErrorsInQuadrature (master, sys);
          } // end loop over cents
        } // end loop over phi

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          master = GetTGAE (h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);

          master = GetTGAE (h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);
        }

        for (short iCent = 1; iCent < numCentBins; iCent++) {
          master = GetTGAE (h_z_trk_zpt_icp[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zpt_icp[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);

          master = GetTGAE (h_z_trk_zxzh_icp[iSpc][iPtZ][iCent]);
          sys = s->GetTGAE (s->h_z_trk_zxzh_icp[iSpc][iPtZ][iCent]);
          if (master && sys) AddErrorsInQuadrature (master, sys);
        } // end loop over cents
      } // end loop over pT^Z bins
    } // end loop over species

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
    for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {

          TH1D* centralVals = nullptr, *highs = nullptr, *lows = nullptr;
          centralVals = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];

          if (!centralVals) continue;

          bool drawn = false;
          short iSys = 0;
          for (Systematic* sys : systematics) {

            if (sys->Name () == string ("bkgSys"))
              continue;

            TH1D* h = sys->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
            highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
            lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
            
            SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

            highs->GetXaxis ()->SetMoreLogLabels ();
            if (iCent == 0)
              highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
            else
              highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            highs->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
            highs->GetYaxis ()->SetTitle ("Y (#it{p}_{T}) Relative error");

            highs->SetLineColor (colors[iSys+1]);
            highs->SetLineStyle (2);
            highs->SetLineWidth (5);

            if (!drawn)
              highs->DrawCopy ("][ hist");
            else
              highs->DrawCopy ("][ hist same");
            drawn = true;

            lows->GetXaxis ()->SetMoreLogLabels ();
            if (iCent == 0)
              lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
            else
              lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            lows->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
            lows->GetYaxis ()->SetTitle ("Y (#it{p}_{T}) Relative error");

            lows->SetLineColor (colors[iSys+1]);
            lows->SetLineStyle (2);
            lows->SetLineWidth (5);

            lows->DrawCopy ("][ same hist");

            myText (0.65, 0.89-0.026*iSys, colors[iSys+1], sys->description.c_str (), 0.026);

            //delete errs;
            delete highs, lows;
            iSys++;
          }

          highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
          SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

          highs->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->SetLineColor (kBlack);
          highs->SetLineStyle (1);
          highs->SetLineWidth (3);

          myText (0.65, 0.92, kBlack, "Total", 0.026);
          
          highs->DrawCopy (systematics.size () == 0 ? "][ hist" : "][ same hist");

          lows->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->SetLineColor (kBlack);
          lows->SetLineStyle (1);
          lows->SetLineWidth (3);

          lows->DrawCopy ("][ same hist");

          delete highs, lows;

          myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
          if (iCent == 0)
            myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
          else {
            myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
          }

          if (iPtZ == nPtZBins-1)
            myText (0.24, 0.86, kBlack, Form ("#it{p}_{T} > %g GeV", zPtBins[iPtZ]), 0.04);
          else
            myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

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
    for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {

        TH1D* centralVals = nullptr, *highs = nullptr, *lows = nullptr;
        centralVals = (useTrkPt ? h_z_trk_zpt[iSpc][iPtZ][iCent] : h_z_trk_zxzh[iSpc][iPtZ][iCent]);

        if (!centralVals) continue;

        bool drawn = false;
        short iSys = 0;
        for (Systematic* sys : systematics) {

          if (sys->Name () == string ("bkgSys"))
            continue;

          TH1D* h = (useTrkPt ? sys->h_z_trk_zpt[iSpc][iPtZ][iCent] : sys->h_z_trk_zxzh[iSpc][iPtZ][iCent]);
          highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
          
          SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

          highs->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T} [GeV]" : "#it{x}_{hZ}");
          highs->GetYaxis ()->SetTitle (useTrkPt ? "Y (#it{p}_{T}) Relative error" : "Y (#it{x}_{hZ}) Relative error");

          highs->SetLineColor (colors[iSys+1]);
          highs->SetLineStyle (2);
          highs->SetLineWidth (5);

          if (!drawn)
            highs->DrawCopy ("][ hist");
          else
            highs->DrawCopy ("][ hist same");
          drawn = true;

          lows->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T} [GeV]" : "#it{x}_{hZ}");
          lows->GetYaxis ()->SetTitle (useTrkPt ? "Y (#it{p}_{T}) Relative error" : "Y (#it{x}_{hZ}) Relative error");

          lows->SetLineColor (colors[iSys+1]);
          lows->SetLineStyle (2);
          lows->SetLineWidth (5);

          lows->DrawCopy ("][ same hist");

          myText (0.65, 0.89-0.026*iSys, colors[iSys+1], sys->description.c_str (), 0.026);

          //delete errs;
          delete highs, lows;
          iSys++;
        }

        highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
        lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
        SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

        highs->GetXaxis ()->SetMoreLogLabels ();
        if (iCent == 0)
          highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
        else
          highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        highs->SetLineColor (kBlack);
        highs->SetLineStyle (1);
        highs->SetLineWidth (3);

        myText (0.65, 0.92, kBlack, "Total", 0.026);
        
        highs->DrawCopy (systematics.size () == 0 ? "][ hist" : "][ same hist");

        lows->GetXaxis ()->SetMoreLogLabels ();
        if (iCent == 0)
          lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
        else
          lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        lows->SetLineColor (kBlack);
        lows->SetLineStyle (1);
        lows->SetLineWidth (3);

        lows->DrawCopy ("][ same hist");

        delete highs, lows;

        myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
        if (iCent == 0)
          myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
        else {
          myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
        }

        if (iPtZ == nPtZBins-1)
          myText (0.24, 0.86, kBlack, Form ("#it{p}_{T} > %g GeV", zPtBins[iPtZ]), 0.04);
        else
          myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

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
    for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {

          gPad->Clear ();

          TH1D* centralVals = nullptr, *highs = nullptr, *lows = nullptr;
          centralVals = h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];

          if (!centralVals) continue;

          bool drawn = false;
          short iSys = 0;
          for (Systematic* sys : systematics) {

            TH1D* h = sys->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
            highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
            lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
            
            SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

            highs->GetXaxis ()->SetMoreLogLabels ();
            if (iCent == 0)
              highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
            else
              highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            highs->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
            highs->GetYaxis ()->SetTitle ("Y_{sub} (#it{p}_{T}) Relative error");

            highs->SetLineColor (colors[iSys+1]);
            highs->SetLineStyle (2);
            highs->SetLineWidth (5);

            if (!drawn)
              highs->DrawCopy ("][ hist");
            else
              highs->DrawCopy ("][ hist same");
            drawn = true;

            lows->GetXaxis ()->SetMoreLogLabels ();
            if (iCent == 0)
              lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
            else
              lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            lows->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
            lows->GetYaxis ()->SetTitle ("Y_{sub} (#it{p}_{T}) Relative error");

            lows->SetLineColor (colors[iSys+1]);
            lows->SetLineStyle (2);
            lows->SetLineWidth (5);

            lows->DrawCopy ("][ same hist");

            myText (0.65, 0.89-0.026*iSys, colors[iSys+1], sys->description.c_str (), 0.026);

            delete highs, lows;
            iSys++;
          }


          highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
          SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

          highs->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->SetLineColor (kBlack);
          highs->SetLineStyle (1);
          highs->SetLineWidth (3);

          myText (0.65, 0.92, kBlack, "Total", 0.026);
          
          highs->DrawCopy (systematics.size () == 0 ? "][ hist" : "][ same hist");

          lows->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->SetLineColor (kBlack);
          lows->SetLineStyle (1);
          lows->SetLineWidth (3);

          lows->DrawCopy ("][ same hist");

          delete highs, lows;

          myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
          if (iCent == 0)
            myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
          else {
            myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
          }

          if (iPtZ == nPtZBins-1)
            myText (0.24, 0.86, kBlack, Form ("#it{p}_{T} > %g GeV", zPtBins[iPtZ]), 0.04);
          else
            myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

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

        bool drawn = false;
        short iSys = 0;
        for (Systematic* sys : systematics) {

          TH1D* h = (useTrkPt ? sys->h_z_trk_zpt_sub[iSpc][iPtZ][iCent] : sys->h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
          highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
          
          SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

          highs->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T} [GeV]" : "#it{x}_{hZ}");
          highs->GetYaxis ()->SetTitle (useTrkPt ? "Y_{sub} (#it{p}_{T}) Relative error" : "Y_{sub} (#it{x}_{hZ}) Relative error");

          highs->SetLineColor (colors[iSys+1]);
          highs->SetLineStyle (2);
          highs->SetLineWidth (5);

          if (!drawn)
            highs->DrawCopy ("][ hist");
          else
            highs->DrawCopy ("][ hist same");
          drawn = true;

          lows->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T} [GeV]" : "#it{x}_{hZ}");
          lows->GetYaxis ()->SetTitle (useTrkPt ? "Y_{sub} (#it{p}_{T}) Relative error" : "Y_{sub} (#it{x}_{hZ}) Relative error");

          lows->SetLineColor (colors[iSys+1]);
          lows->SetLineStyle (2);
          lows->SetLineWidth (5);

          lows->DrawCopy ("][ same hist");

          myText (0.65, 0.89-0.026*iSys, colors[iSys+1], sys->description.c_str (), 0.026);

          delete highs, lows;
          iSys++;
        }


        highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
        lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
        SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

        highs->GetXaxis ()->SetMoreLogLabels ();
        if (iCent == 0)
          highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
        else
          highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        highs->SetLineColor (kBlack);
        highs->SetLineStyle (1);
        highs->SetLineWidth (3);

        myText (0.65, 0.92, kBlack, "Total", 0.026);
        
        highs->DrawCopy (systematics.size () == 0 ? "][ hist" : "][ same hist");

        lows->GetXaxis ()->SetMoreLogLabels ();
        if (iCent == 0)
          lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
        else
          lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        lows->SetLineColor (kBlack);
        lows->SetLineStyle (1);
        lows->SetLineWidth (3);

        lows->DrawCopy ("][ same hist");

        delete highs, lows;

        myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
        if (iCent == 0)
          myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
        else {
          myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
        }

        if (iPtZ == nPtZBins-1)
          myText (0.24, 0.86, kBlack, Form ("#it{p}_{T} > %g GeV", zPtBins[iPtZ]), 0.04);
        else
          myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

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
    for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {

          TH1D* centralVals = nullptr, *highs = nullptr, *lows = nullptr;
          centralVals = h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];

          if (!centralVals) continue;

          bool drawn = false;
          short iSys = 0;
          for (Systematic* sys : systematics) {

            TH1D* h = sys->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];
            highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
            lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
            
            SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

            highs->GetXaxis ()->SetMoreLogLabels ();
            if (iCent == 0)
              highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
            else
              highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            highs->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
            highs->GetYaxis ()->SetTitle ("I_{AA} (#it{p}_{T}) Relative error");

            highs->SetLineColor (colors[iSys+1]);
            highs->SetLineStyle (2);
            highs->SetLineWidth (5);

            if (!drawn)
              highs->DrawCopy ("][ hist");
            else
              highs->DrawCopy ("][ hist same");
            drawn = true;

            lows->GetXaxis ()->SetMoreLogLabels ();
            if (iCent == 0)
              lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
            else
              lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

            lows->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
            lows->GetYaxis ()->SetTitle ("I_{AA} (#it{p}_{T}) Relative error");

            lows->SetLineColor (colors[iSys+1]);
            lows->SetLineStyle (2);
            lows->SetLineWidth (5);

            lows->DrawCopy ("][ same hist");

            myText (0.65, 0.89-0.026*iSys, colors[iSys+1], sys->description.c_str (), 0.026);

            delete highs, lows;
            iSys++;
          }

          highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
          SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

          highs->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->SetLineColor (kBlack);
          highs->SetLineStyle (1);
          highs->SetLineWidth (3);

          myText (0.65, 0.92, kBlack, "Total", 0.026);
          
          highs->DrawCopy (systematics.size () == 0 ? "][ hist" : "][ same hist");

          lows->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->SetLineColor (kBlack);
          lows->SetLineStyle (1);
          lows->SetLineWidth (3);

          lows->DrawCopy ("][ same hist");

          delete highs, lows;

          myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
          if (iCent == 0)
            myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
          else {
            myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
          }

          if (iPtZ == nPtZBins-1)
            myText (0.24, 0.86, kBlack, Form ("#it{p}_{T} > %g GeV", zPtBins[iPtZ]), 0.04);
          else
            myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

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
        centralVals = h_z_trk_zpt_iaa[iSpc][iPtZ][iCent];

        if (!centralVals) continue;

        bool drawn = false;
        short iSys = 0;
        for (Systematic* sys : systematics) {

          TH1D* h = sys->h_z_trk_zpt_iaa[iSpc][iPtZ][iCent];
          highs = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysHigh").c_str ());
          lows = (TH1D*) h->Clone ((string (h->GetName ()) + "_relSysLow").c_str ());
          
          SaveRelativeErrors (sys->GetTGAE (h), GetTGAE (centralVals), highs, lows);

          highs->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          highs->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T} [GeV]" : "#it{x}_{hZ}");
          highs->GetYaxis ()->SetTitle (useTrkPt ? "I_{AA} (#it{p}_{T}) Relative error" : "I_{AA} (#it{x}_{hZ}) Relative error");

          highs->SetLineColor (colors[iSys+1]);
          highs->SetLineStyle (2);
          highs->SetLineWidth (5);

          if (!drawn)
            highs->DrawCopy ("][ hist");
          else
            highs->DrawCopy ("][ hist same");
          drawn = true;

          lows->GetXaxis ()->SetMoreLogLabels ();
          if (iCent == 0)
            lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
          else
            lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

          lows->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T} [GeV]" : "#it{x}_{hZ}");
          lows->GetYaxis ()->SetTitle (useTrkPt ? "I_{AA} (#it{p}_{T}) Relative error" : "I_{AA} (#it{x}_{hZ}) Relative error");

          lows->SetLineColor (colors[iSys+1]);
          lows->SetLineStyle (2);
          lows->SetLineWidth (5);

          lows->DrawCopy ("][ same hist");

          myText (0.65, 0.89-0.026*iSys, colors[iSys+1], sys->description.c_str (), 0.026);

          delete highs, lows;
          iSys++;
        }

        highs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysHigh").c_str ());
        lows = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSysLow").c_str ());
        SaveRelativeErrors (GetTGAE (centralVals), GetTGAE (centralVals), highs, lows);

        highs->GetXaxis ()->SetMoreLogLabels ();
        if (iCent == 0)
          highs->GetYaxis ()->SetRangeUser (-0.1, 0.1);
        else
          highs->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        highs->SetLineColor (kBlack);
        highs->SetLineStyle (1);
        highs->SetLineWidth (3);

        myText (0.65, 0.92, kBlack, "Total", 0.026);
        
        highs->DrawCopy (systematics.size () == 0 ? "][ hist" : "][ same hist");

        lows->GetXaxis ()->SetMoreLogLabels ();
        if (iCent == 0)
          lows->GetYaxis ()->SetRangeUser (-0.1, 0.1);
        else
          lows->GetYaxis ()->SetRangeUser (-max_rel_sys, max_rel_sys);

        lows->SetLineColor (kBlack);
        lows->SetLineStyle (1);
        lows->SetLineWidth (3);

        lows->DrawCopy ("][ same hist");

        delete highs, lows;

        myText (0.24, 0.28, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
        if (iCent == 0)
          myText (0.24, 0.22, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
        else {
          myText (0.24, 0.22, kBlack, Form ("Pb+Pb 5.02 TeV, %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
        }

        if (iPtZ == nPtZBins-1)
          myText (0.24, 0.86, kBlack, Form ("#it{p}_{T} > %g GeV", zPtBins[iPtZ]), 0.04);
        else
          myText (0.24, 0.86, kBlack, Form ("%g < #it{p}_{T} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);

        const char* lo = GetPiString (phiLowBins[1]);
        const char* hi = GetPiString (phiHighBins[numPhiBins-1]);
        myText (0.24, 0.80, kBlack, Form ("%s < #left|#Delta#phi#right| < %s", lo, hi), 0.04);

        c->SaveAs (Form ("%s/IAASystematics/%s_%s_iPtZ%i_iCent%i.pdf", plotPath.Data (), spc, useTrkPt ? "pTTrk":"xHZ", iPtZ, iCent));

      } // end loop over centralities

    } // end loop over PtZ
  } // end loop over species
  
}

#endif
