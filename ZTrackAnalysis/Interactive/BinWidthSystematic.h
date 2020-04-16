#ifndef __BinWidthSystematic_h__
#define __BinWidthSystematic_h__

#include "Params.h"
#include "PhysicsAnalysis.h"
#include "Systematic.h"

#include <GlobalParams.h>
#include <Utilities.h>

#include <iostream>

using namespace std;
using namespace atlashi;


class BinWidthSystematic : public Systematic {

  public:

  BinWidthSystematic (FullAnalysis* _nom, const char* _name, const char* _desc) : Systematic () {
    name = _name;
    description = _desc;
    nom = _nom;
    plotAsSyst = true;
    PhysicsAnalysis :: CopyAnalysis (nom, true);
    CalculateTrackMeans (_nom, _nom->h_z_pt);
    CreateSysGraphs (); 
    NullifyErrors ();

    isBinWidthSys = true;
  }

  BinWidthSystematic (FullAnalysis* _nom, const char* _name, const char* _desc, const char* _graphFileName) : Systematic () {
    name = _name;
    description = _desc;
    nom = _nom; 
    plotAsSyst = true;
    PhysicsAnalysis :: CopyAnalysis (nom, true);
    CalculateTrackMeans (_nom, _nom->h_z_pt);
    LoadGraphs (_graphFileName);

    isBinWidthSys = true;
  }

  //virtual void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;
  virtual void LoadGraphs (const char* graphFileName) override;

  virtual void AddVariations () override; // variations add linearly
  virtual void AddVariationsUsingStdDev () override; // error is set to the standard deviation of several variations 
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void BinWidthSystematic :: LoadGraphs (const char* graphFileName) {
  SetupDirectories ("", "ZTrackAnalysis/");

  TDirectory* _gDirectory = gDirectory;
  graphFile = new TFile (Form ("%s/%s", rootPath.Data (), graphFileName), "read");

  TH1D* h = nullptr;
  TGAE* g = nullptr;

  if (isBinWidthSys)
    return; // don't need to create extra graphs where systematic doesn't apply

  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
        for (short iPtch = 0; iPtch < nPtchBins[iPtZ]; iPtch++) {
          h = h_trk_dphi_sub[iSpc][iPtZ][iPtch][iCent];
          g = (TGAE*) graphFile->Get (Form ("g_trk_dphi_sub_%s_iPtZ%i_iPtch%i_iCent%i_%s", spc, iPtZ, iPtch, iCent, name.c_str ()));
          graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        }
        //for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        //  h = h_trk_pt_dphi_sub[iSpc][iPtZ][iPhi][iCent];
        //  g = (TGAE*) graphFile->Get (Form ("g_trk_pt_dphi_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        //  graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        //  h = h_trk_xhz_dphi_sub[iSpc][iPtZ][iPhi][iCent];
        //  g = (TGAE*) graphFile->Get (Form ("g_z_trk_xzh_sub_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
        //  graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        //}
        h = h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];
        g = (TGAE*) graphFile->Get (Form ("g_trk_pt_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
        h = h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
        g = (TGAE*) graphFile->Get (Form ("g_trk_xhz_ptz_sub_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        graphMap.insert (std::pair <TH1D*, TGAE*> (h, g));
      }
    }
  }

  _gDirectory->cd ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Sets the errors in this systematic as a combination of all added variations.
// Takes the maximum error for each point.
// Intended for combining up & down variations, but expandable for additional categories
// (e.g. track quality criteria)
////////////////////////////////////////////////////////////////////////////////////////////////
void BinWidthSystematic :: AddVariations () {

  //NullifyErrors ();
  cout << "Evaluating " << description << " systematics by taking maximum variations" << endl;

  for (PhysicsAnalysis* a : variations) {

    cout << "Adding variation " << a->Name () << " to systematic " << name << endl;

    TGAE* sys = nullptr;
    TH1D* var = nullptr;

    for (short iSpc = 0; iSpc < 3; iSpc++) {

      if (trackMeansCalculated) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          sys = g_trk_avg_pt_ptz[iSpc][iCent];
          TGAE* gvar = a->g_trk_avg_pt_ptz[iSpc][iCent];
          if (sys && gvar) {
            if (meanTrackUncStoredAtCentralValues) CalcSystematics (sys, sys, gvar, gvar, true);
            else {
              double xle, xhe, yle, yhe;
              for (int ix = 0; ix < sys->GetN (); ix++) {
                xle = fmax (sys->GetErrorXlow (ix), gvar->GetErrorXlow (ix));
                xhe = fmax (sys->GetErrorXhigh (ix), gvar->GetErrorXhigh (ix));
                yle = fmax (sys->GetErrorYlow (ix), gvar->GetErrorYlow (ix));
                yhe = fmax (sys->GetErrorYhigh (ix), gvar->GetErrorYhigh (ix));
                sys->SetPointEXlow (ix, xle);
                sys->SetPointEXhigh (ix, xhe);
                sys->SetPointEYlow (ix, yle);
                sys->SetPointEYhigh (ix, yhe);
              } // end loop over ix
            }
          }

          sys = g_trk_avg_xhz_ptz[iSpc][iCent];
          gvar = a->g_trk_avg_xhz_ptz[iSpc][iCent];
          if (sys && gvar) {
            if (meanTrackUncStoredAtCentralValues) CalcSystematics (sys, sys, gvar, gvar, true);
            else {
              double xle, xhe, yle, yhe;
              for (int ix = 0; ix < sys->GetN (); ix++) {
                xle = fmax (sys->GetErrorXlow (ix), gvar->GetErrorXlow (ix));
                xhe = fmax (sys->GetErrorXhigh (ix), gvar->GetErrorXhigh (ix));
                yle = fmax (sys->GetErrorYlow (ix), gvar->GetErrorYlow (ix));
                yhe = fmax (sys->GetErrorYhigh (ix), gvar->GetErrorYhigh (ix));
                sys->SetPointEXlow (ix, xle);
                sys->SetPointEXhigh (ix, xhe);
                sys->SetPointEYlow (ix, yle);
                sys->SetPointEYhigh (ix, yhe);
              } // end loop over ix
            }
          }
        } // end loop over iCent
      }
    } // end loop over iSpc

  } // end loop over variations

  //variations.clear ();
  //variationDirs.clear ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Sets the errors in this systematic as a combination of all added variations.
// Instead of using the maximum variation, the error is set as the standard deviation of the
// variations.
// Intended for systematics derived from several qualitatively different analysis procedures.
// (e.g. editing the event mixing criteria)
////////////////////////////////////////////////////////////////////////////////////////////////
void BinWidthSystematic :: AddVariationsUsingStdDev () {
  //TODO implement me!
  return;
}


#endif
