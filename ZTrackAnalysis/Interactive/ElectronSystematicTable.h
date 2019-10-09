#ifndef __ElectronSystematicTable_h__
#define __ElectronSystematicTable_h__

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

class ElectronSystematicTable {

  private:
  string name;

  TH1D*** relVarPt = Get2DArray <TH1D*> (nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D*** relVarX  = Get2DArray <TH1D*> (nPtZBins, numCentBins); // iSpc, iPtZ, iCent

  public:

  ElectronSystematicTable (const char* _name = "systematics") {
    name = _name;
  }

  ~ElectronSystematicTable () {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        SaferDelete (relVarPt[iPtZ][iCent]);
        SaferDelete (relVarX[iPtZ][iCent]);
      }
    }
    Delete2DArray (relVarPt, nPtZBins, numCentBins);
    Delete2DArray (relVarX, nPtZBins, numCentBins);
  }

  string Name ()              { return name; }
  void SetName (string _name) { name = _name; }

  void GetRelativeVariations (const string inFileName);
  void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);

  //virtual void AddVariations (); // variations add linearly
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void ElectronSystematicTable :: GetRelativeVariations (const string inFileName) {
  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/DataAnalysis/Variations/ElectronPtVariation/systematics_electron_scale.root", rootPath.Data ()), "read");

  cout << "Relative variation histograms found in " << Form ("%s/DataAnalysis/Variations/ElectronPtVariation/systematics_electron_scale.root", rootPath.Data ()) << endl;

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      const char* cent = (iCent == 3 ? "CENT3" : Form ("CENT%i", numCentBins-iCent-1));

      relVarPt[iPtZ][iCent] = (TH1D*) inFile->Get (Form ("h_hpt_relsys_egs_ZPT%i_%s", iPtZ-2, cent));
      relVarX[iPtZ][iCent] = (TH1D*) inFile->Get (Form ("h_xhz_relsys_egs_ZPT%i_%s", iPtZ-2, cent));
    } // end loop over iCent
  } // end loop over iPtZ
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void ElectronSystematicTable :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  for (short iSpc = 2; iSpc < 3; iSpc++) { // only apply to combined results.
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {

        for (TH1D* h : {a->h_z_trk_zpt[iSpc][iPtZ][iCent], a->h_z_trk_zpt_sub[iSpc][iPtZ][iCent], a->h_z_trk_zpt_iaa[iSpc][iPtZ][iCent]}) {
          if (!h)
            continue;
          TH1D* hrelsys =  relVarPt[iPtZ][iCent];
          for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
            if (upVar) {
              h->SetBinContent (ix, h->GetBinContent (ix) * (1+hrelsys->GetBinContent (ix)));
              h->SetBinError (ix, h->GetBinError (ix) * (1+hrelsys->GetBinContent (ix)));
            }
            else {
              h->SetBinContent (ix, h->GetBinContent (ix) * (1-hrelsys->GetBinContent (ix)));
              h->SetBinError (ix, h->GetBinError (ix) * (1-hrelsys->GetBinContent (ix)));
            }
          }
        }

        for (TH1D* h : {a->h_z_trk_zxzh[iSpc][iPtZ][iCent], a->h_z_trk_zxzh_sub[iSpc][iPtZ][iCent], a->h_z_trk_zxzh_iaa[iSpc][iPtZ][iCent]}) {
          if (!h)
            continue;
          TH1D* hrelsys =  relVarX[iPtZ][iCent];
          for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
            if (upVar) {
              h->SetBinContent (ix, h->GetBinContent (ix) * (1+hrelsys->GetBinContent (ix)));
              h->SetBinError (ix, h->GetBinError (ix) * (1+hrelsys->GetBinContent (ix)));
            }
            else {
              h->SetBinContent (ix, h->GetBinContent (ix) * (1-hrelsys->GetBinContent (ix)));
              h->SetBinError (ix, h->GetBinError (ix) * (1-hrelsys->GetBinContent (ix)));
            }
          }
        }
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  return;
}


#endif
