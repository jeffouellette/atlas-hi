#ifndef __ReweightingVariation_h__
#define __ReweightingVariation_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class ReweightingVariation {

  protected:
  string name;

  public:
  TH1D**** relVarPt = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D**** relVarX  = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent

  //ReweightingVariation (PhysicsAnalysis* nom, const char* _name = "systematics", const char* _desc = "systematic") : Systematic (nom, _name, _desc) { }
  ReweightingVariation (const char* _name = "reweightingVariation") {
    name = _name;
  }

  ~ReweightingVariation () {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          SaferDelete (&relVarPt[iSpc][iPtZ][iCent]);
          SaferDelete (&relVarX[iSpc][iPtZ][iCent]);
        } // end loop over iCent
      } // end loop over iPtZ
    } // end loop over iSpc
    Delete3DArray (&relVarPt, 3, nPtZBins, numCentBins);
    Delete3DArray (&relVarX, 3, nPtZBins, numCentBins);
  }

  string Name ()              { return name; }
  void SetName (string _name) { name = _name; }

  void ScaleRelativeVariations (const float scaleFactor = 1.);

  virtual void GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var);
  virtual void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Rescales the relative variation weights by some constant factor
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingVariation :: ScaleRelativeVariations (const float scaleFactor) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
      for (int iPhi = 0; iPhi <= numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (TH1D* h : {relVarPt[iSpc][iPtZ][iCent], relVarX[iSpc][iPtZ][iCent]}) {
            if (!h)
              continue;
            for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
              float relVar = h->GetBinContent (ix);
              relVar = (relVar-1.)*scaleFactor + 1.;
              h->SetBinContent (ix, relVar);
            }
          } // end loop over h
        } // end loop over iCent
      } // end loop over iPhi
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingVariation :: GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var) {
  TH1D* h = nullptr;

  cout << "Calculating reweighting variation on total yields." << endl;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        h = (TH1D*) nominal->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_relVarPt_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        relVarPt[iSpc][iPtZ][iCent] = h;
        h->Divide (var->h_trk_pt_ptz[iSpc][iPtZ][iCent]);

        h = (TH1D*) nominal->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Clone (Form ("h_relVarX_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        relVarX[iSpc][iPtZ][iCent] = h;
        h->Divide (var->h_trk_xhz_ptz[iSpc][iPtZ][iCent]);
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingVariation :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        if (upVar) {
          a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][iCent]);
          a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][iCent]);
          //a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][iCent]);
          //a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][iCent]);
          //a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][iCent]);
          //a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][iCent]);
        }
        else {
          a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][iCent]);
          a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][iCent]);
          //a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][iCent]);
          //a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][iCent]);
          //a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][iCent]);
          //a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][iCent]);
        }
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  return;
}


#endif
