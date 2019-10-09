#ifndef __ReweightingSystematic_h__
#define __ReweightingSystematic_h__

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

class ReweightingSystematic {

  protected:
  string name;

  TH1D***** relVarPt = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins+1, numCentBins); // iSpc, iPtZ, iPhi (or integrated at numPhiBins), iCent
  TH1D***** relVarX  = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins+1, numCentBins); // iSpc, iPtZ, iPhi (or integrated at numPhiBins), iCent

  public:

  //ReweightingSystematic (PhysicsAnalysis* nom, const char* _name = "systematics", const char* _desc = "systematic") : Systematic (nom, _name, _desc) { }
  ReweightingSystematic (const char* _name = "systematics") {
    name = _name;
  }

  ~ReweightingSystematic () {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 
        for (int iPhi = 0; iPhi <= numPhiBins; iPhi++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            SaferDelete (relVarPt[iSpc][iPtZ][iPhi][iCent]);
            SaferDelete (relVarX[iSpc][iPtZ][iPhi][iCent]);
          }
        }
      }
    }
    Delete4DArray (relVarPt, 3, nPtZBins, numPhiBins+1, numCentBins);
    Delete4DArray (relVarX, 3, nPtZBins, numPhiBins+1, numCentBins);
  }

  string Name ()              { return name; }
  void SetName (string _name) { name = _name; }

  void ScaleRelativeVariations (const float scaleFactor = 1.);

  void GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var);
  void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Rescales the relative variation weights by some constant factor
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingSystematic :: ScaleRelativeVariations (const float scaleFactor) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
      for (int iPhi = 0; iPhi <= numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (TH1D* h : {relVarPt[iSpc][iPtZ][iPhi][iCent], relVarX[iSpc][iPtZ][iPhi][iCent]}) {
            if (!h)
              continue;
            for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
              float relVar = h->GetBinContent (ix);
              relVar = (relVar-1.)*scaleFactor + 1.;
              h->SetBinContent (ix, relVar);
            }
          }
        }
      }
    }
  }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingSystematic :: GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var) {
  TH1D* h = nullptr, *eff1 = nullptr, *eff2 = nullptr, *effrat = nullptr;

  cout << "Calculating reweighting variation on total yields." << endl;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
      //  for (short iCent = 0; iCent < numCentBins; iCent++) {
      //    h = (TH1D*) nominal->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_relVarPt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
      //    relVarPt[iSpc][iPtZ][iPhi][iCent] = h;
      //    h->Divide (var->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);

      //    h = (TH1D*) nominal->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_relVarX_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
      //    relVarX[iSpc][iPtZ][iPhi][iCent] = h;
      //    h->Divide (var->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
      //  } // end loop over iCent
      //} // end loop over iPhi

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        h = (TH1D*) nominal->h_z_trk_zpt[iSpc][iPtZ][iCent]->Clone (Form ("h_relVarPt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, numPhiBins, iCent, name.c_str ()));
        relVarPt[iSpc][iPtZ][numPhiBins][iCent] = h;
        h->Divide (var->h_z_trk_zpt[iSpc][iPtZ][iCent]);

        h = (TH1D*) nominal->h_z_trk_zxzh[iSpc][iPtZ][iCent]->Clone (Form ("h_relVarX_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, numPhiBins, iCent, name.c_str ()));
        relVarX[iSpc][iPtZ][numPhiBins][iCent] = h;
        h->Divide (var->h_z_trk_zxzh[iSpc][iPtZ][iCent]);
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  //for (short iSpc = 0; iSpc < 2; iSpc++) {
  //  const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
  //  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

  //    //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
  //    //  for (short iCent = 0; iCent < numCentBins; iCent++) {
  //    //    relVarPt[iSpc][iPtZ][iPhi][iCent] = relVarPt[2][iPtZ][iPhi][iCent];
  //    //    relVarX[iSpc][iPtZ][iPhi][iCent] = relVarX[2][iPtZ][iPhi][iCent];
  //    //  } // end loop over iCent
  //    //} // end loop over iPhi

  //    for (short iCent = 0; iCent < numCentBins; iCent++) {
  //      relVarPt[iSpc][iPtZ][numPhiBins][iCent] = relVarPt[2][iPtZ][numPhiBins][iCent];
  //      relVarX[iSpc][iPtZ][numPhiBins][iCent] = relVarX[2][iPtZ][numPhiBins][iCent];
  //    } // end loop over iCent
  //  } // end loop over iPtZ
  //} // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingSystematic :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        //  if (upVar) {
        //    a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Multiply (relVarPt[iSpc][iPtZ][iPhi][iCent]);
        //    a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Multiply (relVarX[iSpc][iPtZ][iPhi][iCent]);
        //  }
        //  else {
        //    a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Divide (relVarPt[iSpc][iPtZ][iPhi][iCent]);
        //    a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Divide (relVarX[iSpc][iPtZ][iPhi][iCent]);
        //  }
        //} // end loop over iPhi
        if (upVar) {
          a->h_z_trk_zpt[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][numPhiBins][iCent]);
          a->h_z_trk_zxzh[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][numPhiBins][iCent]);
        }
        else {
          a->h_z_trk_zpt[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][numPhiBins][iCent]);
          a->h_z_trk_zxzh[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][numPhiBins][iCent]);
        }
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  return;
}


#endif
