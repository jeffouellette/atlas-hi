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

class ReweightingSystematic : public Systematic {

  TH1D***** relVarPt = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins+1, numCentBins); // iSpc, iPtZ, iPhi (or integrated at numPhiBins), iCent
  TH1D***** relVarX  = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins+1, numCentBins); // iSpc, iPtZ, iPhi (or integrated at numPhiBins), iCent

  public:

  ReweightingSystematic (PhysicsAnalysis* nom, const char* _name = "systematics", const char* _desc = "systematic") : Systematic (nom, _name, _desc) { }


  void GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var);
  void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);

  //virtual void AddVariations (); // variations add linearly
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingSystematic :: GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var) {
  TH1D* h = nullptr, *eff1 = nullptr, *eff2 = nullptr, *effrat = nullptr;

  cout << "Calculating particle composition variation on total yields." << endl;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) { 

      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          h = (TH1D*) nominal->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_pt_partcomp_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          relVarPt[iSpc][iPtZ][iPhi][iCent] = h;
          h->Divide (var->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);

          h = (TH1D*) nominal->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_z_trk_xzh_partcomp_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          relVarX[iSpc][iPtZ][iPhi][iCent] = h;
          h->Divide (var->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over iCent
      } // end loop over iPhi

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        h = (TH1D*) nominal->h_z_trk_zpt[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zpt_partcomp_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, numPhiBins, iCent, name.c_str ()));
        relVarPt[iSpc][iPtZ][numPhiBins][iCent] = h;
        h->Divide (var->h_z_trk_zpt[iSpc][iPtZ][iCent]);

        h = (TH1D*) nominal->h_z_trk_zxzh[iSpc][iPtZ][iCent]->Clone (Form ("h_z_trk_zxzh_partcomp_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, numPhiBins, iCent, name.c_str ()));
        relVarX[iSpc][iPtZ][numPhiBins][iCent] = h;
        h->Divide (var->h_z_trk_zxzh[iSpc][iPtZ][iCent]);
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void ReweightingSystematic :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          if (upVar) {
            a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Multiply (relVarPt[iSpc][iPtZ][iPhi][iCent]);
            a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Multiply (relVarX[iSpc][iPtZ][iPhi][iCent]);
          }
          else {
            a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Divide (relVarPt[iSpc][iPtZ][iPhi][iCent]);
            a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Divide (relVarX[iSpc][iPtZ][iPhi][iCent]);
          }
        }
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
