#ifndef __TrackPuritySystematic_h__
#define __TrackPuritySystematic_h__

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

class TrackPuritySystematic : public Systematic {

  TH1D***** relVarPt = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins+1, numCentBins); // iSpc, iPtZ, iPhi (or integrated at numPhiBins), iCent
  TH1D***** relVarX  = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins+1, numCentBins); // iSpc, iPtZ, iPhi (or integrated at numPhiBins), iCent

  public:

  TrackPuritySystematic (PhysicsAnalysis* nom, const char* _name = "systematics", const char* _desc = "systematic") : Systematic (nom, _name, _desc) { }


  void GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var);
  void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);

  //virtual void AddVariations (); // variations add linearly
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void TrackPuritySystematic :: GetRelativeVariations (PhysicsAnalysis* nom, PhysicsAnalysis* var) {
  TH1D* hnom = nullptr, *hvar = nullptr, *hden = nullptr;

  cout << "Calculating recommended +/- 25\% fake rate variation on total yields." << endl;
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    hnom = (TH1D*) nom->h2_num_trk_pur[iCent]->ProjectionY ("hnom");
    hden = (TH1D*) nom->h2_den_trk_pur[iCent]->ProjectionY ("hden");
    hnom->Divide (hden);
    delete hden;

    hvar = (TH1D*) var->h2_num_trk_pur[iCent]->ProjectionY ("hvar");
    hden = (TH1D*) var->h2_den_trk_pur[iCent]->ProjectionY ("hden");
    hvar->Divide (hden);
    delete hden;

    for (int ix = 1; ix <= hvar->GetNbinsX (); ix++) {
      const float fakeRate = 1.-hvar->GetBinContent (ix);
      fakeRate *= 1.25;
      hvar->SetBinContent (ix, fakeRate);
    }

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {

          relVarPt[iSpc][iPtZ][iPhi][iCent] = (TH1D*) hnom->Clone (Form "h_relVarPt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ());
          relVarPt[iSpc][iPtZ][iPhi][iCent] = 
          var->h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]->Multiply (h
          

          h = (TH1D*) nominal

          h = (TH1D*) nominal->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]-

          h = (TH1D*) nominal->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Clone (Form ("h_relVarX_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()));
          relVarX[iSpc][iPtZ][iPhi][iCent] = h;
          h->Divide (var->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
        } // end loop over iPhi

        h = (TH1D*) nominal->h_z_trk_zpt[iSpc][iPtZ][iCent]->Clone (Form ("h_relVarPt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, numPhiBins, iCent, name.c_str ()));
        relVarPt[iSpc][iPtZ][numPhiBins][iCent] = h;
        h->Divide (var->h_z_trk_zpt[iSpc][iPtZ][iCent]);

        h = (TH1D*) nominal->h_z_trk_zxzh[iSpc][iPtZ][iCent]->Clone (Form ("h_relVarX_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, numPhiBins, iCent, name.c_str ()));
        relVarX[iSpc][iPtZ][numPhiBins][iCent] = h;
        h->Divide (var->h_z_trk_zxzh[iSpc][iPtZ][iCent]);
      } // end loop over iPtZ
    } // end loop over iSpc
  } // end loop over iCent
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void TrackPuritySystematic :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

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
        } // end loop over iPhi
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

