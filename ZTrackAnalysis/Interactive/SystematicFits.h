#ifndef __SystematicFits_h__
#define __SystematicFits_h__

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

class SystematicFits {

  private:
  string name;

  TF1***** relVarPt = Get4DArray <TF1*> (3, nPtZBins, numPhiBins+1, numCentBins); // iSpc, iPtZ, iPhi (or integrated at numPhiBins), iCent
  TF1***** relVarX  = Get4DArray <TF1*> (3, nPtZBins, numPhiBins+1, numCentBins); // iSpc, iPtZ, iPhi (or integrated at numPhiBins), iCent

  public:

  SystematicFits (const char* _name = "systematics") {
    name = _name;
  }

  ~SystematicFits () {
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

  void GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var);
  void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);

  //virtual void AddVariations (); // variations add linearly
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void SystematicFits :: GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* var) {
  TH1D* h = nullptr;
  TF1* f = nullptr;

  cout << "Calculating variation fits on total yields." << endl;

  for (short iSpc = 2; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
      //  for (short iCent = 0; iCent < numCentBins; iCent++) {
      //    h = (TH1D*) var->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Clone ("temp");
      //    h->Divide (nominal->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
      //    f = new TF1 (Form ("f_relVarPt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "[0]+[1]*log(x)", ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]);
      //    f->SetParameter (0, 1);
      //    f->SetParameter (1, 0);
      //    h->Fit (f, "RN0Q");
      //    delete h;
      //    relVarPt[iSpc][iPtZ][iPhi][iCent] = f;
  
      //    h = (TH1D*) var->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Clone ("temp");
      //    h->Divide (nominal->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
      //    f = new TF1 (Form ("f_relVarX_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, iPhi, iCent, name.c_str ()), "[0]+[1]*log(x)", xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[iPtZ]]);
      //    f->SetParameter (0, 1);
      //    f->SetParameter (1, 0);
      //    h->Fit (f, "RN0Q");
      //    delete h;
      //    relVarX[iSpc][iPtZ][iPhi][iCent] = f;
      //  } // end loop over iCent
      //} // end loop over iPhi

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        h = (TH1D*) var->h_z_trk_zpt[iSpc][iPtZ][iCent]->Clone ("temp");
        h->Divide (nominal->h_z_trk_zpt[iSpc][iPtZ][iCent]);
        f = new TF1 (Form ("f_relVarPt_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, numPhiBins, iCent, name.c_str ()), "[0]+[1]*log(x)", ptTrkBins[iPtZ][0], ptTrkBins[iPtZ][nPtTrkBins[iPtZ]]);
        f->SetParameter (0, 1);
        f->SetParameter (1, 0);
        h->Fit (f, "RN0Q");
        delete h;
        relVarPt[iSpc][iPtZ][numPhiBins][iCent] = f;

        h = (TH1D*) var->h_z_trk_zxzh[iSpc][iPtZ][iCent]->Clone ("temp");
        h->Divide (nominal->h_z_trk_zxzh[iSpc][iPtZ][iCent]);
        f = new TF1 (Form ("f_relVarX_%s_iPtZ%i_iPhi%i_iCent%i_%s", spc, iPtZ, numPhiBins, iCent, name.c_str ()), "[0]+[1]*log(x)", xHZBins[iPtZ][0], xHZBins[iPtZ][nXHZBins[iPtZ]]);
        f->SetParameter (0, 1);
        f->SetParameter (1, 0);
        h->Fit (f, "RN0Q");
        delete h;
        relVarX[iSpc][iPtZ][numPhiBins][iCent] = f;
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void SystematicFits :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  for (short iSpc = 2; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        //for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {

        //  for (TH1D* h : {a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent], a->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]}) {//, a->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]}) {
        //    if (!h)
        //      continue;
        //    TF1* f =  relVarPt[iSpc][iPtZ][numPhiBins][iCent];
        //    for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
        //      if (upVar) {
        //        h->SetBinContent (ix, h->GetBinContent (ix) * f->Eval (h->GetBinCenter (ix)));
        //        h->SetBinError (ix, h->GetBinError (ix) * f->Eval (h->GetBinCenter (ix)));
        //      }
        //      else {
        //        h->SetBinContent (ix, h->GetBinContent (ix) * (1./f->GetParameter (0)) * (1. - (f->GetParameter(1)/f->GetParameter(0)) * TMath::Log ((h->GetBinCenter (ix)))));
        //        h->SetBinError (ix, h->GetBinError (ix) * (1./f->GetParameter (0)) * (1. - (f->GetParameter(1)/f->GetParameter(0)) * TMath::Log ((h->GetBinCenter (ix)))));
        //      }
        //    }
        //  }

        //  for (TH1D* h : {a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent], a->h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]}) {//, a->h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]}) {
        //    if (!h)
        //      continue;
        //    TF1* f =  relVarX[iSpc][iPtZ][numPhiBins][iCent];
        //    for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
        //      if (upVar) {
        //        h->SetBinContent (ix, h->GetBinContent (ix) * f->Eval (h->GetBinCenter (ix)));
        //        h->SetBinError (ix, h->GetBinError (ix) * f->Eval (h->GetBinCenter (ix)));
        //      }
        //      else {
        //        h->SetBinContent (ix, h->GetBinContent (ix) * (1./f->GetParameter (0)) * (1. - (f->GetParameter(1)/f->GetParameter(0)) * TMath::Log ((h->GetBinCenter (ix)))));
        //        h->SetBinError (ix, h->GetBinError (ix) * (1./f->GetParameter (0)) * (1. - (f->GetParameter(1)/f->GetParameter(0)) * TMath::Log ((h->GetBinCenter (ix)))));
        //      }
        //    }
        //  }
        //} // end loop over iPhi

        for (TH1D* h : {a->h_z_trk_zpt[iSpc][iPtZ][iCent], a->h_z_trk_zpt_sub[iSpc][iPtZ][iCent]}) {//, a->h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]}) {
          if (!h)
            continue;
          TF1* f =  relVarPt[iSpc][iPtZ][numPhiBins][iCent];
          for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
            if (upVar) {
              h->SetBinContent (ix, h->GetBinContent (ix) * f->Eval (h->GetBinCenter (ix)));
              h->SetBinError (ix, h->GetBinError (ix) * f->Eval (h->GetBinCenter (ix)));
            }
            else {
              h->SetBinContent (ix, h->GetBinContent (ix) * (1./f->GetParameter (0)) * (1. - (f->GetParameter(1)/f->GetParameter(0)) * TMath::Log ((h->GetBinCenter (ix)))));
              h->SetBinError (ix, h->GetBinError (ix) * (1./f->GetParameter (0)) * (1. - (f->GetParameter(1)/f->GetParameter(0)) * TMath::Log ((h->GetBinCenter (ix)))));
            }
          }
        }

        for (TH1D* h : {a->h_z_trk_zxzh[iSpc][iPtZ][iCent], a->h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]}) {//, a->h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]}) {
          if (!h)
            continue;
          TF1* f =  relVarX[iSpc][iPtZ][numPhiBins][iCent];
          for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
            if (upVar) {
              h->SetBinContent (ix, h->GetBinContent (ix) * f->Eval (h->GetBinCenter (ix)));
              h->SetBinError (ix, h->GetBinError (ix) * f->Eval (h->GetBinCenter (ix)));
            }
            else {
              h->SetBinContent (ix, h->GetBinContent (ix) * (1./f->GetParameter (0)) * (1. - (f->GetParameter(1)/f->GetParameter(0)) * TMath::Log ((h->GetBinCenter (ix)))));
              h->SetBinError (ix, h->GetBinError (ix) * (1./f->GetParameter (0)) * (1. - (f->GetParameter(1)/f->GetParameter(0)) * TMath::Log ((h->GetBinCenter (ix)))));
            }
          }
        }
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  return;
}


#endif
