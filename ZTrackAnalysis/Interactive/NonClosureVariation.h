#ifndef __NonClosureVariation_h__
#define __NonClosureVariation_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class NonClosureVariation {

  protected:
  string name;

  public:
  TH1D**** relVarPt = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent
  TH1D**** relVarX  = Get3DArray <TH1D*> (3, nPtZBins, numCentBins); // iSpc, iPtZ, iCent

  NonClosureVariation (const char* _name = "nonClosureVariation") {
    name = _name;
  }

  ~NonClosureVariation () {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          SaferDelete (&relVarPt[iSpc][iPtZ][iCent]);
          SaferDelete (&relVarX[iSpc][iPtZ][iCent]);
        } // end loop over iCent
      } // end loop over iPtZ
    } // end loop over iSpc
    Delete3DArray (&relVarPt, 3, nPtZBins, numCentBins);
    Delete3DArray (&relVarX, 3, nPtZBins, numCentBins);
  }

  virtual void GetRelativeVariations (PhysicsAnalysis* nominal);
  virtual void ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar = true);
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Gets the relative variations from the MC closure study as a function of S/B ratio.
////////////////////////////////////////////////////////////////////////////////////////////////
void NonClosureVariation :: GetRelativeVariations (PhysicsAnalysis* nominal) {

  nominal->CalculateSigToBkg ();
 
  SetupDirectories ("", "ZTrackAnalysis/");

  TDirectory* _gDirectory = gDirectory;

  TFile* f_pTch = new TFile (Form ("%s/ClosureSystematic/closureSystematic_pTch.root", rootPath.Data ()), "read");
  TF1* f_closure_sigToBkg_pTch = (TF1*) f_pTch->Get ("f_closure_sigToBkg_pTch");
  TFile* f_xhZ = new TFile (Form ("%s/ClosureSystematic/closureSystematic_xhZ.root", rootPath.Data ()), "read");
  TF1* f_closure_sigToBkg_xhZ = (TF1*) f_xhZ->Get ("f_closure_sigToBkg_xhZ");

  _gDirectory->cd ();

  TH1D* h = nullptr, *s2b = nullptr;

  cout << "Calculating non-closure variations on yields." << endl;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        h = (TH1D*) nominal->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_nonClosureVarPt_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        relVarPt[iSpc][iPtZ][iCent] = h;
        s2b = nominal->h_trk_pt_ptz_sig_to_bkg[iSpc][iPtZ][iCent];

        for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
          assert (h->GetBinCenter (iX) == s2b->GetBinCenter (iX));
          float rerr = 1;
          if (iCent == 0) rerr = 1.02; // constant 2% error in pp
          else if (!isnan (s2b->GetBinContent (iX)) && s2b->GetBinContent (iX) > 0)
            rerr = f_closure_sigToBkg_pTch->Eval (s2b->GetBinContent (iX));
          h->SetBinContent (iX, 1.+0.5*fabs(rerr-1.));
          h->SetBinError (iX, 0);
        }
        //cout << relVarPt[iSpc][iPtZ][iCent]->GetNbinsX () << endl;

        h = (TH1D*) nominal->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Clone (Form ("h_nonClosureVarX_%s_iPtZ%i_iCent%i_%s", spc, iPtZ, iCent, name.c_str ()));
        relVarX[iSpc][iPtZ][iCent] = h;
        s2b = nominal->h_trk_xhz_ptz_sig_to_bkg[iSpc][iPtZ][iCent];

        for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
          assert (h->GetBinCenter (iX) == s2b->GetBinCenter (iX));
          float rerr = 1;
          if (iCent == 0) rerr = 1.02; // constant 2% error in pp
          else if (!isnan (s2b->GetBinContent (iX)) && s2b->GetBinContent (iX) > 0)
            rerr = f_closure_sigToBkg_xhZ->Eval (s2b->GetBinContent (iX));
          h->SetBinContent (iX, 1.+0.5*fabs(rerr-1.));
          h->SetBinError (iX, 0);
        }
        //cout << relVarX[iSpc][iPtZ][iCent]->GetNbinsX () << endl;
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  f_pTch->Close ();
  f_xhZ->Close ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void NonClosureVariation :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        if (upVar) {
          if (iCent == 0) {
            a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][iCent]);
            //TH1D* h = a->h_trk_pt_ptz[iSpc][iPtZ][iCent];
            //for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            //  assert (h->GetBinCenter (iX) == relVarPt[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            //  h->SetBinContent (iX, h->GetBinContent (iX) * relVarPt[iSpc][iPtZ][iCent]->GetBinContent (iX));
            //}

            a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][iCent]);
            //h = a->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
            //for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            //  assert (h->GetBinCenter (iX) == relVarX[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            //  h->SetBinContent (iX, h->GetBinContent (iX) * relVarX[iSpc][iPtZ][iCent]->GetBinContent (iX));
            //}
          }
          //a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][iCent]);
          TH1D* h = a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            assert (h->GetBinCenter (iX) == relVarPt[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) * relVarPt[iSpc][iPtZ][iCent]->GetBinContent (iX));
          }

          //a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][iCent]);
          h = a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            assert (h->GetBinCenter (iX) == relVarX[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) * relVarX[iSpc][iPtZ][iCent]->GetBinContent (iX));
          }
          //a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Multiply (relVarPt[iSpc][iPtZ][iCent]);
          //a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Multiply (relVarX[iSpc][iPtZ][iCent]);
        }
        else {
          if (iCent == 0) {
            //a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][iCent]);
            TH1D* h = a->h_trk_pt_ptz[iSpc][iPtZ][iCent];
            for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
              assert (h->GetBinCenter (iX) == relVarPt[iSpc][iPtZ][iCent]->GetBinCenter (iX));
              h->SetBinContent (iX, h->GetBinContent (iX) / relVarPt[iSpc][iPtZ][iCent]->GetBinContent (iX));
            }

            //a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][iCent]);
            h = a->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
            for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
              assert (h->GetBinCenter (iX) == relVarX[iSpc][iPtZ][iCent]->GetBinCenter (iX));
              h->SetBinContent (iX, h->GetBinContent (iX) / relVarX[iSpc][iPtZ][iCent]->GetBinContent (iX));
            }
          }

          a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][iCent]);
          TH1D* h = a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            assert (h->GetBinCenter (iX) == relVarPt[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) / relVarPt[iSpc][iPtZ][iCent]->GetBinContent (iX));
          }

          //a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][iCent]);
          h = a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            assert (h->GetBinCenter (iX) == relVarX[iSpc][iPtZ][iCent]->GetBinCenter (iX));
            h->SetBinContent (iX, h->GetBinContent (iX) / relVarX[iSpc][iPtZ][iCent]->GetBinContent (iX));
          }
          //a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent]->Divide (relVarPt[iSpc][iPtZ][iCent]);
          //a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent]->Divide (relVarX[iSpc][iPtZ][iCent]);
        }
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  return;
}


#endif
