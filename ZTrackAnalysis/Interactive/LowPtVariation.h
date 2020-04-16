#ifndef __LowPtVariation_h__
#define __LowPtVariation_h__

#include "Params.h"
#include "PhysicsAnalysis.h"
#include "Systematic.h"
#include "ReweightingVariation.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class LowPtVariation : ReweightingVariation {

  public:

  LowPtVariation (const char* _name = "lowPtVariation") : ReweightingVariation (_name) { }

  virtual void GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* bkg) override;
  virtual void ApplyRelativeVariations (PhysicsAnalysis* nominal, const bool upVar = true) override;
};


void LowPtVariation :: GetRelativeVariations (PhysicsAnalysis* nominal, PhysicsAnalysis* bkg) {

  cout << "Calculating low pT variation weights from subtracted yields." << endl;

  TGAE*** g_pTch = Get2DArray <TGAE*> (2, numCentBins);
  TGAE*** g_xhZ = Get2DArray <TGAE*> (2, numCentBins);

  for (short iSpc : {0, 1}) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      g_pTch[iSpc][iCent] = new TGAE ();
      g_pTch[iSpc][iCent]->SetName (Form ("g_pTch_%s_iCent%i", iSpc == 0 ? "ee" : "mumu", iCent));
      g_xhZ[iSpc][iCent] = new TGAE ();
      g_xhZ[iSpc][iCent]->SetName (Form ("g_xhZ_%s_iCent%i", iSpc == 0 ? "ee" : "mumu", iCent));
    } // end loop over iCent
  } // end loop over iSpc

  TH1D* h_ee = nullptr, *h_mumu = nullptr, *h_comb = nullptr, *h_ee_bkg = nullptr, *h_mumu_bkg = nullptr;
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    const double wee = nominal->h_z_counts[0][2][iCent]->GetBinContent (2) / nominal->h_z_counts[2][2][iCent]->GetBinContent (2);
    const double wmm = nominal->h_z_counts[1][2][iCent]->GetBinContent (2) / nominal->h_z_counts[2][2][iCent]->GetBinContent (2);    

    {
      TGAE* g_ee = g_pTch[0][iCent];
      TGAE* g_mumu = g_pTch[1][iCent];

      h_ee = nominal->h_trk_pt_ptz_sub[0][2][iCent];
      h_mumu = nominal->h_trk_pt_ptz_sub[1][2][iCent];
      h_comb = nominal->h_trk_pt_ptz_sub[2][2][iCent];
      h_ee_bkg = bkg->h_trk_pt_ptz[0][2][iCent];
      h_mumu_bkg = bkg->h_trk_pt_ptz[1][2][iCent];
      assert (h_ee->GetNbinsX () == h_mumu->GetNbinsX () && h_ee->GetNbinsX () == h_ee_bkg->GetNbinsX () && h_ee->GetNbinsX () == h_mumu_bkg->GetNbinsX ());
      int nBinsX = h_ee->GetNbinsX ();

      for (int ix = 1; ix <= nBinsX; ix++) {
        if (h_ee->GetBinCenter (ix) < pTchBins[2][0] || pTchBins[2][nPtchBins[2]] < h_ee->GetBinCenter (ix))
          continue;

        const double yee = h_ee->GetBinContent (ix);
        const double yee_err = sqrt (pow (h_ee->GetBinError (ix), 2) + pow (h_ee_bkg->GetBinError (ix), 2));
        const double ymm = h_mumu->GetBinContent (ix);
        const double ymm_err = sqrt (pow (h_mumu->GetBinError (ix), 2) + pow (h_mumu_bkg->GetBinError (ix), 2));
        const double ycomb = h_comb->GetBinContent (ix);

        assert (yee != 0 && ymm != 0);

        const double fee = (yee-ycomb) / ycomb;
        const double fee_err = sqrt (pow (wmm*ymm*yee_err / pow (ycomb, 2), 2) + pow (-wmm*yee*ymm_err / pow (ycomb, 2), 2));

        g_ee->SetPoint (g_ee->GetN (), h_ee->GetBinCenter (ix), fee);
        g_ee->SetPointEYhigh (g_ee->GetN () - 1, fee_err);
        g_ee->SetPointEYlow (g_ee->GetN () - 1, fee_err);
        g_ee->SetPointEXhigh (g_ee->GetN () - 1, 0.5 * h_ee->GetBinWidth (ix));
        g_ee->SetPointEXlow (g_ee->GetN () - 1, 0.5 * h_ee->GetBinWidth (ix));

        const double fmm = (ymm-ycomb) / ycomb;
        const double fmm_err = sqrt (pow (wee*yee*ymm_err / pow (ycomb, 2), 2) + pow (-wee*ymm*yee_err / pow (ycomb, 2), 2));

        g_mumu->SetPoint (g_mumu->GetN (), h_mumu->GetBinCenter (ix), fmm);
        g_mumu->SetPointEYhigh (g_mumu->GetN () - 1, fmm_err);
        g_mumu->SetPointEYlow (g_mumu->GetN () - 1, fmm_err);
        g_mumu->SetPointEXhigh (g_mumu->GetN () - 1, 0.5 * h_mumu->GetBinWidth (ix));
        g_mumu->SetPointEXlow (g_mumu->GetN () - 1, 0.5 * h_mumu->GetBinWidth (ix));
      } // end loop over ix
    }

    {
      TGAE* g_ee = g_xhZ[0][iCent];
      TGAE* g_mumu = g_xhZ[1][iCent];

      h_ee = nominal->h_trk_xhz_ptz_sub[0][2][iCent];
      h_mumu = nominal->h_trk_xhz_ptz_sub[1][2][iCent];
      h_comb = nominal->h_trk_xhz_ptz_sub[2][2][iCent];
      h_ee_bkg = bkg->h_trk_xhz_ptz[0][2][iCent];
      h_mumu_bkg = bkg->h_trk_xhz_ptz[1][2][iCent];
      assert (h_ee->GetNbinsX () == h_mumu->GetNbinsX () && h_ee->GetNbinsX () == h_ee_bkg->GetNbinsX () && h_ee->GetNbinsX () == h_mumu_bkg->GetNbinsX ());
      int nBinsX = h_ee->GetNbinsX ();

      for (int ix = 1; ix <= nBinsX; ix++) {
        if (h_ee->GetBinCenter (ix) < xhZBins[2][0] || xhZBins[2][nXhZBins[2]] < h_ee->GetBinCenter (ix))
          continue;

        const double yee = h_ee->GetBinContent (ix);
        const double yee_err = sqrt (pow (h_ee->GetBinError (ix), 2) + pow (h_ee_bkg->GetBinError (ix), 2));
        const double ymm = h_mumu->GetBinContent (ix);
        const double ymm_err = sqrt (pow (h_mumu->GetBinError (ix), 2) + pow (h_mumu_bkg->GetBinError (ix), 2));
        const double ycomb = h_comb->GetBinContent (ix);

        assert (yee != 0 && ymm != 0);

        const double fee = (yee-ycomb) / ycomb;
        const double fee_err = sqrt (pow (wmm*ymm*yee_err / pow (ycomb, 2), 2) + pow (-wmm*yee*ymm_err / pow (ycomb, 2), 2));

        g_ee->SetPoint (g_ee->GetN (), h_ee->GetBinCenter (ix), fee);
        g_ee->SetPointEYhigh (g_ee->GetN () - 1, fee_err);
        g_ee->SetPointEYlow (g_ee->GetN () - 1, fee_err);
        g_ee->SetPointEXhigh (g_ee->GetN () - 1, 0.5 * h_ee->GetBinWidth (ix));
        g_ee->SetPointEXlow (g_ee->GetN () - 1, 0.5 * h_ee->GetBinWidth (ix));

        const double fmm = (ymm-ycomb) / ycomb;
        const double fmm_err = sqrt (pow (wee*yee*ymm_err / pow (ycomb, 2), 2) + pow (-wee*ymm*yee_err / pow (ycomb, 2), 2));

        g_mumu->SetPoint (g_mumu->GetN (), h_mumu->GetBinCenter (ix), fmm);
        g_mumu->SetPointEYhigh (g_mumu->GetN () - 1, fmm_err);
        g_mumu->SetPointEYlow (g_mumu->GetN () - 1, fmm_err);
        g_mumu->SetPointEXhigh (g_mumu->GetN () - 1, 0.5 * h_mumu->GetBinWidth (ix));
        g_mumu->SetPointEXlow (g_mumu->GetN () - 1, 0.5 * h_mumu->GetBinWidth (ix));
      } // end loop over ix
    }

  } // end loop over iCent

  TF1*** f_ratio_pTch = Get2DArray <TF1*> (2, numCentBins);
  TF1*** f_ratio_xhZ = Get2DArray <TF1*> (2, numCentBins);

  for (short iSpc : {0, 1}) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      if (iCent == 0) {
        f_ratio_pTch[iSpc][iCent] = new TF1 (Form ("f_ratio_pTch_%s_iCent%i", iSpc == 0 ? "ee" : "mumu", iCent), "[0]+[1]*log(x)+[2]*log(x)*log(x)", pTchBins[2][0], pTchBins[2][nPtchBins[2]]);
        f_ratio_pTch[iSpc][iCent]->SetParameter (0, 0);
        f_ratio_pTch[iSpc][iCent]->SetParameter (1, 0);
        f_ratio_pTch[iSpc][iCent]->SetParameter (2, 0);
        f_ratio_xhZ[iSpc][iCent] = new TF1 (Form ("f_ratio_xhZ_%s_iCent%i", iSpc == 0 ? "ee" : "mumu", iCent), "[0]+[1]*log(x)+[2]*log(x)*log(x)", xhZBins[2][0], xhZBins[2][nXhZBins[2]]);
        f_ratio_xhZ[iSpc][iCent]->SetParameter (0, 0);
        f_ratio_xhZ[iSpc][iCent]->SetParameter (1, 0);
        f_ratio_xhZ[iSpc][iCent]->SetParameter (2, 0);
      }
      else {
        f_ratio_pTch[iSpc][iCent] = new TF1 (Form ("f_ratio_pTch_%s_iCent%i", iSpc == 0 ? "ee" : "mumu", iCent), "[0]", pTchBins[2][0], pTchBins[2][nPtchBins[2]]);
        f_ratio_pTch[iSpc][iCent]->SetParameter (0, 0);
        f_ratio_xhZ[iSpc][iCent] = new TF1 (Form ("f_ratio_xhZ_%s_iCent%i", iSpc == 0 ? "ee" : "mumu", iCent), "[0]+[1]*log(x)", xhZBins[2][0], xhZBins[2][nXhZBins[2]]);
        f_ratio_xhZ[iSpc][iCent]->SetParameter (0, 0);
        f_ratio_xhZ[iSpc][iCent]->SetParameter (1, 0);
      }
      g_pTch[iSpc][iCent]->Fit (f_ratio_pTch[iSpc][iCent], "RN0Q");
      g_xhZ[iSpc][iCent]->Fit (f_ratio_xhZ[iSpc][iCent], "RN0Q");
    } // end loop over iCent
  } // end loop over iSpc

  TH1D* h = nullptr;

  for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      h = (TH1D*) nominal->h_trk_pt_ptz_sub[2][iPtZ][iCent]->Clone (Form ("h_relVarPt_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      relVarPt[2][iPtZ][iCent] = h;
      for (short iX = 1; iX <= h->GetNbinsX (); iX++) {
        const double rUnc = (iPtZ != 2 ? 0 : fmax (fabs (f_ratio_pTch[0][iCent]->Eval (h->GetBinCenter (iX))), fabs (f_ratio_pTch[1][iCent]->Eval (h->GetBinCenter (iX)))) / sqrt (12.));
        h->SetBinContent (iX, rUnc * h->GetBinContent (iX));
        h->SetBinError (iX, 0);
      }

      h = (TH1D*) nominal->h_trk_xhz_ptz_sub[2][iPtZ][iCent]->Clone (Form ("h_relVarX_comb_iPtZ%i_iCent%i_%s", iPtZ, iCent, name.c_str ()));
      relVarX[2][iPtZ][iCent] = h;
      for (short iX = 1; iX <= h->GetNbinsX (); iX++) {
        const double rUnc = (iPtZ != 2 ? 0 : fmax (fabs (f_ratio_xhZ[0][iCent]->Eval (h->GetBinCenter (iX))), fabs (f_ratio_xhZ[1][iCent]->Eval (h->GetBinCenter (iX)))) / sqrt (12.));
        h->SetBinContent (iX, rUnc * h->GetBinContent (iX));
        h->SetBinError (iX, 0);
      }
    } // end loop over iCent
  } // end loop over iPtZ

  //Delete2DArray (f_ratio_pTch, 2, numCentBins);
  //Delete2DArray (f_ratio_xhZ, 2, numCentBins);
  //Delete2DArray (g_pTch, 2, numCentBins);
  //Delete2DArray (g_xhZ, 2, numCentBins);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void LowPtVariation :: ApplyRelativeVariations (PhysicsAnalysis* a, const bool upVar) {
  TH1D* rVar = nullptr, *h = nullptr;
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        rVar = relVarPt[2][iPtZ][iCent];
        h = a->h_trk_pt_ptz[iSpc][iPtZ][iCent];
        if (h) {
          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * rVar->GetBinContent (iX));
          } // end loop over iX
        }
        h = a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];
        if (h) {
          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * rVar->GetBinContent (iX));
          } // end loop over iX
        }
        h = a->h_trk_pt_ptz_iaa[iSpc][iPtZ][iCent];
        if (h) {
          TH1D* hpp = a->h_trk_pt_ptz_sub[iSpc][iPtZ][0];
          TH1D* varpp = relVarPt[2][iPtZ][0];
          TH1D* hPbPb = a->h_trk_pt_ptz_sub[iSpc][iPtZ][iCent];
          TH1D* varPbPb = relVarPt[2][iPtZ][iCent];

          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            const double unc = fabs (h->GetBinContent (iX)) * sqrt (pow (varpp->GetBinContent (iX) / hpp->GetBinContent (iX), 2) + pow (varPbPb->GetBinContent (iX) / hPbPb->GetBinContent (iX), 2));
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * unc);
          } // end loop over iX
        }

        rVar = relVarX[2][iPtZ][iCent];
        h = a->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
        if (h) {
          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * rVar->GetBinContent (iX));
          } // end loop over iX
        }
        h = a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
        if (h) {
          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * rVar->GetBinContent (iX));
          } // end loop over iX
        }
        h = a->h_trk_xhz_ptz_iaa[iSpc][iPtZ][iCent];
        if (h) {
          TH1D* hpp = a->h_trk_xhz_ptz_sub[iSpc][iPtZ][0];
          TH1D* varpp = relVarX[2][iPtZ][0];
          TH1D* hPbPb = a->h_trk_xhz_ptz_sub[iSpc][iPtZ][iCent];
          TH1D* varPbPb = relVarX[2][iPtZ][iCent];

          assert (h->GetNbinsX () == rVar->GetNbinsX ());
          for (int iX = 1; iX <= h->GetNbinsX (); iX++) {
            const double unc = fabs (h->GetBinContent (iX)) * sqrt (pow (varpp->GetBinContent (iX) / hpp->GetBinContent (iX), 2) + pow (varPbPb->GetBinContent (iX) / hPbPb->GetBinContent (iX), 2));
            h->SetBinContent (iX, h->GetBinContent (iX) + (upVar ? 1 : -1) * unc);
          } // end loop over iX
        }
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc

  return;
}


#endif
