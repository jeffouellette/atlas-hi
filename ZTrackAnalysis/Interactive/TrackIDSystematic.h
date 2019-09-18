#ifndef __TrackIDSystematic_h__
#define __TrackIDSystematic_h__

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

class TrackIDSystematic : public Systematic {

  public:

  float**** relVar = Get4DArray <float> (3, nPtZBins, numPhiBins+1, numCentBins); // iSpc, iPtZ, iPhi (or integrated at numPhiBins), iCent
  TrackIDSystematic (PhysicsAnalysis* nom, const char* _name = "systematics", const char* _desc = "systematic") : Systematic (nom, _name, _desc) { }


  void GetRelativeVariation (PhysicsAnalysis* nominal, PhysicsAnalysis* var);

  //virtual void AddVariations (); // variations add linearly
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Calculates the relative variation between the nominal and variation results on the 
// raw (unscaled) hadron yields.
////////////////////////////////////////////////////////////////////////////////////////////////
void TrackIDSystematic :: GetRelativeVariation (PhysicsAnalysis* nominal, PhysicsAnalysis* var) {
  TH1D* hnom = nullptr, *hvar = nullptr, *hrat = nullptr, *eff1 = nullptr, *eff2 = nullptr, *effrat = nullptr, *purrat = nullptr;

  cout << "Calculating track ID variation on total yields." << endl;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) { 

      for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          //cout << "Getting relative variation for iSpc=" << iSpc << ", iPtZ=" << iPtZ << ", iPhi=" << iPhi << ", iCent=" << iCent << endl;

          hnom = nominal->h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent];
          hvar = var->h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent];

          hrat = (TH1D*) hvar->Clone ("rattemp");
          hrat->Divide (hnom);
          for (int ix = 1; ix <= hrat->GetNbinsX (); ix++)
            hrat->SetBinError (ix, sqrt ((hrat->GetBinContent (ix)) * (1-hrat->GetBinContent (ix)) / hnom->GetBinContent (ix)));

          nominal->LoadTrackingEfficiencies (true);
          var->LoadTrackingEfficiencies (true);
          nominal->LoadTrackingPurities (true);
          var->LoadTrackingPurities (true);

          eff1 = (TH1D*) nominal->h2_num_trk_effs[iCent]->ProjectionY ("1");
          effrat = (TH1D*) nominal->h2_den_trk_effs[iCent]->ProjectionY ("2");
          eff1->Divide (effrat);
          delete effrat;
          eff2 = (TH1D*) var->h2_num_trk_effs[iCent]->ProjectionY ("3");
          effrat = (TH1D*) var->h2_den_trk_effs[iCent]->ProjectionY ("4");
          eff2->Divide (effrat);
          delete effrat;

          effrat = (TH1D*) eff2->Clone ("efftemp");
          effrat->Divide (eff1);
          delete eff1, eff2;

          eff1 = (TH1D*) nominal->h2_num_trk_purs[iCent]->ProjectionY ("1");
          eff1->Divide ((TH1D*) nominal->h2_den_trk_purs[iCent]->ProjectionY ("2"));
          eff2 = (TH1D*) var->h2_num_trk_purs[iCent]->ProjectionY ("3");
          eff2->Divide ((TH1D*) var->h2_den_trk_purs[iCent]->ProjectionY ("4"));

          purrat = (TH1D*) eff2->Clone ("purtemp");
          purrat->Divide (eff1);
          delete eff1, eff2;

          const int bin1 = effrat->FindBin (hrat->GetBinCenter (1));
          const int bin2 = purrat->FindBin (hrat->GetBinCenter (1));
          for (int iy = 1; iy <= hrat->GetNbinsX (); iy++)
            hrat->SetBinContent (iy, hrat->GetBinContent (iy) * purrat->GetBinContent (iy+bin2-1) / effrat->GetBinContent (iy+bin1-1));

          TF1* fit = new TF1 ("fit", "[0]", allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
          fit->SetParameter (0, 1);
          hrat->Fit (fit, "RQN0");
          relVar[iSpc][iPtZ][iPhi][iCent] = fit->GetParameter (0);
          delete hrat, effrat, purrat, fit;
        } // end loop over iCent
      } // end loop over iPhi

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        hnom = (TH1D*) nominal->h_z_trk_raw_pt[iSpc][iPtZ][1][iCent]->Clone ("hnom");
        hvar = (TH1D*) var->h_z_trk_raw_pt[iSpc][iPtZ][1][iCent]->Clone ("hvar");
        hnom->Add (nominal->h_z_trk_raw_pt[iSpc][iPtZ][2][iCent]);
        hvar->Add (var->h_z_trk_raw_pt[iSpc][iPtZ][2][iCent]);

        hrat = (TH1D*) hvar->Clone ("rattemp");
        hrat->Divide (hnom);
        for (int ix = 1; ix <= hrat->GetNbinsX (); ix++)
          hrat->SetBinError (ix, sqrt ((hrat->GetBinContent (ix)) * (1-hrat->GetBinContent (ix)) / hnom->GetBinContent (ix)));

        nominal->LoadTrackingEfficiencies (true);
        var->LoadTrackingEfficiencies (true);
        nominal->LoadTrackingPurities (true);
        var->LoadTrackingPurities (true);

        eff1 = (TH1D*) nominal->h2_num_trk_effs[iCent]->ProjectionY ("1");
        effrat = (TH1D*) nominal->h2_den_trk_effs[iCent]->ProjectionY ("2");
        eff1->Divide (effrat);
        delete effrat;
        eff2 = (TH1D*) var->h2_num_trk_effs[iCent]->ProjectionY ("3");
        effrat = (TH1D*) var->h2_den_trk_effs[iCent]->ProjectionY ("4");
        eff2->Divide (effrat);
        delete effrat;

        effrat = (TH1D*) eff2->Clone ("efftemp");
        effrat->Divide (eff1);
        delete eff1, eff2;

        eff1 = (TH1D*) nominal->h2_num_trk_purs[iCent]->ProjectionY ("1");
        eff1->Divide ((TH1D*) nominal->h2_den_trk_purs[iCent]->ProjectionY ("2"));
        eff2 = (TH1D*) var->h2_num_trk_purs[iCent]->ProjectionY ("3");
        eff2->Divide ((TH1D*) var->h2_den_trk_purs[iCent]->ProjectionY ("4"));

        purrat = (TH1D*) eff2->Clone ("purtemp");
        purrat->Divide (eff1);
        delete eff1, eff2;

        const int bin1 = effrat->FindBin (hrat->GetBinCenter (1));
        const int bin2 = purrat->FindBin (hrat->GetBinCenter (1));
        for (int iy = 1; iy <= hrat->GetNbinsX (); iy++)
          hrat->SetBinContent (iy, hrat->GetBinContent (iy) * purrat->GetBinContent (iy+bin2-1) / effrat->GetBinContent (iy+bin1-1));

        TF1* fit = new TF1 ("fit", "[0]", allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
        fit->SetParameter (0, 1);
        hrat->Fit (fit, "RQN0");
        relVar[iSpc][iPtZ][numPhiBins][iCent] = fit->GetParameter (0);
        delete hrat, effrat, purrat, fit;
        delete hnom, hvar;
      } // end loop over iCent
    } // end loop over iPtZ
  } // end loop over iSpc
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Sets the errors in this systematic as a combination of all added variations.
// Takes the maximum error for each point.
// Intended for combining up & down variations, but expandable for additional categories
// (e.g. track quality criteria)
////////////////////////////////////////////////////////////////////////////////////////////////
/*void TrackIDSystematic :: AddVariations () {

  NullifyErrors ();

  for (PhysicsAnalysis* a : variations) {

    cout << "Adding variation " << a->Name () << " to systematic " << name << endl;

    //a->ApplyRelativeVariation (relVars, variationDirs[a]);
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
            if (sys && var) CalcSystematics (sys, var, true);

            sys = GetTGAE (h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, true);
            sys = GetTGAE (h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, true);
            sys = GetTGAE (h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, true);

            sys = GetTGAE (h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, true);
            sys = GetTGAE (h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, true);
            sys = GetTGAE (h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent]);
            var = a->h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var, true);

          } // end loop over cents
        } // end loop over phi

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            sys = GetTGAE (h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]);
            var = a->h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent];
            if (sys && var) CalcSystematics (sys, var, true);
            sys = GetTGAE (h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent]);
            var = a->h_z_trk_phi_sub[iSpc][iPtZ][iPtTrk][iCent];
            if (sys && var) CalcSystematics (sys, var, true);
          } // end loop over cents
        } // end loop over iPtTrk

        for (short iCent = 0; iCent < numCentBins; iCent++) {
          sys = GetTGAE (h_z_trk_zpt[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zpt[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, true);
          sys = GetTGAE (h_z_trk_zpt_sub[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zpt_sub[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, true);
          sys = GetTGAE (h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zpt_sig_to_bkg[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, true);

          sys = GetTGAE (h_z_trk_zxzh[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zxzh[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, true);
          sys = GetTGAE (h_z_trk_zxzh_sub[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zxzh_sub[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, true);
          sys = GetTGAE (h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent]);
          var = a->h_z_trk_zxzh_sig_to_bkg[iSpc][iPtZ][iCent];
          if (sys && var) CalcSystematics (sys, var, true);
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
}*/

   
#endif
