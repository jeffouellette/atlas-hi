#ifndef __Systematics_h__
#define __Systematics_h__

#include "Params.h"
#include "Analysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class Systematics : public Analysis {

  private:
  vector<Analysis*> variations;

  void NullifyErrors ();

  public:
  Systematics (Analysis* nom) : Analysis (){
    name = "systematics";
    directory = "Systematics/";
    SetupDirectories (directory, "ZTrackAnalysis/");

    CopyAnalysis (nom, true);
    NullifyErrors ();
  }

  vector<Analysis*>& GetVariations () { return variations; }
  void AddVariation (Analysis* a);
  void CombineErrors ();
  void AddSystematic (Systematics* s);

  void PlotTrkYieldSystematics (const short pSpc = 2, const short pPtZ = nPtZBins-1);

};


////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: AddVariation (Analysis* a) {
  variations.push_back (a);
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Sets the errors in this systematic as a combination of all added variations.
// Takes the maximum error for each point.
// Intended for combining up & down variations, but expandable for additional categories
// (e.g. track quality criteria)
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: CombineErrors () {

  NullifyErrors ();

  for (Analysis* a : variations) {

    a->SubtractBackground ();
    a->CalculateIAA ();
    a->CalculateICP ();

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

        // Hadron yield systematics, signal & signal+bkg levels
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
              TH1D* sys = h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              TH1D* var = a->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              if (sys && var) CalcSystematics (sys, var);

              
              sys = h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              var = a->h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              if (sys && var) CalcSystematics (sys, var);

              sys = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              var = a->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              if (sys && var) CalcSystematics (sys, var);
            } // end loop over xztrk bins
          } // end loop over cents
        } // end loop over phi

        // IAA, ICP systematics
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 1; iCent < numCentBins; iCent++) {
            TH1D* sys = h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];
            TH1D* var = a->h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
          } // end loop over cents
          for (short iCent = 2; iCent < numCentBins; iCent++) {
            TH1D* sys = h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];
            TH1D* var = a->h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
          } // end loop over cents
        } // end loop over phi

      } // end loop over pT^Z bins
    } // end loop over species

  } // end loop over variations
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Addition of independent systematics in quadrature; adds the errors of s to this systematic.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: NullifyErrors () {

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

      // Hadron yield systematics, signal & signal+bkg levels
      for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
            TH1D* master = h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master) ResetHistErrors (master);

            master = h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master) ResetHistErrors (master);

            master = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master) ResetHistErrors (master);
          } // end loop over xztrk bins
        } // end loop over cents
      } // end loop over phi

      // IAA, ICP systematics
      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];
          if (master) ResetHistErrors (master);
        } // end loop over cents
        for (short iCent = 2; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];
          if (master) ResetHistErrors (master);
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Addition of independent systematics in quadrature; adds the errors of s to this systematic.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: AddSystematic (Systematics* s) {

  for (Analysis* a : s->GetVariations ()) {
    variations.push_back (a);
  }

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

      // Hadron yield systematics, signal & signal+bkg levels
      for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
            TH1D* master = h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            TH1D* sys = s->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master && sys) AddSystematics (master, sys);

            master = h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            sys = s->h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master && sys) AddSystematics (master, sys);

            master = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            sys = s->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master && sys) AddSystematics (master, sys);
          } // end loop over xztrk bins
        } // end loop over cents
      } // end loop over phi

      // IAA, ICP systematics
      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];
          TH1D* sys = s->h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];
          if (master && sys) AddSystematics (master, sys);
        } // end loop over cents
        for (short iCent = 2; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];
          TH1D* sys = s->h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];
          if (master && sys) AddSystematics (master, sys);
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot this set of systematics
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: PlotTrkYieldSystematics (const short pSpc, const short pPtZ) {
  const char* canvasName = Form ("c_TrkYieldSys");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 600, 600);
    gDirectory->Add (c);
    gPad->SetLogx ();
  }
  c->cd ();

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made
      for (short iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {

          TH1D* source = nullptr, *errs = nullptr;
          for (short iSys = 0; iSys < variations.size (); iSys++) {
            source = variations[iSys]->h_z_trk_pt[iSpc][iPtZ][0][iPhi][iCent];
            errs = (TH1D*) source->Clone ( (string (source->GetName ()) + "_relSys").c_str ());
            SaveRelativeErrors (errs, source);

            errs->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
            errs->GetYaxis ()->SetTitle ("Relative error");

            errs->SetLineColor (colors[iSys]);
            errs->SetLineStyle (2);
            errs->SetLineWidth (2);

            errs->DrawCopy (iSys == 0 ? "hist" : "same hist");

            delete errs;
          }

          source = h_z_trk_pt[iSpc][iPtZ][0][iPhi][iCent];
          errs = (TH1D*) source->Clone ((string (source->GetName ()) + "_relSys").c_str ());
          SaveRelativeErrors (errs, source);

          errs->SetLineColor (kBlack);
          errs->SetLineStyle (1);
          errs->SetLineWidth (2);
          
          errs->DrawCopy ("same hist");
          delete errs;

          c->SaveAs (Form ("%s/TrkYieldSystematics/%s_iPtZ%i_iPhi%i_iCent%i.pdf", plotPath.Data (), spc, iPtZ, iPhi, iCent));

        } // end loop over centralities

      } // end loop over Phi bins

    } // end loop over PtZ
  } // end loop over species
  
}


#endif
