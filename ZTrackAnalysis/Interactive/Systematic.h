#ifndef __Systematic_h__
#define __Systematic_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>
#include <algorithm>

using namespace std;
using namespace atlashi;

class Systematic : public PhysicsAnalysis {

  protected:
  vector<PhysicsAnalysis*> variations;
  vector<Systematic*> systematics;

  void NullifyErrors ();

  public:

  string description;

  Systematic (PhysicsAnalysis* nom, const char* _name = "systematics", const char* _desc = "systematic") : PhysicsAnalysis (){
    name = _name;
    directory = "Systematics/";
    description = _desc;
    SetupDirectories (directory, "ZTrackAnalysis/");

    CopyAnalysis (nom, true);
    NullifyErrors ();
  }

  vector<PhysicsAnalysis*>& GetVariations ()  { return variations;  }
  vector<Systematic*>&      GetSystematics () { return systematics; }

  virtual TGAE* GetTGAE (TH1* h) override;

  void AddVariation (PhysicsAnalysis* a);
  void AddSystematic (Systematic* a);

  void AddVariations (); // variations add linearly
  void AddSystematics (); // systematics add in quadrature

  void PlotTrkYieldSystematics (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  void PlotSignalTrkYieldSystematics (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  void PlotIAASystematics (const short pSpc = 2, const short pPtZ = nPtZBins-1);

};




////////////////////////////////////////////////////////////////////////////////////////////////
// Returns a TGraphAsymmErrors corresponding to this systematic
////////////////////////////////////////////////////////////////////////////////////////////////
TGAE* Systematic :: GetTGAE (TH1* h) {
  TGAE* g = make_graph (h);

  return g;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Adds a variation to consider
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: AddVariation (PhysicsAnalysis* a) {
  if (find (variations.begin (), variations.end (), a) == variations.end ())
    variations.push_back (a);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Adds a systematic error set to consider
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: AddSystematic (Systematic* a) {
  if (find (systematics.begin (), systematics.end (), a) == systematics.end ())
    systematics.push_back (a);
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Sets the errors in this systematic as a combination of all added variations.
// Takes the maximum error for each point.
// Intended for combining up & down variations, but expandable for additional categories
// (e.g. track quality criteria)
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: AddVariations () {

  NullifyErrors ();

  for (PhysicsAnalysis* a : variations) {

    cout << "Adding variation " << a->Name () << " to systematic " << name << endl;

    a->SubtractBackground ();
    a->CalculateIAA ();
    a->CalculateICP ();

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

        // Hadron yield systematics, signal & signal+bkg levels
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            TH1D* sys = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
            TH1D* var = a->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
            sys = h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
            var = a->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
            
            sys = h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
            var = a->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
            sys = h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];
            var = a->h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);

            sys = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            var = a->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
            sys = h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            var = a->h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
          } // end loop over cents
        } // end loop over phi

        // IAA, ICP systematics
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 1; iCent < numCentBins; iCent++) {
            TH1D* sys = h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];
            TH1D* var = a->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
            sys = h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent];
            var = a->h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
          } // end loop over cents
          for (short iCent = 2; iCent < numCentBins; iCent++) {
            TH1D* sys = h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent];
            TH1D* var = a->h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
            sys = h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent];
            var = a->h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent];
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
void Systematic :: NullifyErrors () {

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

      // Hadron yield systematics, signal & signal+bkg levels
      for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
          if (master) ResetHistErrors (master);
          master = h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
          if (master) ResetHistErrors (master);

          master = h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
          if (master) ResetHistErrors (master);
          master = h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];
          if (master) ResetHistErrors (master);

          master = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
          if (master) ResetHistErrors (master);
          master = h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
          if (master) ResetHistErrors (master);
        } // end loop over cents
      } // end loop over phi

      // IAA, ICP systematics
      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];
          if (master) ResetHistErrors (master);
          master = h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent];
          if (master) ResetHistErrors (master);
        } // end loop over cents
        for (short iCent = 2; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent];
          if (master) ResetHistErrors (master);
          master = h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent];
          if (master) ResetHistErrors (master);
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Addition of independent systematics in quadrature; adds the errors of s to this systematic.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: AddSystematics () {

  NullifyErrors ();

  for (PhysicsAnalysis* s : systematics) {

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

        // Hadron yield systematics, signal & signal+bkg levels
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            TH1D* master = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
            TH1D* sys = s->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
            sys = s->h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
            if (master && sys) AddErrorsInQuadrature (master, sys);

            master = h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
            sys = s->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];
            sys = s->h_z_trk_xzh_sub[iSpc][iPtZ][iPhi][iCent];
            if (master && sys) AddErrorsInQuadrature (master, sys);

            master = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            sys = s->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            sys = s->h_z_trk_xzh_sig_to_bkg[iSpc][iPtZ][iPhi][iCent];
            if (master && sys) AddErrorsInQuadrature (master, sys);
          } // end loop over cents
        } // end loop over phi

        // IAA, ICP systematics
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 1; iCent < numCentBins; iCent++) {
            TH1D* master = h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];
            TH1D* sys = s->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent];
            sys = s->h_z_trk_xzh_iaa[iSpc][iPtZ][iPhi][iCent];
            if (master && sys) AddErrorsInQuadrature (master, sys);
          } // end loop over cents
          for (short iCent = 2; iCent < numCentBins; iCent++) {
            TH1D* master = h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent];
            TH1D* sys = s->h_z_trk_pt_icp[iSpc][iPtZ][iPhi][iCent];
            if (master && sys) AddErrorsInQuadrature (master, sys);
            master = h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent];
            sys = s->h_z_trk_xzh_icp[iSpc][iPtZ][iPhi][iCent];
            if (master && sys) AddErrorsInQuadrature (master, sys);
          } // end loop over cents
        } // end loop over phi
      } // end loop over pT^Z bins
    } // end loop over species

  }

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot this set of systematics
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: PlotTrkYieldSystematics (const short pSpc, const short pPtZ) {
  const char* canvasName = Form ("c_TrkYieldSys");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
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

          TH1D* centralVals = nullptr, *errs = nullptr;

          centralVals = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];

          bool drawn = false;
          short iSys = 0;
          for (Systematic* sys : systematics) {

            if (sys->Name () == string ("bkgSys"))
              continue;

            TH1D* h = sys->h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
            errs = (TH1D*) h->Clone ( (string (h->GetName ()) + "_relSys").c_str ());
            SaveRelativeErrors (errs, centralVals);

            errs->GetXaxis ()->SetMoreLogLabels ();
            errs->GetYaxis ()->SetRangeUser (-0.5, 0.5);

            errs->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
            errs->GetYaxis ()->SetTitle ("Relative error");

            errs->SetLineColor (colors[iSys+1]);
            errs->SetLineStyle (2);
            errs->SetLineWidth (5);

            if (!drawn)
              errs->DrawCopy ("][ hist");
            else
              errs->DrawCopy ("][ hist same");
            drawn = true;

            errs->Scale (-1);

            errs->DrawCopy ("][ same hist");

            myText (0.65, 0.86-0.04*iSys, colors[iSys+1], sys->description.c_str (), 0.04);

            delete errs;
            iSys++;
          }

          errs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSys").c_str ());
          SaveRelativeErrors (errs, centralVals);

          errs->GetXaxis ()->SetMoreLogLabels ();
          errs->GetYaxis ()->SetRangeUser (-0.5, 0.5);

          errs->SetLineColor (kBlack);
          errs->SetLineStyle (1);
          errs->SetLineWidth (3);

          myText (0.65, 0.91, kBlack, "Total", 0.04);
          
          errs->DrawCopy (systematics.size () == 0 ? "][ hist" : "][ same hist");

          errs->Scale (-1);

          errs->DrawCopy ("][ same hist");

          delete errs;

          myText (0.72, 0.34, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
          if (iCent == 0)
            myText (0.72, 0.28, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
          else {
            myText (0.72, 0.28, kBlack, Form ("Pb+Pb %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
            myText (0.72, 0.22, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
          }

          c->SaveAs (Form ("%s/TrkYieldSystematics/%s_iPtZ%i_iPhi%i_iCent%i.pdf", plotPath.Data (), spc, iPtZ, iPhi, iCent));

        } // end loop over centralities

      } // end loop over Phi bins

    } // end loop over PtZ
  } // end loop over species
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot this set of systematics on the signal track yield
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: PlotSignalTrkYieldSystematics (const short pSpc, const short pPtZ) {
  const char* canvasName = Form ("c_SignalTrkYieldSys");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
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

          gPad->Clear ();

          TH1D* centralVals = nullptr, *errs = nullptr;

          centralVals = h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];

          bool drawn = false;
          short iSys = 0;
          for (Systematic* sys : systematics) {

            TH1D* h = sys->h_z_trk_pt_sub[iSpc][iPtZ][iPhi][iCent];
            errs = (TH1D*) h->Clone ( (string (h->GetName ()) + "_relSys").c_str ());
            SaveRelativeErrors (errs, centralVals);

            //if (iPhi == 1 && iCent == 3) {
            //  cout << sys->description << endl;
            //  for (int ix = 1; ix <= errs->GetNbinsX (); ix++) {
            //    cout << errs->GetBinContent (ix) << endl;
            //  }
            //}

            errs->GetXaxis ()->SetMoreLogLabels ();
            errs->GetYaxis ()->SetRangeUser (-0.4, 0.4);

            errs->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
            errs->GetYaxis ()->SetTitle ("Relative error");

            errs->SetLineColor (colors[iSys+1]);
            errs->SetLineStyle (2);
            errs->SetLineWidth (5);

            if (!drawn)
              errs->DrawCopy ("][ hist");
            else
              errs->DrawCopy ("][ hist same");
            drawn = true;

            errs->Scale (-1);

            errs->DrawCopy ("][ same hist");

            myText (0.65, 0.86-0.04*iSys, colors[iSys+1], sys->description.c_str (), 0.04);

            delete errs;
            iSys++;
          }

          errs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSys").c_str ());

            //if (iPhi == 1 && iCent == 3) {
            //  cout << "Total" << endl;
            //  for (int ix = 1; ix <= errs->GetNbinsX (); ix++) {
            //    cout << errs->GetBinError (ix) << endl;
            //  }
            //}

          SaveRelativeErrors (errs, errs);

          //if (iPhi == 1 && iCent == 3) {
          //  cout << "Total" << endl;
          //  for (int ix = 1; ix <= errs->GetNbinsX (); ix++) {
          //    cout << errs->GetBinContent (ix) << endl;
          //  }
          //}

          errs->GetXaxis ()->SetMoreLogLabels ();
          errs->GetYaxis ()->SetRangeUser (-0.4, 0.4);

          errs->SetLineColor (kBlack);
          errs->SetLineStyle (1);
          errs->SetLineWidth (3);

          myText (0.65, 0.91, kBlack, "Total", 0.04);
          
          errs->DrawCopy (systematics.size () == 0 ? "][ hist" : "][ same hist");

          errs->Scale (-1);

          errs->DrawCopy ("][ same hist");

          delete errs;

          myText (0.72, 0.34, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
          if (iCent == 0)
            myText (0.72, 0.28, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
          else {
            myText (0.72, 0.28, kBlack, Form ("Pb+Pb %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
            myText (0.72, 0.22, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
          }

          c->SaveAs (Form ("%s/TrkYieldSignalSystematics/%s_iPtZ%i_iPhi%i_iCent%i.pdf", plotPath.Data (), spc, iPtZ, iPhi, iCent));

        } // end loop over centralities

      } // end loop over Phi bins

    } // end loop over PtZ
  } // end loop over species
  
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plot this set of systematics on the signal track yield
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematic :: PlotIAASystematics (const short pSpc, const short pPtZ) {
  const char* canvasName = Form ("c_iaasys");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 800, 600);
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
        for (short iCent = 1; iCent < numCentBins; iCent++) {

          TH1D* centralVals = nullptr, *errs = nullptr;

          centralVals = h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];
          if (!centralVals) continue;

          bool drawn = false;
          short iSys = 0;
          for (Systematic* sys : systematics) {

            TH1D* h = sys->h_z_trk_pt_iaa[iSpc][iPtZ][iPhi][iCent];

            if (!h) continue;

            errs = (TH1D*) h->Clone ( (string (h->GetName ()) + "_relSys").c_str ());
            SaveRelativeErrors (errs, centralVals);

            errs->GetXaxis ()->SetMoreLogLabels ();
            errs->GetYaxis ()->SetRangeUser (-0.4, 0.4);

            errs->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
            errs->GetYaxis ()->SetTitle ("Relative error");

            errs->SetLineColor (colors[iSys+1]);
            errs->SetLineStyle (2);
            errs->SetLineWidth (5);

            if (!drawn)
              errs->DrawCopy ("][ hist");
            else
              errs->DrawCopy ("][ hist same");
            drawn = true;

            errs->Scale (-1);

            errs->DrawCopy ("][ same hist");

            myText (0.65, 0.86-0.04*iSys, colors[iSys+1], sys->description.c_str (), 0.04);

            delete errs;
            iSys++;
          }

          errs = (TH1D*) centralVals->Clone ((string (centralVals->GetName ()) + "_relSys").c_str ());
          SaveRelativeErrors (errs, centralVals);

          errs->GetXaxis ()->SetMoreLogLabels ();
          errs->GetYaxis ()->SetRangeUser (-0.4, 0.4);

          errs->SetLineColor (kBlack);
          errs->SetLineStyle (1);
          errs->SetLineWidth (3);

          myText (0.65, 0.91, kBlack, "Total", 0.04);
          
          errs->DrawCopy (systematics.size () == 0 ? "][ hist" : "][ same hist");

          errs->Scale (-1);

          errs->DrawCopy ("][ same hist");

          delete errs;

          myText (0.72, 0.34, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
          if (iCent == 0)
            myText (0.72, 0.28, kBlack, "#it{pp}, #sqrt{s} = 5.02 TeV", 0.04);
          else {
            myText (0.72, 0.28, kBlack, Form ("Pb+Pb %i-%i%%", centCuts[iCent], centCuts[iCent-1]), 0.04);
            myText (0.72, 0.22, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
          }

          c->SaveAs (Form ("%s/IAASystematics/%s_iPtZ%i_iPhi%i_iCent%i.pdf", plotPath.Data (), spc, iPtZ, iPhi, iCent));

        } // end loop over centralities

      } // end loop over Phi bins

    } // end loop over PtZ
  } // end loop over species
  
}

#endif
