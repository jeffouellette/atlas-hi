#ifndef __Signal_h__
#define __Signal_h__

#include "Params.h"
#include "Analysis.h"
#include "ZTrackUtilities.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class Signal : public Analysis {

  private:
  Analysis* obs;
  Analysis* bkg;

  public:
  //TH1D***** h_z_trk_pt_sub = nullptr;
  //TH1D***** h_z_trk_pt_iaa = nullptr;
  //TH1D***** h_z_trk_pt_icp = nullptr;

  Signal (Analysis* _obs, Analysis* _bkg);

  void SubtractBackground () override;

  //void CalculateIAA ();
  //void CalculateICP ();

  //void PlotTrkYield (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  //void PlotIAARatios (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  //void PlotICPRatios (const short pSpc = 2, const short pPtZ = nPtZBins-1);
};


////////////////////////////////////////////////////////////////////////////////////////////////
// Construct a signal object
////////////////////////////////////////////////////////////////////////////////////////////////
Signal::Signal (Analysis* _obs, Analysis* _bkg) : Analysis (_obs) {

  name = "signal";
  directory = "Signal/";
  plotFill = false;

  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");
  assert (!_obs);
  assert (!_bkg);

  obs = _obs;
  bkg = _bkg;

  h_z_trk_pt_sub = Get5DArray <TH1D*> (3, nPtZBins, nXZTrkBins, numPhiBins, numCentBins); // iSpc, iPtZ, iXZTrk, iPhi, iCent
  h_z_trk_pt_sig_to_bkg  = Get5DArray <TH1D*> (3, nPtZBins, nXZTrkBins, numPhiBins, numCentBins); // iSpc, iPtZ, iXZTrk, iPhi, iCent
  h_z_trk_pt_iaa = Get5DArray <TH1D*> (3, nPtZBins, nXZTrkBins, numPhiBins, numCentBins); // iSpc, iPtZ, iXZTrk, iPhi, iCent
  h_z_trk_pt_icp = Get5DArray <TH1D*> (3, nPtZBins, nXZTrkBins, numPhiBins, numCentBins); // iSpc, iPtZ, iXZTrk, iPhi, iCent

  backgroundSubtracted = false;
  SubtractBackground ();
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Subtracts a background from the observed track yields
////////////////////////////////////////////////////////////////////////////////////////////////
void Signal::SubtractBackground () {
  if (backgroundSubtracted)
    return;
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            TH1D* h = new TH1D (Form ("h_z_trk_pt_sub_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_signal_%s_%s", spc, iPtZ, iXZTrk, iPhi, iCent, obs->Name (), bkg->Name ()), "", nPtTrkBins, ptTrkBins);
            h->Sumw2 ();

            h->Add (obs->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]);
            h->Add (bkg->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent], -1);

            h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent] = h;

            h = new TH1D (Form ("h_z_trk_pt_sigToBkg_%s_iPtZ%i_iXZTrk%i_iPhi%i_iCent%i_%s", spc, iPtZ, iXZTrk, iPhi, iCent, name), "", nPtTrkBins, ptTrkBins);
            h->Sumw2 ();

            h->Add (obs->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]);
            h->Add (bkg->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent], -1);
            h->Divide (bkg->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]);

            h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent] = h;
          } // end loop over phi
        } // end loop over centralities
      } // end loop over xZ^Trk
    } // end loop over pT^Z bins
  } // end loop over species
  backgroundSubtracted = true;
}

#endif
