#ifndef __MixingVariation_h__
#define __MixingVariation_h__

#include "Params.h"
#include "PhysicsAnalysis.h"

////////////////////////////////////////////////////////////////////////////////////////////////
// Applies a blanket relative variation on the background hadron yield.
////////////////////////////////////////////////////////////////////////////////////////////////
void ApplyBkgVariation (PhysicsAnalysis* a, const float relErr) {

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) { 

      for (int iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          a->h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Scale (1+relErr);
        } // end loop over cents
      } // end loop over iPtTrk

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        a->h_trk_pt_ptz[iSpc][iPtZ][iCent]->Scale (1+relErr);
        a->h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Scale (1+relErr);
      } // end loop over cents
    } // end loop over pT^Z bins
  } // end loop over species

}


#endif
