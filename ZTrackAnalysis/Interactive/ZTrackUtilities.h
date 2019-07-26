#ifndef __ZTrackUtilities_h__
#define __ZTrackUtilities_h__

#include "Params.h"

#include <TBox.h>
#include <TPad.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>

short GetidPhi (const float dphi) {
  short idPhi = 0;
  while (idPhi < numPhiBins) {
    if (phiLowBins[idPhi] < dphi && dphi < phiHighBins[idPhi])
      break;
    else
      idPhi++;
  }
  return idPhi;
}


short GetiZH (const float zH) {
  short iZH = 0;
  while (iZH < nZHBins) {
    if (zHBins[iZH+1] < zH)
      iZH++;
    else
      break;
  }
  return iZH;
}


TBox* TBoxNDC (const double x1, const double y1, const double x2, const double y2) {
  TPad* p;
  if (gDirectory->Get ("box_pad"))
    p = (TPad*)gDirectory->Get ("box_pad");
  else {
    p = new TPad ("box_pad", "box_pad", 0., 0., 1., 1.);
    p->SetFillStyle (0);
  } 
  p->Draw ();
  p->cd ();
  TBox* b = new TBox (x1, y1, x2, y2);
  return b;
}


void RescaleWithError (TH1* h, double sf, double err) {
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    double y = h->GetBinContent (ix), yerr = h->GetBinError (ix);
    h->SetBinContent (ix, y / sf);
    h->SetBinError (ix, sqrt (pow (yerr / sf, 2) + pow (y * err / (sf * sf), 2)));
  }
  return;
}

#endif
