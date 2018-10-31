#ifndef __RooZfit_h__
#define __RooZfit_h__

#include <GlobalParams.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooVoigtian.h>
#include <RooExponential.h>
#include <RooConstVar.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>

using namespace atlashi;
using namespace RooFit;

namespace JetCalibration {

struct RooZfit {
  public:
   RooRealVar* mass;
   RooRealVar* mean;
   RooRealVar* sigma;
   RooRealVar* width;
   RooRealVar* coeff;
   RooRealVar* sigFrac;
 
   RooVoigtian* model_sig;
   RooExponential* model_bkg;

   RooDataHist* hist_data;
 
   RooAddPdf* model;
   double reducedChi2;
 
   RooPlot* plot;

   RooZfit (TH1F* hist, Color_t color, bool verbose);
   ~RooZfit ();
};

} // end namespace

#endif
