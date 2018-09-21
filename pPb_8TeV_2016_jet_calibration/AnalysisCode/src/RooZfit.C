#include "RooZfit.h"

namespace pPb8TeV2016JetCalibration {

RooZfit::RooZfit (TH1F* hist, Color_t color = kBlack, bool verbose = true) {
  mass = new RooRealVar ("mass","m_{ee}", 60, 110, "GeV");
  mean = new RooRealVar ("mean","Zee Mass mean", Z_mass, 85, 95);
  sigma = new RooRealVar ("sigma","Zee reco resolution", 0.5*Z_width, 0, 10);
  width = new RooRealVar ("width","Zee decay width", Z_width, 0, 20);

  coeff = new RooRealVar ("exp","exp", 0);
  sigFrac = new RooRealVar ("sigFrac","signal fraction", 1);

  model_sig = new RooVoigtian ("model_sig", "sig pdf", *mass, *mean, *width, *sigma);
  model_bkg = new RooExponential ("model_bkg", "bkg pdf", *mass, *coeff);

  model = new RooAddPdf ("model", "model", RooArgSet (*model_sig, *model_bkg), *sigFrac);

  hist_data = new RooDataHist ("hist_data", "binned data", *mass, Import (*hist));

  RooFitResult* fr = model->fitTo (*hist_data, Save (kTRUE));
  if (verbose) fr->Print ("v"); 

  plot = mass->frame (Bins (50), Title ("Zee mass"));
  hist_data->plotOn (plot,Name ("m_data"), MarkerStyle (24));
  model->plotOn (plot,Name ("m_model"),LineColor (color));
  reducedChi2 = plot->chiSquare ("m_model","m_data",0);

//  RooHist* hpull = m_mframe->pullHist ();
//  RooPlot* m_pullframe = mass.frame (Bins (50), Title ("Zee mass fit pull"));
//  m_pullframe->addPlotable (hpull,"P") ;

}

RooZfit::~RooZfit () {
  if (mass) delete mass;
  if (mean) delete mean;
  if (sigma) delete sigma;
  if (width) delete width;
  if (coeff) delete coeff;
  if (sigFrac) delete sigFrac;
  if (model_sig) delete model_sig;
  if (model_bkg) delete model_bkg;
  if (model) delete model;
  if (hist_data) delete hist_data;
}

} // end namespace
