// covarience sensitive fourier fit
//#include "upcUtils.C"
// author: Blair Seidlitz, Jeff Ouellette

#include <TH1D.h>
#include <TH2D.h>
#include <Fit/Fitter.h>
#include <TMatrixD.h>

class FitCov {

  private:

  TH1D* h = nullptr;
  TH2D* hcov = nullptr;

  static double xmin;
  static double xmax;
  static int nBins;
  static double* z;
  static double* x;
  static double** covInv;
  static int verbosity;
  //static double* errorZ = nullptr;

  public:

  ROOT::Fit::FitResult result;

  bool noCov = 0;

  static double GetXMin () { return xmin; }
  static double GetXMax () { return xmax; }
  static void SetXMin (double _xmin) { xmin = _xmin; return; }
  static void SetXMax (double _xmax) { xmax = _xmax; return; }

  FitCov (TH1D* _h, TH2D* _hcov, int nPar, double* par, double* parCov);

  static void fcn (int &nPar, double* gin, double &f, double* par, int iFlag);
  static double PowerLaw (double xval, double *par);
};


double FitCov::xmin = DBL_MIN;
double FitCov::xmax = DBL_MAX;

int FitCov::nBins = -1;
double* FitCov::z = nullptr;
double* FitCov::x = nullptr;
double** FitCov::covInv = nullptr;
int FitCov::verbosity = 0;




/**
 * Calculates chisquare at given set of parameters.
 */
void FitCov :: fcn (int &nPar, double* gin, double &f, double* par, int iFlag) {
  double chisq = 0;
  for (int i = 0; i < nBins; i++) {
    if (x[i] < xmin || xmax < x[i]) continue;
    for (int j = 0; j < nBins; j++) {
      if (x[j] < xmin || xmax < x[j]) continue;
      chisq += (z[i]-PowerLaw (x[i], par)) * (z[j]-PowerLaw (x[j], par)) * (covInv[i][j]);
    }
  }
  f = chisq;
}




FitCov :: FitCov (TH1D* _h, TH2D* _hcov, int nPar, double* par, double* parCov) {

  assert (_h != nullptr && _hcov != nullptr);
  assert (_h->GetNbinsX () == _hcov->GetNbinsX ());

  h = _h;
  hcov = _hcov;

  if (nBins != -1) {
    if (z) delete[] z;
    if (x) delete[] x;
    if (covInv) {
      for (int i = 0; i < nBins; i++)
        delete[] covInv[i];
      delete[] covInv;
    }
  }

  nBins = h->GetNbinsX ();

  x = new double[nBins];
  z = new double[nBins];
  covInv = new double*[nBins];
  for (int i = 0; i < nBins; i++)
    covInv[i] = new double[nBins];

  double cov[nBins*nBins];
  for (int i = 0; i < nBins; i++) {
    for (int j = 0; j < nBins; j++) {
      cov[i*nBins+j] = hcov->GetBinContent(j+1, i+1);
      //cout << cov[i*j+j] << " ";
      if (i != j && noCov) cov[i*nBins+j] = 0;
      if (i == j)          cov[i*nBins+j] = 1;
    }
    //cout << endl;
  }

  TMatrixD covMat = TMatrixD (nBins, nBins, cov);
  //for (int i = 0; i < nBins; i++) {
  //    for (int j = 0; j < nBins; j++) {
  //      cout << (float) covMat[i][j] << " ";
  //    }
  //    cout << endl;
  //}
  covMat.Invert();
  for (int i = 0; i < nBins; i++) {
    for (int j = 0; j < nBins; j++) {
      covInv[i][j] = covMat[i][j];
    }
  }

  for (int i = 0; i < nBins; i++) {
    x[i] = (double) h->GetBinCenter (i+1);
    z[i] = (double) h->GetBinContent (i+1);
    //errorZ[i] = (double) h->GetBinError (i+1);
  }

  ROOT::Fit::Fitter fitter;
  fitter.Config ().SetParamsSettings (nPar, par);
  fitter.Config ().MinimizerOptions ().SetPrintLevel (verbosity);
  fitter.Config ().SetMinimizer ("Minuit2", "Migrad");

  for (int iPar = 0; iPar < nPar; iPar++) {
    fitter.Config ().ParSettings (iPar).SetName (Form ("par%i", iPar));
    par[iPar] = 0;
  }

  double chi2;

  fitter.SetFCN (&fcn, nPar, par, nBins, true);
  fitter.FitFCN (&fcn, nPar, par, nBins, true);
  result = fitter.Result ();

  // save results
  for (int iPar = 0; iPar < nPar; iPar++){
    par[iPar] = result.Value (iPar);
    for (int jPar = 0; jPar < nPar; jPar++) {
      parCov[iPar*nPar+jPar] = result.CovMatrix (iPar, jPar);
    }
  }

}




double FitCov :: PowerLaw (double xval, double* par) {
  double fitval = par[0] * pow (xval, par[1]);
  return fitval;
}

