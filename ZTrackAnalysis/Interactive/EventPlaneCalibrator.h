#ifndef __EventPlaneCalibrator_h__
#define __EventPlaneCalibrator_h__

#include <TFile.h>
#include <TH1D.h>

#include <string>

struct EventPlaneCalibrator {

  private:
  // Event plane calibration requires private information

  // Correlation functions of qx and qy
  float _qx_a = 0;
  float _qy_a = 0;
  float _qx_c = 0;
  float _qy_c = 0;
  float _qxqxhat_a = 0;
  float _qxqyhat_a = 0;
  float _qyqyhat_a = 0;
  float _qxqxhat_c = 0;
  float _qxqyhat_c = 0;
  float _qyqyhat_c = 0;

  // Derived quantities
  float _D_a = 0;
  float _D_c = 0;
  float _N_a = 0;
  float _N_c = 0;

  public:

  EventPlaneCalibrator () {}
  EventPlaneCalibrator (const char* _fname) {
    LoadCalibration (_fname);
  }
  ~EventPlaneCalibrator () {}

  void LoadCalibration (const char* file);

  float CalibrateQ2XA (float qx_a, float qy_a);
  float CalibrateQ2YA (float qx_a, float qy_a);
  float CalibrateQ2XC (float qx_c, float qy_c);
  float CalibrateQ2YC (float qx_c, float qy_c);
};


////////////////////////////////////////////////////////////////////////////////////////////////
// Loads calibration matrix from a specified file
////////////////////////////////////////////////////////////////////////////////////////////////
void EventPlaneCalibrator :: LoadCalibration (const char* file) {
  TFile* f = new TFile (file, "read");

  _qx_a = ((TH1D*) f->Get ("h_qx_a"))->GetMean ();
  _qy_a = ((TH1D*) f->Get ("h_qx_a"))->GetMean ();
  _qx_a = ((TH1D*) f->Get ("h_qx_a"))->GetMean ();
  _qy_a = ((TH1D*) f->Get ("h_qx_a"))->GetMean ();
  const float _qxqx_a = ((TH1D*) f->Get ("h_qxqx_a"))->GetMean ();
  const float _qxqy_a = ((TH1D*) f->Get ("h_qxqy_a"))->GetMean ();
  const float _qyqy_a = ((TH1D*) f->Get ("h_qyqy_a"))->GetMean ();
  const float _qxqx_c = ((TH1D*) f->Get ("h_qxqx_c"))->GetMean ();
  const float _qxqy_c = ((TH1D*) f->Get ("h_qxqy_c"))->GetMean ();
  const float _qyqy_c = ((TH1D*) f->Get ("h_qyqy_c"))->GetMean ();

  f->Close ();
  SaferDelete (f);

  _qxqxhat_a = _qxqx_a - _qx_a*_qx_a;
  _qxqyhat_a = _qxqy_a - _qx_a*_qy_a;
  _qyqyhat_a = _qyqy_a - _qy_a*_qy_a;

  _qxqxhat_c = _qxqx_c - _qx_c*_qx_c;
  _qxqyhat_c = _qxqy_c - _qx_c*_qy_c;
  _qyqyhat_c = _qyqy_c - _qy_c*_qy_c;

  _D_a = _qxqxhat_a * _qyqyhat_a - _qxqyhat_a * _qxqyhat_a;
  _D_c = _qxqxhat_c * _qyqyhat_c - _qxqyhat_c * _qxqyhat_c;

  _N_a = _D_a * (_qxqxhat_a + _qyqyhat_a + 2*_D_a);
  _N_c = _D_c * (_qxqxhat_c + _qyqyhat_c + 2*_D_c);
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Calibrates the qx value in the A side
////////////////////////////////////////////////////////////////////////////////////////////////
float EventPlaneCalibrator :: CalibrateQ2XA (float qx_a, float qy_a) {
  qx_a = qx_a - _qx_a;
  qy_a = qy_a - _qy_a;
  return (1./sqrt (_N_a)) * ((_qyqyhat_a+_D_a)*(qx_a) + (-_qxqyhat_a)*(qy_a));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Calibrates the qy value in the A side
////////////////////////////////////////////////////////////////////////////////////////////////
float EventPlaneCalibrator :: CalibrateQ2YA (float qx_a, float qy_a) {
  qx_a = qx_a - _qx_a;
  qy_a = qy_a - _qy_a;
  return (1./sqrt (_N_a)) * ((-_qxqyhat_a)*(qx_a) + (_qxqxhat_a+_D_a)*(qy_a));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Calibrates the qx value in the C side
////////////////////////////////////////////////////////////////////////////////////////////////
float EventPlaneCalibrator :: CalibrateQ2XC (float qx_c, float qy_c) {
  qx_c = qx_c - _qx_c;
  qy_c = qy_c - _qy_c;
  return (1./sqrt (_N_c)) * ((_qyqyhat_c+_D_c)*(qx_c) + (-_qxqyhat_c)*(qy_c));
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Calibrates the qy value in the C side
////////////////////////////////////////////////////////////////////////////////////////////////
float EventPlaneCalibrator :: CalibrateQ2YC (float qx_c, float qy_c) {
  qx_c = qx_c - _qx_c;
  qy_c = qy_c - _qy_c;
  return (1./sqrt (_N_c)) * ((-_qxqyhat_c)*(qx_c) + (_qxqxhat_c+_D_c)*(qy_c));
}

#endif
