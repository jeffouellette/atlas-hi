#ifndef __Jet_h__
#define __Jet_h__

#include "JetType.h"

struct Jet {
  public:
    JetType jetType = FullCalibrated;
    float Pt = 0;
    float Eta = 0;
    float Phi = 0;
    float E = 0;
    float EvtWeight = 0;

    Jet () {}
    Jet (const float _Pt, const float _Eta, const float _Phi, const float _E, const float _EvtWeight, const JetType _jetType = FullCalibrated) {
      Pt = _Pt;
      Eta = _Eta;
      Phi = _Phi;
      E = _E;
      EvtWeight = _EvtWeight;
      jetType = _jetType;
    }
    Jet (const Jet* jet) {
      Pt = jet->Pt;
      Eta = jet->Eta;
      Phi = jet->Phi;
      E = jet->E;
      EvtWeight = jet->EvtWeight;
      jetType = jet->jetType;
    }

    ~Jet () {}
};

#endif
