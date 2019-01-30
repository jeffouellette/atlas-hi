#ifndef __Track_h__
#define __Track_h__

struct Track {
  public:
    float Charge = 0;
    float Pt = 0;
    float Eta = 0;
    float Phi = 0;

    Track () {}
    Track (float _Pt, float _Eta, float _Phi, float _Charge) {
      Pt = _Pt;
      Eta = _Eta;
      Phi = _Phi;
      Charge = _Charge;
    }

    ~Track () {}
};

#endif
