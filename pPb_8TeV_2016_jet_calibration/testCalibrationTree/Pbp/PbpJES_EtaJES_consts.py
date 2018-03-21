## Jet energy calibration
##   [0]+[1]*log(x)+[2]*log(x)^2+...
energyCorrDict[ 'AntiKt4HIJets_HI' ] = [
  [   9.7478e-01,  4.1136e-03,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  0 -0.9<eta<-0.8
  [   8.3368e-01,  8.5786e-02, -1.4770e-02,  8.4738e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  1 -0.8<eta<-0.7
  [   1.0091e+00, -1.8037e-02,  4.8401e-03, -3.4642e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  2 -0.7<eta<-0.6
  [   9.8731e-01,  1.8280e-03,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  3 -0.6<eta<-0.5
  [   9.7262e-01,  7.6290e-03, -5.4118e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  4 -0.5<eta<-0.4
  [   8.6563e-01,  6.4004e-02, -1.0381e-02,  5.6478e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  5 -0.4<eta<-0.3
  [   9.9285e-01, -1.0212e-02,  3.9256e-03, -3.3690e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  6 -0.3<eta<-0.2
  [   1.6823e+00, -5.7702e-01,  1.7467e-01, -2.2732e-02,  1.0824e-03,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  7 -0.2<eta<-0.1
  [   8.6248e-01,  6.8041e-02, -1.1393e-02,  6.3622e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  8 -0.1<eta<-0.0
  [   1.0527e+00, -4.2936e-02,  9.5919e-03, -6.5291e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  9  0.0<eta< 0.1
  [   9.3891e-01,  1.9885e-02, -1.6518e-03,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 10  0.1<eta< 0.2
  [   9.6049e-01,  1.1844e-02, -9.4939e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 11  0.2<eta< 0.3
  [   1.2508e+00, -2.3910e-01,  7.9279e-02, -1.1146e-02,  5.6744e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 12  0.3<eta< 0.4
  [   9.5199e-01,  1.5738e-02, -1.3486e-03,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 13  0.4<eta< 0.5
  [   1.6200e+00, -5.4755e-01,  1.7199e-01, -2.3101e-02,  1.1276e-03,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 14  0.5<eta< 0.6
  [   7.6202e-01,  1.2532e-01, -2.2071e-02,  1.2841e-03,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 15  0.6<eta< 0.7
  [   1.1026e+00, -1.2959e-01,  4.9355e-02, -7.5694e-03,  4.0792e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 16  0.7<eta< 0.8
  [   1.4090e+00, -3.4820e-01,  1.0435e-01, -1.3409e-02,  6.2864e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ]  ## bin 17  0.8<eta< 0.9
]

## Jet eta correction
##   [0]+[1]*log(x)+[2]*log(x)^2+...
etaCorrDict[ 'AntiKt4HIJets_HI' ] = [
  [   2.5152e-02, -1.0139e-02,  1.4738e-03, -7.2060e-05,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  0 -0.9<eta<-0.8
  [  -4.4556e-02,  2.4255e-02, -4.0276e-03,  2.1447e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  1 -0.8<eta<-0.7
  [  -1.8911e-02,  1.1657e-02, -2.0290e-03,  1.1259e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  2 -0.7<eta<-0.6
  [   3.8588e-03, -2.5626e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  3 -0.6<eta<-0.5
  [  -5.4604e-02,  3.1655e-02, -5.6372e-03,  3.2492e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  4 -0.5<eta<-0.4
  [  -3.9972e-02,  2.3945e-02, -4.2774e-03,  2.4529e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  5 -0.4<eta<-0.3
  [  -3.0714e-02,  1.8471e-02, -3.2767e-03,  1.8825e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  6 -0.3<eta<-0.2
  [  -1.8015e-02,  1.1043e-02, -1.8674e-03,  1.0302e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  7 -0.2<eta<-0.1
  [  -5.4530e-04,  1.4877e-03, -1.4518e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  8 -0.1<eta<-0.0
  [   4.3807e-03, -2.4278e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  9  0.0<eta< 0.1
  [   4.5342e-03, -2.2159e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 10  0.1<eta< 0.2
  [   3.9666e-02, -1.9997e-02,  3.5291e-03, -2.0325e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 11  0.2<eta< 0.3
  [   3.8777e-02, -1.8842e-02,  3.2855e-03, -1.9017e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 12  0.3<eta< 0.4
  [   6.6699e-02, -3.5072e-02,  6.2253e-03, -3.5969e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 13  0.4<eta< 0.5
  [   2.8421e-03,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 14  0.5<eta< 0.6
  [   2.5476e-03,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 15  0.6<eta< 0.7
  [   5.3236e-02, -2.6352e-02,  4.4487e-03, -2.4585e-04,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 16  0.7<eta< 0.8
  [  -1.4869e-02,  7.6368e-03, -1.1122e-03,  5.3265e-05,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ]  ## bin 17  0.8<eta< 0.9
]

## Jet pT correction (input EtaJESpT)
##   [0]+[1]*log(x)+[2]*log(x)^2+...
pTCorrDict[ 'AntiKt4HIJets_HI' ] = [
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  0 -0.9<eta<-0.8
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  1 -0.8<eta<-0.7
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  2 -0.7<eta<-0.6
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  3 -0.6<eta<-0.5
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  4 -0.5<eta<-0.4
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  5 -0.4<eta<-0.3
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  6 -0.3<eta<-0.2
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  7 -0.2<eta<-0.1
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  8 -0.1<eta<-0.0
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin  9  0.0<eta< 0.1
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 10  0.1<eta< 0.2
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 11  0.2<eta< 0.3
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 12  0.3<eta< 0.4
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 13  0.4<eta< 0.5
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 14  0.5<eta< 0.6
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 15  0.6<eta< 0.7
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ], ## bin 16  0.7<eta< 0.8
  [   0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00 ]  ## bin 17  0.8<eta< 0.9
]

