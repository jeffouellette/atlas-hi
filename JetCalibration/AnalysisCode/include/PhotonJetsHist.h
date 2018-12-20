#ifndef __PhotonJetsHist_h__
#define __PhotonJetsHist_h__

namespace JetCalibration {

const double plos[12] = {0.85, 0.7, 0.7, 0.7, 0.7, 0.7, 0.75, 0.8, 0.8, 0.75, 0.75, 0.75};
const double phis[12] = {1.3, 1.2, 1.25, 1.25, 1.3, 1.25, 1.2, 1.2, 1.2, 1.25, 1.3, 1.25};
const double plos_mc[12] = {0.85, 0.7, 0.7, 0.7, 0.6, 0.7, 0.75, 0.8, 0.8, 0.75, 0.75, 0.75};
const double phis_mc[12] = {1.3, 1.2, 1.25, 1.25, 1.3, 1.25, 1.2, 1.2, 1.2, 1.25, 1.3, 1.25};
const double etalos[7] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
const double etahis[7] = {1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};

void PhotonJetsHist ();

} // end namespace

#endif
