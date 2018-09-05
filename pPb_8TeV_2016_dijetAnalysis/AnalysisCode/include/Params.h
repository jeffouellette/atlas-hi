#ifndef __Params_h__
#define __Params_h__

#include <vector>

#include <TString.h>

using namespace std;

namespace pPb8TeV2016DijetAnalysis {

extern const int useDataVersion;

extern const bool runPeriodA;
extern const bool runPeriodB;
extern const bool runData;
extern const bool runMC;

extern const bool fillTracksTree;

extern const bool scaleBackHEC;
extern const bool scaleAnalyses;
extern const TString fittedFunctionType;

extern const double dijetPtRatioCut;
extern const double dijetMinimumPt;
extern const double thirdJetMaximumPt;
extern const bool highPtJetsOnly;

extern const vector<int> runNumbers;
extern const int& numruns;
extern const vector<TString> mcSamples;
extern const int& nummcs;


extern const double pbins[68];
extern const int numpbins;
extern const double etabins[9];
//extern const double etabins[2];
extern const int numetabins;

} // end namespace

#endif

