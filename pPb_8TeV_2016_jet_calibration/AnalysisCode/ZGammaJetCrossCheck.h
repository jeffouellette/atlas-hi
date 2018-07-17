#include "../Params.C"
#include "../../Initialization.C"
#include "../TreeVariables.C"

//int eventNumber = 0;
//double eventWeight = 0;
//
//int nvert = 0;
//vector<float>* vert_z = NULL;
//vector<int>* vert_ntrk = NULL;
//vector<int>* vert_type = NULL;
//
//int jet_n = 0; 
//vector<float>* jet_pt = NULL;
//vector<float>* jet_eta = NULL;
//vector<float>* jet_phi = NULL;
//vector<float>* jet_e = NULL;
//vector<float>* init_jet_pt = NULL;
//
//int muon_n = 0;
//vector<float>* muon_pt = NULL;
//vector<float>* muon_eta = NULL;
//vector<float>* muon_phi = NULL;
//vector<int>* muon_quality = NULL;
//vector<int>* muon_charge = NULL;
//vector<bool>* muon_tight = NULL;
//vector<bool>* muon_loose = NULL;
//
//int electron_n = 0;
//vector<float>* electron_pt = NULL;
//vector<float>* electron_eta = NULL;
//vector<float>* electron_phi = NULL;
//vector<int>* electron_charge = NULL;
//vector<bool>* electron_loose = NULL;
//vector<bool>* electron_tight = NULL;
//
//int photon_n = 0;
//vector<float>* photon_pt = NULL;
//vector<float>* photon_eta = NULL;
//vector<float>* photon_phi = NULL;
//vector<bool>* photon_tight = NULL;
//vector<unsigned int>* photon_isem = NULL;
//vector<int>* photon_convFlag = NULL;
//vector<float>* photon_Rconv = NULL;
//vector<float>* photon_topoetcone40 = NULL;

static const short electronTrigLength = 3;
static const char* electronTriggerNames[electronTrigLength] = {
 "HLT_e20_lhloose",
 "HLT_e22_lhloose",
 "HLT_e24_lhloose"
};
static const float electronTriggerMinPtCuts[electronTrigLength] = {20, 22, 24};
static const float electronTriggerMaxPtCuts[electronTrigLength] = {8000, 8000, 8000};

static const short muonTrigLength = 3;
static const char* muonTriggerNames[muonTrigLength] = {
 "HLT_mu15",
 "HLT_mu18",
 //"HLT_mu20",
 "HLT_mu20_L1MU15"
};
static const float muonTriggerMinPtCuts[muonTrigLength] = {15, 18, 20};
static const float muonTriggerMaxPtCuts[muonTrigLength] = {18, 20, 8000};

static const short photonTrigLength = 7;
static const char* photonTriggerNames[photonTrigLength] = {
 "HLT_g10_loose",
 "HLT_g15_loose",
 "HLT_g20_loose",
 "HLT_g25_loose",
 "HLT_g30_loose",
 "HLT_g35_loose",
 "HLT_g60_loose"
};
static const float photonTriggerMinPtCuts[photonTrigLength] = {15, 20, 25, 30, 35, 40, 65};
static const float photonTriggerMaxPtCuts[photonTrigLength] = {20, 25, 30, 35, 40, 65, 8000};

//vector<Trigger*> electronTriggers = {};
//vector<Trigger*> muonTriggers = {};
//vector<Trigger*> photonTriggers = {};

/**
 * Calculates the delta R between two 4-vectors given their eta, phi coordinates.
 * eta1: eta of 1st vector
 * eta2: eta of 2nd vector
 * phi1: phi of 1st vector
 * phi2: phi of 2nd vector
 */
double deltaR (const double eta1, const double eta2, const double phi1, const double phi2 );

/**
 * Calculates the original systematic error on this jet from the cross calib.
 * jpt: pt of the jet
 * jeta: eta of the jet
 */
double GetXCalibSystematicError(const double jpt, const double jeta);

/**
 * Calculates the additional systematic on jet pt given the jet eta and the
 * reference vector boson pt^ref.
 * jeta: eta of the jet
 * refpt: reference pt of the vector boson
 */
double GetNewXCalibSystematicError(TFile* file, const double jeta, const double refpt);

/**
 * Primary macro.
 * dataSet: Data set identifier. This should be a run number for data or some other identifier for MC (e.g., slice number).
 * luminosity: Integrated luminosity for this run. Presumed constant over the run period. Meaningless for MC.
 * isMC: is data/MC flag.
 * isMCperiodAflag: flag that is raised for MC (meaningless if isMC is false)
 * inFileName: Input root file name where tree is stored; if == "" code will try to guess file name based on other info
 */
void ZGammaJetCrossCheck (const int dataSet,
                          const double luminosity = 0, 
                          const bool isMC = false,
                          const bool isMCperiodAflag = false, 
                          const string inFileName = "");

