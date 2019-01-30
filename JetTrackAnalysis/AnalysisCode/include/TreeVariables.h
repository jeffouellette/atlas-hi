#ifndef __TreeVariables_h__
#define __TreeVariables_h__

#include <JetType.h>

#include <TTree.h>
#include <vector>

/** 
 * This file defines the TreeVariables struct used by this analysis.
 * Author: Jeff Ouellette
 * Dated: 8/20/2018
 */

using namespace std;

namespace JetTrackAnalysis {

struct TreeVariables {
  private:
    TTree* tree;

    bool isMC = false;
    bool getMCInfo = true;
    bool getCollisionRateInfo = false;
    bool getVertices = false;
    bool getFCals = false;
    bool getTracks = false;
    bool getJets = false;
    JetType jetType = FullCalibrated;
    bool getTruthJets = false;

  public:

    // public functions
    TreeVariables (TTree* t, const bool _isMC = false);
    ~TreeVariables ();    
    void SetBranchAddresses ();
    void PrintAll (const long long entry);

    // setter functions
    void SetGetMCInfo (const bool _getMCInfo = true, const double _crossSection_microbarns = 0, const double _filterEfficiency = 0, const int _numberEvents = 0);
    void SetGetCollisionRateInfo (const bool _getCollisionRateInfo = true);
    void SetGetVertices (const bool _getVertices = true);
    void SetGetFCals (const bool _getFCals = true);
    void SetGetTracks (const bool _getTracks = true);
    void SetGetJets (const bool _getJets = true, const JetType _jetType = FullCalibrated);
    void SetGetTruthJets (const bool _getTruthJets = true);

    // public variables
    int eventNumber = 0;
    int runNumber = 0;
    unsigned int lumiBlock = 0;

    int numberEvents = 0;
    double crossSection_microbarns = 0;
    double filterEfficiency = 0;

    float actualInteractionsPerCrossing = 0;
    float averageInteractionsPerCrossing = 0;

    int nvert = 0;
    vector<float>* vert_x = NULL;
    vector<float>* vert_y = NULL;
    vector<float>* vert_z = NULL;
    vector<int>* vert_ntrk = NULL;
    vector<int>* vert_type = NULL;

    // fcal energy
    float fcalA_et;
    float fcalC_et;
 
    // tracking info (0th or primary vertex only)
    int ntrk = 0;
    vector<bool>* trk_quality_4 = NULL;
    vector<float>* trk_d0 = NULL;
    vector<float>* trk_z0 = NULL;
    vector<float>* trk_theta = NULL;
    vector<float>* trk_charge = NULL;
    vector<float>* trk_pt = NULL;
    vector<float>* trk_eta = NULL;
    vector<float>* trk_phi = NULL;

    int jet_n = 0; 
    vector<float>* jet_pt = NULL;
    vector<float>* jet_eta = NULL;
    vector<float>* jet_phi = NULL;
    vector<float>* jet_e = NULL;

    int truth_jet_n = 0;
    vector<float>* truth_jet_pt = NULL;
    vector<float>* truth_jet_eta = NULL;
    vector<float>* truth_jet_phi = NULL;
    vector<float>* truth_jet_e = NULL;
   
};

} // end namespace

#endif
