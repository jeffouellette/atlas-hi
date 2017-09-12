#ifndef jetAnalysis_jetTreeMaker_H
#define jetAnalysis_jetTreeMaker_H

#include <EventLoop/Algorithm.h>

// root includes
#include <TTree.h>

#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

#include <vector>

// initial tree array lengths
static const int max_njet = 60; //!
static const int max_nvertex = 100; //!

class jetTreeMaker : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;


  float m_jet_pt_cut;
  std::string outputName;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!



  // this is a standard constructor
  jetTreeMaker ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(jetTreeMaker, 1);


private:

  TTree* m_tree; //!

  // Tree branches 
  int m_b_runNumber; //!
  int m_b_eventNumber; //!
  unsigned int m_b_lumiBlock; //! 
  float m_b_actualInteractionsPerCrossing; //!
  float m_b_averageInteractionsPerCrossing; //!

  // global event info
  int m_b_ntrk; //!
  float m_b_fcalA_et; //!
  float m_b_fcalC_et; //!
  int m_b_nvert; //!
  float m_b_vert_x[max_nvertex]; //!
  float m_b_vert_y[max_nvertex]; //!
  float m_b_vert_z[max_nvertex]; //!
  int m_b_vert_ntrk[max_nvertex]; //!
  float m_b_sumpT[max_nvertex]; //!

  // triggers
  //
  static const int trigLength = 24; //!
  const char* m_trig_string[trigLength] = 
  {
      "HLT_j30_0eta490_L1TE10",
      "HLT_j30_ion_0eta490_L1TE10",
      "HLT_j40_L1J5",
      "HLT_j40_ion_L1J5",
      "HLT_j50_L1J10",
      "HLT_j50_ion_L1J10",
      "HLT_j60",
      "HLT_j60_ion_L1J20",
      "HLT_j75_L1J20",
      "HLT_j75_ion_L1J20",
      "HLT_j100_L1J20",
      "HLT_j100_ion_L1J20",
      "HLT_j15_p320eta490_L1MBTS_1_1",
      "HLT_j15_ion_p320eta490_L1MBTS_1_1",
      "HLT_j55_p320eta490",
      "HLT_j55_ion_p320eta490",
      "HLT_j45_p200eta320",
      "HLT_j45_ion_p200eta320",
      "HLT_j65_p200eta320",
      "HLT_j65_ion_p200eta320",
      "HLT_2j10_p320eta490_L1TE10",
      "HLT_2j10_ion_p320eta490_L1TE10",
      "HLT_2j30_p320eta490",
      "HLT_2j30_ion_p320eta490"
  };
  bool m_trig_bool[trigLength]; //!
  float m_trig_prescale[trigLength]; //!

  // jet info
  int m_b_njet; //!
  float m_b_j_pt[max_njet]; //!
  float m_b_j_eta[max_njet]; //!
  float m_b_j_phi[max_njet]; //!
  float m_b_j_e[max_njet]; //!

  GoodRunsListSelectionTool* m_grl; //!
  Trig::TrigDecisionTool* m_trigDecisionTool; //!
  TrigConf::xAODConfigTool* m_trigConfigTool; //!


};

#endif
