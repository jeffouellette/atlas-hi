#ifndef __Util_h__
#define __Util_h__

#include "Params.h"

#include <Initialization.h>
#include <Trigger.h>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>

#include <vector>

using namespace atlashi;

namespace pPb8TeV2016DijetAnalysis {

//=========================================================================================================
/** Variable declarations **/

// Trigger vectors used for pt and eta binning
double** kinematicLumiVec; // total exposed luminosity in a kinematic bin. Units: ub^-1
Trigger*** kinematicTriggerVec; // best trigger to use in a kinematic bin.


//=========================================================================================================
// General functions


/**
 * Returns xp for the event.
 */
double get_xp(double jpt0, double jpt1, double jeta0, double jeta1, bool periodA);


/**
 * Returns xa for the event.
 */
double get_xa(double jpt0, double jpt1, double jeta0, double jeta1, bool periodA);


/**
 * Returns the momentum transfer ("hardness") Q for the event.
 */
double get_q2(double xp, double je, double jpt);


/**
 * Returns the dijet mass Mjj.
 */
double get_mjj(TLorentzVector jet0, TLorentzVector jet1);


/**
 * Returns the bin number of pt in pbins.
 * Returns -1 iff pt < pbins[0] or (numpbins+1) iff pbins[numpbins] <= pt.
 * Otherwise returns n iff pbins[n] <= jpt and jpt < pbins[n+1], etc.
 */
int getPbin (const double pt);


/**
 * Returns the bin number of eta in etabins.
 * Returns -1 if eta < etabins[0] and numetabins if etabins[numetabins] <= eta.
 * Otherwise returns n iff etabins[n] <= eta && eta < etabins[n+1], etc.
 */
int getEtabin (const double eta);


/**
 * Determines whether analysis of this run should be skipped.
 */
bool SkipRun (const bool pA);


/**
 * Determines whether analysis of this MC sample should be skipped.
 */
bool SkipMC (const bool pA);


/**
 * Returns a vector that is a sublist of triggers used in the run number given.
 */
vector<Trigger*>* GetTriggerSubvector(const int runNumber);


/**
 * Returns a list of trigger efficiencies.
 */
vector<TF1*>* GetTriggerEfficiencyFunctions();


/**
 * Returns a list of trigger prescale corrected luminosities.
 * Format: lumi(trig, run #) = getTriggerLuminosities()[(run # index) + (trigger index)*(number of runs)]
 * Returned luminosity is in units of ub^-1.
 * To delete, need to run a loop over the 1st index, e.g.
 *   for (int i = 0; i < numtrig; i++) {
 *    delete [] triggerLuminosities[i];
 *   }
 *   delete [] triggerLuminosities;
 */
double** GetTriggerLuminosities();


/**
 * Stores the best triggers for the run referenced by rnIndex in kinematicTriggerVec.
 */
void SetBestTriggers(const int rnIndex);


/**
 * Returns true if this dataSet is listed in the run list.
 */
bool GetIsRun (const int dataSet);


/**
 * Returns a data type string based on dType.
 */
const char* GetDataStr (const short dType);


/**
 * Returns a period string based on pType.
 */
const char* GetPeriodStr (const short pType);


/**
 * Returns a file identifier for the dataset specified.
 */
TString GetIdentifier (const int dataSet, const bool isMC, const bool pA);


/**
 * Creates a new TFile that points at the relevant data set.
 */
TFile* GetTFile(const TString fileIdentifier);


/**
 * Initializes triggers complete with momentum and pseudorapidity cutoffs.
 * Also (optionally) initializes a trigger selection scheme for each kinematic bin in a given run.
 * This includes calculating the relevant luminosity based on the selection method:
 dealRpPbAnalysisHist*    Lumi (run #, pt, eta) = Lumi (Trigger (pt, eta), run #)
 * where Trigger (pt, eta) is a function defining the selection scheme.
 *
 * IMPORTANT: Note that initialize should be called as early as possible to ensure that all global variables are in working order.
 */
void Initialize (const int dataSetId, const bool isMC, const bool initTriggerMaps=true);

} // end namespace

#endif
