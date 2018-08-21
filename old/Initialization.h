#include <vector>
#include <GlobalParams.h>
#include <Trigger.h>

#ifndef __Initialization_h__
#define __Initialization_h__

/**
 * This class initializes important jet trigger lists and used during these analyses.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

using namespace std;

vector<Trigger*> triggerVec(0); // total list of triggers used.
const int trigThres = 0; // Additional threshold requirement for triggers
const string minbiasTriggerName = "HLT_mb_mbts_L1MBTS_1";


/**
 * Returns true iff trigName is already a trigger in triggerVec.
 */
bool InTriggerVec(const string trigName);


/**
 * Initializes triggers complete with momentum and pseudorapidity cutoffs.
 * Triggers are stored in a vector TriggerVec.
 */
void SetupTriggers (const int runNumber=0);

#endif
