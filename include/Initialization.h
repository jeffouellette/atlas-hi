#ifndef __Initialization_h__
#define __Initialization_h__

#include <vector>
#include "Trigger.h"

using namespace std;

namespace atlashi {

/**
 * This class initializes important jet trigger lists and used during these analyses.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

extern vector<Trigger*> triggerVec; // total list of triggers used.
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

} // end namespace

#endif
