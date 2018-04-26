/**
 * This class initializes important jet trigger lists and used during these analyses.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */


#include "Trigger.C"

using namespace std;

vector<Trigger*> triggerVec(0); // total list of triggers used.
const int trigThres = 0; // Additional threshold requirement for triggers
const string minbiasTriggerName = "HLT_mb_mbts_L1MBTS_1";


/**
 * Returns true iff trigName is already a trigger in triggerVec.
 */
bool inTriggerVec(const string trigName) {
  for (Trigger* trig : triggerVec) {
    if (trig->name == trigName) {
      return true;
    }
  }
  return false;
}


/**
 * Initializes triggers complete with momentum and pseudorapidity cutoffs.
 * Triggers are stored in a vector TriggerVec.
 */
void setupTriggers (const int runNumber=0) {
  const string triggerListTxt = "/Users/jeffouellette/Research/atlas-hi/pPbJetTriggerList.txt";
  ifstream triggerListFile (triggerListTxt.c_str());
  if (!triggerListFile.is_open()) {
    cout << "Status: In Initialization (37): Cannot find pPbJetTriggerList.txt!" << endl;
    throw runtime_error("ifstream file not found");
  }

  string trigName, referenceTriggerName;
  double trigLetaFloat, trigUetaFloat; 
  int trigMinPt, trigThresholdPt, trigLowerRunNumber, trigUpperRunNumber = 0;

  Trigger* minbiasTrigger = new Trigger(minbiasTriggerName, 0, -4.9, 4.9, 0, INT_MAX);
  minbiasTrigger->referenceTrigger = minbiasTrigger;
  minbiasTrigger->disabled = false;
  triggerVec.push_back(minbiasTrigger);

  while (triggerListFile) {
    triggerListFile >> trigName;
    triggerListFile >> trigThresholdPt; // listed threshold pt for trigger - e.g. 15 for a j15 trigger
    triggerListFile >> trigLetaFloat;
    triggerListFile >> trigUetaFloat;
    triggerListFile >> referenceTriggerName;
    triggerListFile >> trigLowerRunNumber;
    triggerListFile >> trigUpperRunNumber;
    triggerListFile >> trigMinPt; // pt required above turn on for efficiency ~1
    assert (trigUetaFloat >= trigLetaFloat);

    const bool isPhysicsTrigger = (trigLowerRunNumber <= runNumber && runNumber < trigUpperRunNumber) || (runNumber == 0);

    /**
    *** If we want to use the trigger for this run, then
    *** disabled = [if this trigger is not an active physics trigger for the run]
    **/

    if (!inTriggerVec(trigName)) {
      Trigger* newTrig = new Trigger(trigName, trigThresholdPt, trigLetaFloat, trigUetaFloat, !isPhysicsTrigger, trigLowerRunNumber, trigUpperRunNumber, trigThresholdPt+trigMinPt+trigThres);
      triggerVec.push_back(newTrig);

      if (referenceTriggerName == "MINBIAS") {
        newTrig->referenceTrigger = minbiasTrigger;
      }
      else {
        for (Trigger* trig : triggerVec) {
          if (trig->name == referenceTriggerName) {
            newTrig->referenceTrigger = trig;
            break;
          }
        }
      }
    }
  }
  if (debugStatements) {
    for (Trigger* trig : triggerVec) {
      cout << "Status: In Initialization.C (87): Trig name:\t" << trig->name << "\tReference trig:\t" << trig->referenceTrigger->name << endl;
    }
  }
  triggerListFile.close();
  numtrigs = triggerVec.size();

  /**** Assign indices for tree branching. ****/
  for (int i = 0; i < numtrigs; i++) {
    triggerVec[i]->index = i;
  }
  /**** End assign indices ****/
  return;
}
