#include "Initialization.h"
#include "GlobalParams.h"
#include "Utilities.h"

#include <iostream>
#include <fstream>

namespace atlashi {

vector<Trigger*> triggerVec (0);

bool InTriggerVec(const string trigName) {
  for (Trigger* trig : triggerVec) {
    if (trig->name == trigName) {
      return true;
    }
  }
  return false;
}


void SetupTriggers (const int runNumber) {
  const string triggerListTxt = "/Users/jeffouellette/Research/atlas-hi/pPbJetTriggerList.txt";
  ifstream triggerListFile (triggerListTxt.c_str());
  if (!triggerListFile.is_open()) {
    cout << "Status: In Initialization (37): Cannot find pPbJetTriggerList.txt!" << endl;
    throw runtime_error("ifstream file not found");
  }

  string trigName, referenceTriggerName;
  double trigLetaFloat, trigUetaFloat; 
  int trigMinPt, trigThresholdPt, trigLowerRunNumber, trigUpperRunNumber = 0;

  Trigger* minbiasTrigger = new Trigger(minbiasTriggerName, 0, -4.9, 4.9, 0, 2147483647);
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

    if (!InTriggerVec(trigName)) {
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
  int i = 0;
  for (Trigger* trig : triggerVec) {
    trig->index = i++;
  }
  /**** End assign indices ****/
  return;
}

} // end namespace
