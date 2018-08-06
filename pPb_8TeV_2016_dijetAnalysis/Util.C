#include <vector>
#include <fstream>
#include <iostream>
#include "Params.C"
#include "../Initialization.C"
//#include <TFile.h>
//#include <TTree.h>
//#include <TLorentzVector.h>
//#include <TF1.h>
//#include <TSystemDirectory.h>
//#include <TSystemFile.h>

//=========================================================================================================
/** Variable declarations **/

// Trigger vectors used for pt and eta binning
double kinematicLumiVec[numpbins][numetabins]; // total exposed luminosity in a kinematic bin. Units: ub^-1
Trigger* kinematicTriggerVec[numpbins][numetabins]; // best trigger to use in a kinematic bin.


//=========================================================================================================
// General functions


/**
 * Returns xp for the event.
 */
static double get_xp(double jpt0, double jpt1, double jeta0, double jeta1, bool periodA) {
  double prefactor = TMath::Sqrt(Z/A) / sqrt_s_nn;
  if (!periodA) return prefactor * (jpt0 * TMath::Exp(jeta0) + jpt1 * TMath::Exp(jeta1));
  return prefactor * (jpt0 * TMath::Exp(-jeta0) + jpt1 * TMath::Exp(-jeta1));
}


/**
 * Returns xa for the event.
 */
static double get_xa(double jpt0, double jpt1, double jeta0, double jeta1, bool periodA) {
  double prefactor = TMath::Sqrt(A/Z) / sqrt_s_nn;
  if (!periodA) return prefactor * (jpt0 * TMath::Exp(-jeta0) + jpt1 * TMath::Exp(-jeta1));
  return prefactor * (jpt0 * TMath::Exp(jeta0) + jpt1 * TMath::Exp(jeta1));
}


/**
 * Returns the momentum transfer ("hardness") Q for the event.
 */
static double get_q2(double xp, double je, double jpt) {
  return (double)TMath::Sqrt(A/Z)*sqrt_s_nn*xp*(je-TMath::Sqrt(je*je-jpt*jpt));
}


/**
 * Returns the dijet mass Mjj.
 */
static double get_mjj(TLorentzVector jet0, TLorentzVector jet1) {
  return (jet0+jet1).Mag();
}


/**
 * Returns the bin number of pt in pbins.
 * Returns -1 iff pt < pbins[0] or (numpbins+1) iff pbins[numpbins] <= pt.
 * Otherwise returns n iff pbins[n] <= jpt and jpt < pbins[n+1], etc.
 */
static int getPbin (const double pt) {
  int bin = 0;
  while (bin < numpbins+1 && pbins[bin] < pt) bin++;
  return bin - 1;
}


/**
 * Returns the bin number of eta in etabins.
 * Returns -1 if eta < etabins[0] and numetabins if etabins[numetabins] <= eta.
 * Otherwise returns n iff etabins[n] <= eta && eta < etabins[n+1], etc.
 */
static int getEtabin (const double eta) {
  int bin = 0;
  while (bin < numetabins+1 && etabins[bin] < eta) bin++;
  return bin - 1;
}


//bool isPeriodA (const int dataset, const bool isMC) {
//  if (useDataVersion == 0) {
//    TSystemDirectory dir(dataPath.c_str(), dataPath.c_str());
//    TList* sysfiles = dir.GetListOfFiles();
//    if (sysfiles) {
//      TSystemFile* sysfile;
//      TString fname;
//      TIter next(sysfiles);
//
//      while ((sysfile = (TSystemFile*)next())) {
//        fname = sysfile->GetName();
//        if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
//          if (fname.Contains(to_string(dataset))) {
//            return fname.Contains("e6394") || fname.Contains("e6518");
//            //pa = fname.Contains("e6395") || fname.Contains("e6519");
//          }
//        }
//      }
//    }
//    delete sysfiles; // this line feels wrong to write...
//  }
//  else {
//      return dataset < 313500;
//  }
//  return false;
//}


/**
 * Determines whether analysis of this run should be skipped.
 */
static bool SkipRun (const bool pA) {
  if (!runData) return true;
  if (pA && !runPeriodA) return true;
  if (!pA && !runPeriodB) return true;
  return false;
}


/**
 * Determines whether analysis of this MC sample should be skipped.
 */
bool SkipMC (const bool pA) {
  if (!runMC) return true;
  if (pA && !runPeriodA) return true;
  if (!pA && !runPeriodB) return true;
  return false;
}


/**
 * Returns a vector that is a sublist of triggers used in the run number given.
 */
vector<Trigger*>* GetTriggerSubvector(const int runNumber) {
  vector<Trigger*>* triggerSublist = new vector<Trigger*>(0);
  for (Trigger* trig : triggerVec) {
   if (trig->lowerRunNumber <= runNumber && runNumber < trig->upperRunNumber) triggerSublist->push_back(trig);
  }
  return triggerSublist;
}


/**
 * Returns a list of trigger efficiencies.
 */
vector<TF1*>* GetTriggerEfficiencyFunctions() {
  TFile* thisFile = new TFile((effPath+"allEfficiencyHistograms.root").c_str(), "READ");
  vector<TF1*>* triggerEfficiencyFunctions = new vector<TF1*>(numtrigs);
  for (Trigger* trig : triggerVec) {
   int index = trig->index; 
   if (trig->referenceTrigger == trig) {
    (*triggerEfficiencyFunctions)[index] = new TF1(Form("%s_%s", trig->name.c_str(), fittedFunctionType.c_str()), "1", pbins[0], pbins[numpbins]);
    continue;
   }
   TGraphAsymmErrors* thisGraph = (TGraphAsymmErrors*)thisFile->Get(Form("%s_t_graph", trig->name.c_str()));
   (*triggerEfficiencyFunctions)[index] = (TF1*)(thisGraph->GetFunction(Form("%s_%s", trig->name.c_str(), fittedFunctionType.c_str())));
  }
  return triggerEfficiencyFunctions;
}


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
double** GetTriggerLuminosities() {

  double** triggerLuminosities = new double*[numtrigs];
  for (int trig_itr = 0; trig_itr < numtrigs; trig_itr++) {
   triggerLuminosities[trig_itr] = new double[numruns];
  }

  ifstream luminositiesTxt((homePath + "pPbJetTriggerLuminosities.txt").c_str());

  string thisTriggerName, thisLumiUnits;
  int thisRN;
  double thisLumi, thisConversionFactor;
  Trigger* thisTrigger;
  getline(luminositiesTxt, thisTriggerName); // skip the first line (table layout)

  while(luminositiesTxt) {
    thisTrigger = NULL;
    luminositiesTxt >> thisTriggerName;
    luminositiesTxt >> thisLumiUnits;
    for (Trigger* trig : triggerVec) {
      if (trig->name == thisTriggerName) {
        thisTrigger = trig;
        break;
      }
    }
    if (thisTrigger == NULL) {
      if (debugStatements) cout << "Warning: In Util.C (277): " << thisTriggerName << " is not a registered trigger! Ignoring." << endl;
      for (int rn_itr = 0; rn_itr < 60; rn_itr++) {
       luminositiesTxt >> thisLumi;
      } // just read through the rest of the line.
      continue;
    }

    if (thisLumiUnits == "ub") thisConversionFactor = 1e0;
    else if (thisLumiUnits == "mb") thisConversionFactor = 1e-3;
    else thisConversionFactor = 1e3;

    for (int rn_itr = 0; rn_itr < 30; rn_itr++) {
      luminositiesTxt >> thisRN;
      luminositiesTxt >> thisLumi;
      if (debugStatements) cout << "Status: In Util.C (289): adding " << thisTriggerName << "\tRun: " << thisRN << "\tLumi: " << thisLumi*thisConversionFactor << " ub^-1" << endl;
      if (!SkipRun(thisRN<313500)) { // ensures that thisRN is in the list of run numbers
        int thisRNIndex = 0;
        while (runNumbers[thisRNIndex] != thisRN) thisRNIndex++; // don't need to check array bounds since we are not skipping thisRN
        if (thisRNIndex >= numruns) {
          cout << "Error: In Util.C (294): run number index runs past bounds of array!" << endl;
          throw out_of_range("Index out of bounds");
        }
        triggerLuminosities[thisTrigger->index][thisRNIndex] = thisLumi * thisConversionFactor;
      }
    }
  }
  luminositiesTxt.close();
  return triggerLuminosities;
}


/**
 * Stores the best triggers for the run referenced by rnIndex in kinematicTriggerVec.
 */
void SetBestTriggers(const int rnIndex) {
  /**
  *** Now calculate the best trigger to use for each bin.
  *** The best trigger is the one satisfying the kinematic cuts of the bin with the highest luminosity.
  **/
  //vector<int>* runNumbers = getRunNumbers();
  const int thisRunNumber = runNumbers[rnIndex];
  //const int numruns = runNumbers->size();
  double** totalLumiVec = GetTriggerLuminosities(); // luminosity a given trigger sees at a given pbin, etabin in a given run

  // Get list of triggers used in this run
  vector<Trigger*>* triggerSubvector = GetTriggerSubvector(thisRunNumber);

  short actetabin;
  double bestTrigLumi, thisTrigLumi;
  Trigger* bestTrigger;
  for (short etabin = 0; etabin < numetabins; etabin++) { 
    if (thisRunNumber < 313500) actetabin = numetabins - etabin - 1;
    else actetabin = etabin;
    bestTrigger = NULL;
    bestTrigLumi = 0;

    for (short pbin = 0; pbin < numpbins; pbin++) {
      for (Trigger* trig : (*triggerSubvector)) {
        if (!(trig->lower_eta <= etabins[etabin] && etabins[etabin+1] <= trig->upper_eta && trig->min_pt <= pbins[pbin])) continue;
        thisTrigLumi = totalLumiVec[trig->index][rnIndex];
        if ((bestTrigger == NULL) || (bestTrigLumi < thisTrigLumi) || (bestTrigLumi == thisTrigLumi && bestTrigger->min_pt < trig->min_pt)) {
          bestTrigLumi = thisTrigLumi;
          bestTrigger = trig;
        }
      }
      kinematicTriggerVec[pbin][etabin] = bestTrigger;
    }
  }
  delete[] totalLumiVec;
  return;
}


/**
 * Returns true if this dataSet is listed in the run list.
 */
bool GetIsRun (const int dataSet) {
  for (int rn : runNumbers) {
    if (dataSet == rn) return true;
  }
  return false;
}


/**
 * Returns a data type string based on dType.
 */
const char* GetDataStr (const short dType) {
  if (dType == 0) return "data";
  else return "mc";
}

/**
 * Returns a period string based on pType.
 */
const char* GetPeriodStr (const short pType) {
  switch (pType) {
   case 0: return "periodA";
   case 1: return "periodB";
   case 2: return "periodAB";
  }
  return "";
}

/**
 * Returns a file identifier for the dataset specified.
 */
TString GetIdentifier (const int dataSet, const bool isMC, const bool pA) {
  if (!isMC) return TString(dataSet);
  else {
   if (dataSet == 6) { //then it is a JZ2 validation sample
    if (pA) return "jetjet.valid.JZ2R04.pPb";
    else return "jetjet.valid.JZ2R04.Pbp";
   }
   TString id = Form("jetjet.JZ%iR04.", dataSet);
   if (pA) {
    id = id + "pPb"; // period A MC samples have Pbp in the file id
   }
   else {
    id = id + "Pbp"; // likewise period B MC samples have pPb in the file id
   }
   return id;
  }
}


/**
 * Creates a new TFile that points at the relevant data set.
 */
TFile* GetTFile(const TString fileIdentifier) {
  TSystemDirectory dir(dataPath.c_str(), dataPath.c_str());
  TList* sysfiles = dir.GetListOfFiles();
  TFile* file = NULL;
  if (sysfiles) {
    TSystemFile* sysfile;
    TString fname;
    TIter next(sysfiles);

    while ((sysfile = (TSystemFile*)next())) {
      fname = sysfile->GetName();
      if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
        if (debugStatements) cout << "Status: In Util::GetTFile: Found " << fname.Data() << endl; 
        if (fname.Contains(fileIdentifier)) {
          file = new TFile(dataPath+fname, "READ");
          if (!file) {
           cout << "Error: In Util::GetTFile: Cannot create TFile." << endl;
           break;
          }
          break;
        }
      }
    }
    if (!file)
     cout << "Warning: In Util::GetTFile: Did not find file matching specificiation." << endl;
  }
  else {
    cout << "Error: In Util::GetTFile: Cannot obtain list of files." << endl;
  }
  if (sysfiles) delete sysfiles;
  return file;
}


/**
 * Initializes triggers complete with momentum and pseudorapidity cutoffs.
 * Also (optionally) initializes a trigger selection scheme for each kinematic bin in a given run.
 * This includes calculating the relevant luminosity based on the selection method:
 dealRpPbAnalysisHist*    Lumi (run #, pt, eta) = Lumi (Trigger (pt, eta), run #)
 * where Trigger (pt, eta) is a function defining the selection scheme.
 *
 * IMPORTANT: Note that initialize should be called as early as possible to ensure that all global variables are in working order.
 */
void Initialize (const int dataSetId, const bool isMC, const bool initTriggerMaps=true) {

  //assert (useDataVersion == 8 || useDataVersion == 0);
  if (debugStatements) cout << Form("Status: In Util.C (362): Initializing trigger system for dataSetId %i...", dataSetId) << endl;

  setupDirectories(Form("v%i/", useDataVersion), "pPb_8TeV_2016_dijetAnalysis/");

  /**** If Monte Carlo (MC), we don't need to get trigger information. ****/
  if (isMC) return;

  /**** Create an array of triggers ****/
  SetupTriggers(dataSetId);  

  /**** Debugging statements ****/
  if (debugStatements) {
   cout << "Status: In Util.C (380): Triggers initialized for dataSetId " << dataSetId << endl;
   cout << "Status: In Util.C (381): Num trigs = " << numtrigs << endl;
   cout << "Status: In Util.C (382): Trigger pt bin path is " << trigPath << endl;    
   cout << "Status: In Util.C (383): Saving plots to " << plotPath << endl;
  }
  /**** End debugging ****/


  /**** Instantiate the pseudorapidity interval trigger vectors if required. ****/
  if (initTriggerMaps) {

    if (debugStatements) cout << "Status: In Util.C (391): Generating triggering scheme..." << endl;

    /**** Local variable array declarations ****/
    //vector<int>* runNumbers = getRunNumbers();
    //const int numruns = runNumbers->size();
    double** totalLumiVec = GetTriggerLuminosities(); // luminosity a particular trigger sees at a given pbin, etabin in a given run
    for (short pbin = 0; pbin < numpbins; pbin++) {
     for (short etabin = 0; etabin < numetabins; etabin++) {
      kinematicLumiVec[pbin][etabin] = 0;
      kinematicTriggerVec[pbin][etabin] = NULL;
     }
    }
    /**** End local variable declarations ****/


    /**
    *** Combine trigger data from all runs into one histogram for each trigger. 
    *** If the trigger never fired in a run, assume it wasn't on so don't add its luminosity.
    **/

    /**** Loop over runs to find the best trigger for each phase space bin, then add its luminosity to the effective luminosity in that bin ****/
    for (int thisRunNumber : runNumbers) {

     if (SkipRun(thisRunNumber<313500)) { // Only allow desired runs to be considered
      continue;
     }
     int thisRNIndex = 0;
     while (runNumbers[thisRNIndex] != thisRunNumber) thisRNIndex++; // don't have to check array bounds since we are not skipping this run
  
     vector<Trigger*>* triggerSubvector = GetTriggerSubvector(thisRunNumber); 

     /**
     *** Now calculate the best trigger to use for each bin.
     *** The best trigger is the one satisfying the kinematic cuts of the bin with the highest luminosity.
     **/
     int actetabin;
     double bestTrigLumi, thisTrigLumi;
     Trigger* bestTrigger;
     for (int etabin = 0; etabin < numetabins; etabin++) { 
      if (thisRunNumber < 313500) actetabin = numetabins - etabin - 1;
      else actetabin = etabin;
      bestTrigger = NULL;
      bestTrigLumi = 0;

      for (int pbin = 0; pbin < numpbins; pbin++) {
       for (Trigger* trig : (*triggerSubvector)) {
        if (!(trig->lower_eta <= etabins[etabin] && etabins[etabin+1] <= trig->upper_eta && trig->min_pt <= pbins[pbin])) continue;
        thisTrigLumi = totalLumiVec[thisRNIndex][trig->index];
        if ((bestTrigger == NULL) || (bestTrigLumi < thisTrigLumi) || (bestTrigLumi == thisTrigLumi && bestTrigger->min_pt < trig->min_pt)) {
         bestTrigLumi = thisTrigLumi;
         bestTrigger = trig;
        }
       }
       kinematicLumiVec[pbin][actetabin] += bestTrigLumi;
       if (dataSetId == thisRunNumber) {
        kinematicTriggerVec[pbin][etabin] = bestTrigger;
       }
      }
     }
    }
    /**** End loop over runs ****/


    /**** Debug statements ****/
    if (debugStatements) cout << "Status: In Util.C (453): Triggering scheme completed for run number " << dataSetId << endl;  
    if (dataSetId == 313063 && debugStatements) {
      cout << "Status: In Util.C (455): Example trigger assignment (run " << dataSetId << "):" << endl << endl;
      for (int etabin = 0; etabin < numetabins; etabin++) {
        int actetabin;
        if (dataSetId < 313500) actetabin = numetabins - etabin - 1;
        else actetabin = etabin;
        for (int pbin = 0; pbin < numpbins; pbin++) {
          if (kinematicTriggerVec[pbin][etabin] != NULL) cout << Form("etabin=%i,\tpbin=%i, \teta=(%.1f, %.1f),\tp=(%i, %i),\ttrig=%s,\tlumi=%f", etabin, pbin, etabins[etabin], etabins[etabin+1], (int)pbins[pbin], (int)pbins[pbin+1], kinematicTriggerVec[pbin][etabin]->name.c_str(), kinematicLumiVec[pbin][actetabin]) << endl;
        }
      } 
    }
    /**** End debugging ****/

    /**** Memory cleanup ****/
    delete [] totalLumiVec;
    /**** End memory cleanup ****/
  }
  /**** End instantiate trigger maps ****/

  return;
}
