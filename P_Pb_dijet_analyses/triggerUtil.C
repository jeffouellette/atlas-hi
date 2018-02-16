#include <vector>
#include <fstream>
#include <iostream>
#include "trigger.C"

using namespace std;

//=========================================================================================================
/** Variable declarations **/

// Trigger vectors used for pt and eta binning
vector<Trigger*> triggerVec(0); // total list of triggers used.
double kinematicLumiVec[numpbins*numetabins]; // total exposed luminosity in a kinematic bin. Units: ub^-1
Trigger* kinematicTriggerVec[numpbins*numetabins]; // best trigger to use in a kinematic bin.
//double kinematicEfficiencyVec[numpbins*numetabins]; // best trigger efficiency in a kinematic bin.


//=========================================================================================================
// General functions

/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
static double* linspace(double lo, double hi, int num) {
    double* arr = new double[num+1];
    double delta = ((double)(hi)-(double)(lo))/(double)(num);
    for (int i = 0; i <= num; i++) {
        arr[i] = lo + i * delta;
    }
    return arr;
}


/**
 * Returns a logarithmically spaced array, where the 0th element is lo and the num-th element is hi.
 */
static double* logspace(double lo, double hi, int num) {
    double loghi = TMath::Log2(hi);
    if (lo == 0) {
        double* arr = linspace(TMath::Log2(hi/(100*num)), loghi, num);
        for (int i = 0; i <= num; i++) {
            arr[i] = TMath::Power(2, arr[i]);
        }
        return arr;
    } else {
        double loglo = TMath::Log2(lo);
        double* arr = linspace(loglo, loghi, num);
        for (int i = 0; i <= num; i++) {
            arr[i] = TMath::Power(2, arr[i]);
        }
        return arr;
    }
}


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
static int getPbin (double pt) {
    int bin = 0;
    while (bin < numpbins+1 && pbins[bin++] < pt);
    return bin - 2;
}


/**
 * Returns the bin number of eta in etabins.
 * Returns -1 if eta < etabins[0] and (numetabins+1) if etabins[numetabins] <= eta.
 * Otherwise returns n iff etabins[n] <= eta && eta < etabins[n+1], etc.
 */
static int getEtabin (double eta) {
    int bin = 0;
    while (bin < numetabins+1 && etabins[bin++] < eta);
    return bin - 2;
}


/**
 * Determines whether analysis of this run should be skipped.
 */
static bool skipRun (int rn) {
    if (rn < 313500 && !runPeriodA) return true;
    if (rn > 313500 && !runPeriodB) return true;
    bool contains_rn = false;
    int i = 0;
    switch (useDataVersion) {
        case 7: {
            while (i < sizeof(run_list_v7)/sizeof(int) && !contains_rn) {
                contains_rn = run_list_v7[i] == rn;
                i++;
            }
            break;
        }
    }
    return !contains_rn;
}

/**
 * Returns a pointer to an std:vector<int> of the run numbers currently being processed.
 */ 
static vector<int>* getRunNumbers() {
    vector<int>* rns = new vector<int>(0);
    switch (useDataVersion) {
        case 7: {
            for (int i = 0; i < sizeof(run_list_v7)/sizeof(int); i++) {
                rns->push_back(run_list_v7[i]);
            }
            break;
        }
    }
    return rns;
}


/**
 * Returns true iff trigName is already a trigger in triggerVec.
 */
bool inTriggerVec(string trigName) {
    for (Trigger trig : triggerVec) {
        if (trig.name == trigName) {
            return true;
        }
    }
    return false;
}


/**
 * Returns a vector that is a sublist of triggers used in the run number given.
 */
vector<Trigger*>* getTriggerSubvector(int runNumber) {
    vector<Trigger*>* triggerSublist = new vector<Trigger*>(0);
    for (Trigger* trig : triggerVec) {
        if (trig->lowerRunNumber <= runNumber && runNumber < trig->upperRunNumber) triggerSublist->push_back(trig);
    }
    return triggerSublist;
}


/**
 * Returns a list of trigger efficiencies.
 */
vector<TF1*>* getTriggerEfficiencyFunctions() {

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
 */
double* getTriggerLuminosities() {

    vector<int>* runNumbers = getRunNumbers();
    const int numruns = runNumbers->size();
    double* triggerLuminosities = new double[numtrigs*numruns];
    ifstream luminositiesTxt((workPath + "luminosities.txt").c_str());

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
            if (debugStatements) cout << "Warning: In triggerUtil.C (200): " << thisTriggerName << " is not a registered trigger! Ignoring." << endl;
            for (int rn_itr = 0; rn_itr < 60; rn_itr++) { luminositiesTxt >> thisLumi; } // just read through the rest of the line.
            continue;
        }

        if (thisLumiUnits == "ub") thisConversionFactor = 1e0;
        else if (thisLumiUnits == "mb") thisConversionFactor = 1e-3;
        else thisConversionFactor = 1e3;

        for (int rn_itr = 0; rn_itr < 30; rn_itr++) {
            luminositiesTxt >> thisRN;
            luminositiesTxt >> thisLumi;
            if (debugStatements) cout << "Status: In triggerUtil.C (212): adding " << thisTriggerName << "\tRun: " << thisRN << "\tLumi: " << thisLumi*thisConversionFactor << " ub^-1" << endl;
            if (!skipRun(thisRN)) { // ensures that thisRN is in the list of run numbers
                int thisRNIndex = 0;
                while ((*runNumbers)[thisRNIndex] != thisRN) thisRNIndex++; // don't need to check array bounds since we are not skipping thisRN
                if (thisRNIndex >= numruns) {
                    cout << "Error: In triggerUtil.C (217): run number index runs past bounds of array!" << endl;
                    throw out_of_range("Index out of bounds");
                }
                triggerLuminosities[thisRNIndex + (thisTrigger->index)*numruns] = thisLumi * thisConversionFactor;
            }
        }
    }
    luminositiesTxt.close();
    delete runNumbers;
    return triggerLuminosities;
}


/**
 * Stores the best triggers for the run referenced by rnIndex in kinematicTriggerVec.
 */
void setBestTriggers(int rnIndex) {
    /**
    *** Now calculate the best trigger to use for each bin.
    *** The best trigger is the one satisfying the kinematic cuts of the bin with the highest luminosity.
    **/
    vector<int>* runNumbers = getRunNumbers();
    const int thisRunNumber = (*runNumbers)[rnIndex];
    const int numruns = runNumbers->size();
    double* totalLumiVec = getTriggerLuminosities(); // luminosity a given trigger sees at a given pbin, etabin in a given run

    // Get list of triggers used in this run
    vector<Trigger*>* triggerSubvector = getTriggerSubvector(thisRunNumber);

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
                thisTrigLumi = totalLumiVec[rnIndex + (trig->index)*numruns];
                if ((bestTrigger == NULL) || (bestTrigLumi < thisTrigLumi) || (bestTrigLumi == thisTrigLumi && bestTrigger->min_pt < trig->min_pt)) {
                    bestTrigLumi = thisTrigLumi;
                    bestTrigger = trig;
                }
            }
            kinematicTriggerVec[pbin + etabin*numpbins] = bestTrigger;
        }
    }
    delete[] totalLumiVec;
    delete runNumbers;
    return;
}


/**
 * Initializes triggers complete with momentum and pseudorapidity cutoffs.
 * Also (optionally) initializes a trigger selection scheme for each kinematic bin in a given run.
 * This includes calculating the relevant luminosity based on the selection method:
 dealRpPbAnalysisHist*    Lumi (run #, pt, eta) = Lumi (Trigger (pt, eta), run #)
 * where Trigger (pt, eta) is a function defining the selection scheme.
 */
void initialize (int runNumber=0, bool initTriggerMaps=true) {

    assert (useDataVersion >= 7);
    if (debugStatements) cout << Form("Status: In triggerUtil.C (248): Initializing trigger system for run %i...", runNumber) << endl;

    /**** Reset directory information for correct versioning ****/ 
    {
        string versionString;
        versionString = "v" + to_string(useDataVersion) + "/";
        rootPath = rootPath + versionString;
        dataPath = dataPath + versionString;
        plotPath = plotPath + versionString;
        ptPath = rootPath + "ptData/";
        trigPath = rootPath + "trigData/";
        effPath = rootPath + "effData/";
        xPath = rootPath + "xData/";
        RpPbPath = rootPath + "RpPbData/";
    }
    /**** End reset direction information ****/


    /**** Create an array of triggers ****/
    {
        string triggerListTxt = workPath + "fullTriggerList.txt";
        ifstream triggerListFile (triggerListTxt.c_str());
        if (!triggerListFile.is_open()) {
            cout << "Status: In triggerUtil (277): Cannot find triggerList.txt!" << endl;
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
            
            bool isPhysicsTrigger = (trigLowerRunNumber <= runNumber && runNumber < trigUpperRunNumber) || (runNumber == 0);

            /**
            *** If we want to use the trigger for this run, then
            *** disabled = [if this trigger is not an active physics trigger for the run]
            **/

            if (!inTriggerVec(trigName)) {
                Trigger* newTrig = new Trigger(trigName, trigThresholdPt, trigLetaFloat, trigUetaFloat, !isPhysicsTrigger, trigLowerRunNumber, trigUpperRunNumber, trigThresholdPt+trigMinPt);
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
                cout << "Status: In triggerUtil.C (327): Trig name:\t" << trig->name << "\tReference trig:\t" << trig->referenceTrigger->name << endl;
            }
        }
        triggerListFile.close();
        numtrigs = triggerVec.size();

        /**** Assign indices for tree branching. ****/
        for (int i = 0; i < numtrigs; i++) {
            triggerVec[i]->index = i;
        }
        /**** End assign indices ****/
    }
    /**** End create a list of triggers ****/


    /**** Debugging statements ****/
    if (debugStatements) {
        cout << "Status: In triggerUtil.C (344): Triggers initialized for run " << runNumber << endl;
        cout << "Status: In triggerUtil.C (345): Num trigs = " << numtrigs << endl;
        cout << "Status: In triggerUtil.C (346): Trigger pt bin path is " << trigPath << endl;    
        cout << "Status: In triggerUtil.C (347): Saving plots to " << plotPath << endl;
    }
    /**** End debugging ****/


    /**** Instantiate the pseudorapidity interval trigger vectors if required. ****/
    if (initTriggerMaps) {

        if (debugStatements) cout << "Status: In triggerUtil.C (355): Generating triggering scheme..." << endl;

        /**** Local variable array declarations ****/
        vector<int>* runNumbers = getRunNumbers();
        const int numruns = runNumbers->size();
        double* totalLumiVec = getTriggerLuminosities(); // luminosity a particular trigger sees at a given pbin, etabin in a given run
        for (int n = 0; n < numpbins*numetabins; n++) {
            kinematicLumiVec[n] = 0;
//            kinematicEfficiencyVec[n] = 0.;
            kinematicTriggerVec[n] = NULL;
        }
        /**** End local variable declarations ****/


        /**
        *** Combine trigger data from all runs into one histogram for each trigger. 
        *** If the trigger never fired in a run, assume it wasn't on so don't add its luminosity.
        **/

        /**** Loop over runs to find the best trigger for each phase space bin, then add its luminosity to the effective luminosity in that bin ****/
        for (int thisRunNumber : *runNumbers) {

            if (skipRun(thisRunNumber)) { // Only allow desired runs to be considered
                continue;
            }
            int thisRNIndex = 0;
            while ((*runNumbers)[thisRNIndex] != thisRunNumber) thisRNIndex++; // don't have to check array bounds since we are not skipping this run
   
            vector<Trigger*>* triggerSubvector = getTriggerSubvector(thisRunNumber); 

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
                        thisTrigLumi = totalLumiVec[thisRNIndex + (trig->index)*numruns];
                        if ((bestTrigger == NULL) || (bestTrigLumi < thisTrigLumi) || (bestTrigLumi == thisTrigLumi && bestTrigger->min_pt < trig->min_pt)) {
                            bestTrigLumi = thisTrigLumi;
                            bestTrigger = trig;
                        }
                    }
                    kinematicLumiVec[pbin + actetabin*numpbins] += bestTrigLumi;
                    if (runNumber == thisRunNumber) {
                        kinematicTriggerVec[pbin + etabin*numpbins] = bestTrigger;
                    }
                }
            }
        }
        /**** End loop over runs ****/


        /**** Debug statements ****/
        if (debugStatements) cout << "Status: In triggerUtil.C (436): Triggering scheme completed for run number " << runNumber << endl;  
        if (runNumber == 313063 && debugStatements) {
            cout << "Status: In triggerUtil.C (438): Example trigger assignment (run " << runNumber << "):" << endl << endl;
            for (int etabin = 0; etabin < numetabins; etabin++) {
                int actetabin;
                if (runNumber < 313500) actetabin = numetabins - etabin - 1;
                else actetabin = etabin;
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    if (kinematicTriggerVec[pbin + etabin*numpbins] != NULL) cout << Form("etabin=%i,\tpbin=%i, \teta=(%.1f, %.1f),\tp=(%i, %i),\ttrig=%s,\tlumi=%f", etabin, pbin, etabins[etabin], etabins[etabin+1], (int)pbins[pbin], (int)pbins[pbin+1], kinematicTriggerVec[pbin + etabin*numpbins]->name.c_str(), kinematicLumiVec[pbin + actetabin*numpbins]) << endl;
                }
            } 
        }
        /**** End debugging ****/

        /**** Memory cleanup ****/
        delete runNumbers;
        delete [] totalLumiVec;
        /**** End memory cleanup ****/
    }
    /**** End instantiate trigger maps ****/

    return;
}
