#include <vector>
#include <fstream>
#include <iostream>
#include "trigger.C"

using namespace std;

//=========================================================================================================
/** Variable declarations **/

// Global parameters
int runNumber;
bool periodA;
double totalLuminosity;

// Trigger vectors used for pt and eta binning
vector<Trigger*> triggerVec(0);
double kinematicLumiVec[numpbins*numetabins];
Trigger* kinematicTriggerVec[numpbins*numetabins];
double kinematicEfficiencyVec[numpbins*numetabins];

// Run list information
//const int run_list_v3[30] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313688, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170}; // full run list for future reference
const int run_list_v5[23]= {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313574, 313575, 313629, 313630, 313688, 313695, 313929, 313935, 313984, 314014/*, 314077, 314112*/, 314157, 314170};
const int run_list_v6[25] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313574, 313629, 313630, 313688, 313695, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170};


//=========================================================================================================
// General functions

/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
double* linspace(double lo, double hi, int num) {
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
double* logspace(double lo, double hi, int num) {
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
double get_xp(double jpt0, double jpt1, double jeta0, double jeta1, bool periodA) {
    double prefactor = TMath::Sqrt(Z/A) / sqrt_s_nn;
    if (!periodA) return prefactor * (jpt0 * TMath::Exp(jeta0) + jpt1 * TMath::Exp(jeta1));
    return prefactor * (jpt0 * TMath::Exp(-jeta0) + jpt1 * TMath::Exp(-jeta1));
}


/**
 * Returns xa for the event.
 */
double get_xa(double jpt0, double jpt1, double jeta0, double jeta1, bool periodA) {
    double prefactor = TMath::Sqrt(A/Z) / sqrt_s_nn;
    if (!periodA) return prefactor * (jpt0 * TMath::Exp(-jeta0) + jpt1 * TMath::Exp(-jeta1));
    return prefactor * (jpt0 * TMath::Exp(jeta0) + jpt1 * TMath::Exp(jeta1));
}


/**
 * Returns the momentum transfer ("hardness") Q for the event.
 */
double get_q2(double xp, double je, double jpt) {
    return (double)TMath::Sqrt(TMath::Sqrt(A/Z)*sqrt_s_nn*xp*(je-TMath::Sqrt(je*je-jpt*jpt)));
}


/**
 * Returns the dijet mass Mjj.
 */
double get_mjj(TLorentzVector jet0, TLorentzVector jet1) {
    return (jet0+jet1).Mag();
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
        case 6: {
            while (i < sizeof(run_list_v6)/sizeof(int) && !contains_rn) {
                contains_rn = run_list_v6[i] == rn;
                i++;
            }
            break;
        }
        case 5: {
            while (i < sizeof(run_list_v5)/sizeof(int) && !contains_rn) {
                contains_rn = run_list_v5[i] == rn;
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
        case 6: {
            for (int i = 0; i < sizeof(run_list_v6)/sizeof(int); i++) {
                rns->push_back(run_list_v6[i]);
            }
            break;
        }
        case 5: {
            for (int i = 0; i < sizeof(run_list_v5)/sizeof(int); i++) {
                rns->push_back(run_list_v5[i]);
            }
            break;
        }
    }
    return rns;
}


/**
 * Returns true iff trigName is already a trigger in triggerVec.
 */
bool in_triggerVec(string trigName) {
    for (Trigger trig : triggerVec) {
        if (trig.name == trigName) {
            return true;
        }
    }
    return false;
}


/**
 * Returns a list of trigger efficiencies.
 */
double* getTriggerEfficiencies() {

    TFile* thisFile = new TFile((effPath+"allEfficiencyHistograms.root").c_str(), "READ");
    double* triggerEfficiencies = new double[numtrigs*numpbins];
    for (Trigger* trig : triggerVec) {
        int index = trig->index; 
        TH1F* thisHist = (TH1F*)thisFile->Get(Form("%s_efficiency", trig->name.c_str()));
        for (int pbin = 0; pbin < numpbins; pbin++) {
            if (trig->min_pt + 50 < pbins[pbin]) triggerEfficiencies[index + pbin*numtrigs] = 1.;
            else triggerEfficiencies[index + pbin*numtrigs] = thisHist->GetBinContent(pbin+1);
        }
    }
    thisFile->Close();
    delete thisFile;
    return triggerEfficiencies;
}


/**
 * Returns a list of trigger prescale corrected luminosities.
 * Format: lumi(trig, run #, pbin, etabin) = getTriggerLuminosities()[(run # index) + ((trigger index) + (pbin + etabin*numpbins)*numtrigs)*numruns]
 */
double* getTriggerLuminosities() {

    vector<int>* runNumbers = getRunNumbers();
    const int numruns = runNumbers->size();
    double* triggerLuminosities = new double[numtrigs*numruns*numpbins*numetabins];
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
            if (debugStatements) cout << "Warning: In triggerUtil.C (347): " << thisTriggerName << " is not a registered trigger! Ignoring." << endl;
            for (int rn_itr = 0; rn_itr < 60; rn_itr++) { luminositiesTxt >> thisLumi; } // just read through the rest of the line.
            continue;
        }

        if (thisLumiUnits == "ub") thisConversionFactor = 1e-3;
        else if (thisLumiUnits == "mb") thisConversionFactor = 1e-6;
        else thisConversionFactor = 1.;

        for (int rn_itr = 0; rn_itr < 30; rn_itr++) {
            luminositiesTxt >> thisRN;
            luminositiesTxt >> thisLumi;
            if (debugStatements) cout << "Status: In triggerUtil.C (329): adding " << thisTriggerName << "\tRun: " << thisRN << "\tLumi: " << thisLumi*thisConversionFactor << " nb^-1" << endl;
            if (!skipRun(thisRN)) { // ensures that thisRN is in the list of run numbers
                int thisRNIndex = 0;
                while ((*runNumbers)[thisRNIndex] != thisRN) thisRNIndex++; // don't need to check array bounds since we are not skipping thisRN
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    for (int etabin = 0; etabin < numetabins; etabin++) {
                        if (thisTrigger->lower_eta <= etabins[etabin] && etabins[etabin+1] <= thisTrigger->upper_eta && thisTrigger->min_pt <= pbins[pbin]) {
                            triggerLuminosities[thisRNIndex + ((thisTrigger->index) + (pbin + etabin*numpbins)*numtrigs)*numruns] = thisLumi * thisConversionFactor;
                        }
                    }
                }
            }
        }
    }
    delete runNumbers;
    return triggerLuminosities;
}




/**
 * Initializes triggers complete with momentum and pseudorapidity cutoffs.
 * Also (optionally) initializes a trigger selection scheme for each kinematic bin in a given run.
 * This includes calculating the relevant luminosity based on the selection method:
 *    Lumi (run #, pt, eta) = Lumi (Trigger (pt, eta), run #)
 * where Trigger (pt, eta) is a function defining the selection scheme.
 */
void initialize (int rn=0, bool initTriggerMaps=true, bool initTriggerEfficiencies=true) {

    assert (useDataVersion == 5 || useDataVersion == 6);
    if (debugStatements) cout << Form("Status: In triggerUtil.C (176): Initializing trigger system for run %i...", rn) << endl;

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
    }
    /**** End reset direction information ****/


    /**** Store run information in global variables ****/
    runNumber = rn;
    if (rn < 313500) periodA = true;
    else periodA = false;
    /**** End store run information in global variables ****/


    /**** Create an array of triggers ****/
    {
        string triggerListTxt = workPath + "fullTriggerList.txt";
        ifstream triggerListFile (triggerListTxt.c_str());
        if (!triggerListFile.is_open()) {
            cout << "Status: In triggerUtil (189): Cannot find triggerList.txt!" << endl;
            throw runtime_error("ifstream file not found");
        }

        string trigName, referenceTriggerName;
        float trigPtFloat, trigLetaFloat, trigUetaFloat; 
        int trigLowerRunNumber, trigUpperRunNumber = 0;

        Trigger* minbiasTrigger = new Trigger(minbiasTriggerName, 0, -4.9, 4.9, 0, INT_MAX);
        minbiasTrigger->referenceTrigger = minbiasTrigger;
        minbiasTrigger->disabled = false;
        triggerVec.push_back(minbiasTrigger);

        while (triggerListFile) {

            triggerListFile >> trigName;
            triggerListFile >> trigPtFloat;
            triggerListFile >> trigLetaFloat;
            triggerListFile >> trigUetaFloat;
            triggerListFile >> referenceTriggerName;
            triggerListFile >> trigLowerRunNumber;
            triggerListFile >> trigUpperRunNumber;
            
            bool isPhysicsTrigger = (trigLowerRunNumber <= runNumber && runNumber < trigUpperRunNumber) || (runNumber == 0);

            /**
            *** If we want to use the trigger for this run, then
            *** iontrigger = [whether ion is in the name of the trigger]
            *** disabled = [if this trigger is not an active physics trigger for the run]
            **/

            if (!in_triggerVec(trigName)) {
                Trigger* newTrig = new Trigger(trigName, trigPtFloat, trigLetaFloat, trigUetaFloat, !isPhysicsTrigger, trigLowerRunNumber, trigUpperRunNumber);
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
                cout << "Status: In triggerUtil.C (242): Trig name:\t" << trig->name << "\tReference trig:\t" << trig->referenceTrigger->name << endl;
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
        cout << "Status: In triggerUtil.C (267): Triggers initialized for run " << rn << endl;
        cout << "Status: In triggerUtil.C (268): Num trigs = " << numtrigs << endl;
        cout << "Status: In triggerUtil.C (269): Trigger pt bin path is " << trigPath << endl;    
        cout << "Status: In triggerUtil.C (270): Saving plots to " << plotPath << endl;
    }
    /**** End debugging ****/


    /**** Instantiate the pseudorapidity interval trigger vectors if required. ****/
    if (initTriggerMaps) {

        if (debugStatements) cout << "Status: In triggerUtil.C (285): Generating triggering scheme..." << endl;

        /**** Local variable array declarations ****/
        const int numruns = getRunNumbers()->size();
        double* totalLumiVec = getTriggerLuminosities(); // luminosity a particular trigger sees at a given pbin, etabin in a given run
        for (int n = 0; n < numpbins*numetabins; n++) {
            kinematicLumiVec[n] = 0;
            kinematicEfficiencyVec[n] = 0.;
            kinematicTriggerVec[n] = NULL;
        }
        /**** End local variable declarations ****/


        /**
        *** Combine trigger data from all runs into one histogram for each trigger. 
        *** If the trigger never fired in a run, assume it wasn't on so don't add its luminosity.
        **/

        /**** Loop over root files for each run generated by trigger_pt_counts module ****/
        for (int thisRunNumber : *(getRunNumbers())) {

            if (skipRun(thisRunNumber)) { // Only allow desired runs to be considered
                continue;
            }
            int thisRNIndex = 0;
            while ((*getRunNumbers())[thisRNIndex] != thisRunNumber) thisRNIndex++; // don't have to check array bounds since we are not skipping this run
    
            vector<Trigger*> triggerSubList(0);
            for (Trigger* trig : triggerVec) {
                if (trig->lowerRunNumber <= thisRunNumber && thisRunNumber < trig->upperRunNumber && trig->name != minbiasTriggerName) {
                    triggerSubList.push_back(trig);
                }
            }

            /**
            *** Now calculate the best trigger to use for each bin.
            *** The best trigger is the one satisfying the kinematic cuts of the bin with the highest luminosity.
            **/
            for (int pbin = 0; pbin < numpbins; pbin++) {
                for (int etabin = 0; etabin < numetabins; etabin++) {
                    double maxTrigLumi = 0;
                    double thisTrigLumi;
                    Trigger* bestTrigger = NULL;
                    for (Trigger* trig : triggerSubList) {
                        thisTrigLumi = totalLumiVec[thisRNIndex + (trig->index + (pbin + etabin*numpbins)*numtrigs)*numruns];
                        if (maxTrigLumi < thisTrigLumi) {
                            maxTrigLumi = thisTrigLumi;
                            bestTrigger = trig;
                        } else if (maxTrigLumi == thisTrigLumi && bestTrigger != NULL && bestTrigger->min_pt < trig->min_pt) {
                            bestTrigger = trig;
                        }
                    }
                    kinematicLumiVec[pbin + etabin*numpbins] += maxTrigLumi;
                    if (rn == thisRunNumber) {
                        kinematicTriggerVec[pbin + etabin*numpbins] = bestTrigger;
                    }
                }
            }
        }
        /**** End loop over root files ****/


        /**** Instantiate trigger efficiency information ****/
        if (initTriggerEfficiencies) {
           
            if (debugStatements) cout << "Status: In triggerUtil.C (269): Gathering trigger efficiency information..." << endl; 

            double* triggerEfficiencies = getTriggerEfficiencies();
            for (int pbin = 0; pbin < numpbins; pbin++) {
                for (int etabin = 0; etabin < numetabins; etabin++) {
                    Trigger* bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins];
                    if (bestTrigger == NULL) continue; // for some reason no trigger covered this kinematic bin, so move on

                    kinematicEfficiencyVec[pbin + etabin*numpbins] = triggerEfficiencies[bestTrigger->index + pbin*numtrigs];
                }
            }
        }
        /**** End instantiate trigger efficiency information ****/


        /**** Debug statements ****/
        if (debugStatements) cout << "Status: In triggerUtil.C (462): Triggering scheme completed for run number " << runNumber << endl;  
        if (runNumber == 313063 && debugStatements) {
            cout << "Status: In triggerUtil.C (464): Example trigger assignment (run " << runNumber << "):" << endl << endl;
            for (int etabin = 0; etabin < numetabins; etabin++) {
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    if (kinematicTriggerVec[pbin + etabin*numpbins] != NULL) cout << Form("etabin=%i,\tpbin=%i, \teta=(%.1f, %.1f),\tp=(%i, %i),\ttrig=%s,\tlumi=%f", etabin, pbin, etabins[etabin], etabins[etabin+1], (int)pbins[pbin], (int)pbins[pbin+1], kinematicTriggerVec[pbin + etabin*numpbins]->name.c_str(), kinematicLumiVec[pbin + etabin*numpbins]) << endl;
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

