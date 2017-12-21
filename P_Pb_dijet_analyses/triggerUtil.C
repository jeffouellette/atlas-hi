#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

//=========================================================================================================
/** Variable declarations **/

// ANALYSIS DATA SELECTION - SET THESE VARIABLES FOR DESIRED DATA SELECTION CHOICE
const int useDataVersion = 4;
const bool runPeriodA = false;
const bool runPeriodB = true;
const bool debugStatements = false;
const bool considerDisabledTriggers = false;

// Transverse momentum and pseudorapidity binning
const double pbins[66] = {25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145., 150., 155., 160., 165., 170., 175., 180., 185., 190., 195., 200., 205., 210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 320., 340., 360., 380., 400., 425., 450., 475., 500., 550., 600., 700., 800., 1000., 1250., 1500., 2000., 2500., 6000.};
const int numpbins = sizeof(pbins)/sizeof(pbins[0]) - 1;
const double etabins[9] = {-4.9, -3.2, -2., -1., 0, 1, 2., 3.2, 4.9};
//const double etabins[2] = {-4.9, 4.9}; // Used for avoiding eta-binning.
const int numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;
int numtrigs; // Total number of triggers

// Global parameters
int runNumber;
bool periodA;
bool useIonTrigs;
const double dijet_pt_ratio_cutoff = 0.7; // Minimum subleading-to-leading jet ratio for the event to be considered a dijet.

const int global_max_pt = 6000; // Maximum allowed transverse momentum
const double min_eta = -4.9; // Minimum detectable pseudorapidity in hadronic calorimeter
const double max_eta = 4.9; // Maximum detectable pseudorapidity
//const int trigthres = 10; // Jet pt threshold for triggers
const int trigthres[numetabins] = {10, 10, 10, 10, 10, 10, 10, 10};

// Directory information
const string workPath = "/Users/jeffouellette/Research/atlas-hi/P_Pb_dijet_analyses/";
string rootPath = workPath + "rootFiles/";
string dataPath = workPath + "rundata/";
string plotPath = workPath + "Plots/";
string ptPath = workPath + "rootFiles/pt_data/";
string trigPath = workPath + "rootFiles/trig_data/";
const int run_list_v1[27] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314112, 314157, 314170};
const int run_list_v2[14] = {313063, 313107, 313136, 313259, 313603, 313630, 313688, 313695, 313878, 313929, 314014, 314105, 314157, 314170};
const int run_list_v3[30] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313688, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170};

// Useful constants
const float Z = 82;   // value of Z for Pb
const float A = 208;  // value of A for Pb
const float sqrt_s_nn = 8160; // Collision energy in CoM frame (GeV)


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
 * Trigger stores information about a trigger, including a jet momenta range, pseudorapidity interval,
 * name, and a (unique) branching index for event analysis.
 */
class Trigger {
    
    public:
    string name;

    int* min_pt;
    //int min_pt;
    int max_pt;
    double lower_eta;    
    double upper_eta;
    int index;
    bool enabled;
    bool iontrigger;
    bool disabled;

    Trigger(string, int, double, double, bool);
    Trigger(string, int, double, double, bool, bool);
    Trigger(const Trigger* t);

};

/**
 * Creates a Trigger object. By default, the maximum momentum and branching index are both 0. It is
 * expected that these values will be nonzero by the time the object is used purposefully.
 */
Trigger::Trigger(string thisname, int thismin_pt, double etal, double etau, bool thisiontrigger) {
    name = thisname;
    //min_pt = thismin_pt + trigthres;
    min_pt = new int[numetabins];
    for (int ebin = 0; ebin < numetabins; ebin++) {
        min_pt[ebin] = thismin_pt + trigthres[ebin];
    }
    lower_eta = etal;
    upper_eta = etau;
    max_pt = 0;
    index = 0;
    iontrigger = thisiontrigger;
    disabled = false;
}

/**
 * Creates a Trigger object. By default, the maximum momentum and branching index are both 0. It is
 * expected that these values will be nonzero by the time the object is used purposefully.
 */
Trigger::Trigger(string thisname, int thismin_pt, double etal, double etau, bool thisiontrigger, bool thisdisabled) {
    name = thisname;
    //min_pt = thismin_pt + trigthres;
    min_pt = new int[numetabins];
    for (int ebin = 0; ebin < numetabins; ebin++) {
        min_pt[ebin] = thismin_pt + trigthres[ebin];
    }
    lower_eta = etal;
    upper_eta = etau;
    max_pt = 0;
    index = 0;
    iontrigger = thisiontrigger;
    disabled = thisdisabled && !considerDisabledTriggers;
}

/**
 * Creates a copy of trigger t.
 */
Trigger::Trigger(const Trigger* t) {
    name = t->name;
    min_pt = t->min_pt;
    max_pt = t->max_pt;
    lower_eta = t->lower_eta;
    upper_eta = t->upper_eta;
    index = t->index;
    iontrigger = t->iontrigger;
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
        case 3: {
            while (i < sizeof(run_list_v3)/sizeof(int) && !contains_rn) {
                contains_rn = run_list_v3[i] == rn;
                i++;
            }
            return !contains_rn;
        }
        case 2: {
            while (i < sizeof(run_list_v2)/sizeof(int) && !contains_rn) {
                contains_rn = run_list_v2[i] == rn;
                i++;
            }
            return !contains_rn;
        }
        default: {
            while (i < sizeof(run_list_v1)/sizeof(int) && !contains_rn) {
                contains_rn = run_list_v1[i] == rn;
                i++;
            }
            return !contains_rn;
        } 
    }
}

/**
 * Returns a pointer to an std:vector<int> of the run numbers currently being processed.
 */ 
static std::vector<int>* getRunNumbers() {
    std::vector<int>* rns = new std::vector<int>(0);
    switch (useDataVersion) {
        case 3: {
            for (int i = 0; i < sizeof(run_list_v3)/sizeof(int); i++) {
                rns->push_back(run_list_v3[i]);
            }
            break;
        }
        case 2: {
            for (int i = 0; i < sizeof(run_list_v2)/sizeof(int); i++) {
                rns->push_back(run_list_v2[i]);
            }
            break;
        }
        default: {
            for (int i = 0; i < sizeof(run_list_v1)/sizeof(int); i++) {
                rns->push_back(run_list_v1[i]);
            }
            break;
        }
    }
    return rns;
}


/** Trigger vectors used for pt and eta binning. **/
std::vector<Trigger*> trigger_vec(0);
int best_bins[numpbins*numetabins];
double kinematic_lumi_vec[numpbins*numetabins];


/**
 * Adds triggers relevant to period A analysis.
 */
void add_period_A_triggers() {
    // Use ION triggers for period A, but use NON-ION triggers for period B.
    // Ion triggers
    bool useLatestTriggers = useDataVersion >= 2;

    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j10_ion_p320eta490_L1MBTS_1_1", 10, 3.2, 4.9, true));
                           trigger_vec.push_back(new Trigger("HLT_j15_ion_p320eta490_L1MBTS_1_1", 15, 3.2, 4.9, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j25_ion_p320eta490_L1TE5", 25, 3.2, 4.9, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j25_ion_p320eta490_L1TE10", 25, 3.2, 4.9, true));
                           trigger_vec.push_back(new Trigger("HLT_j30_ion_0eta490_L1TE10", 30, -4.9, 4.9, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j30_ion_L1J5", 30, -3.2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j35_ion_p320eta490_L1TE10", 35, 3.2, 4.9, true));
                           trigger_vec.push_back(new Trigger("HLT_j40_ion_L1J5", 40, -3.2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j40_ion_L1J10", 40, -3.2, 3.2, true));
                           trigger_vec.push_back(new Trigger("HLT_j45_ion_p200eta320", 45, 2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j45_ion_n200eta320", 45, -3.2, -2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j45_ion_p320eta490", 45, 3.2, 4.9, true));
                           trigger_vec.push_back(new Trigger("HLT_j50_ion_L1J10", 50, -3.2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j55_ion_p200eta320", 55, 2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j55_ion_n200eta320", 55, -3.2, -2, true));
                           trigger_vec.push_back(new Trigger("HLT_j55_ion_p320eta490", 55, 3.2, 4.9, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j60_ion_L1J15", 60, -3.2, 3.2, true));
                           trigger_vec.push_back(new Trigger("HLT_j60_ion_L1J20", 60, -3.2, 3.2, true));
                           trigger_vec.push_back(new Trigger("HLT_j65_ion_p200eta320", 65, 2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j65_ion_n200eta320", 65, -3.2, -2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j65_ion_p320eta490", 65, 3.2, 4.9, true));
                           trigger_vec.push_back(new Trigger("HLT_j75_ion_L1J20", 75, -3.2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j75_ion_p200eta320", 75, 2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j75_ion_n200eta320", 75, -3.2, -2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j75_ion_p320eta490", 75, 3.2, 4.9, true));
                           trigger_vec.push_back(new Trigger("HLT_j100_ion_L1J20", 100, -3.2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j110_ion_L1J30", 110, -3.2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j125_ion_L1J30", 125, -3.2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j150_ion_L1J30", 150, -3.2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j175_ion_L1J50", 175, -3.2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j200_ion_L1J50", 200, -3.2, 3.2, true));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j250_ion_L1J50", 250, -3.2, 3.2, true));
    return;
}


/**
 * Adds triggers relevant to period B analysis.
 */
void add_period_B_triggers() {
    // Alternate p-p triggers (non-ion)
    bool useLatestTriggers = useDataVersion >= 3;
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j10_p320eta490_L1MBTS_1_1", 10, 3.2, 4.9, false));
                           trigger_vec.push_back(new Trigger("HLT_j15_p320eta490_L1MBTS_1_1", 15, 3.2, 4.9, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j25_p320eta490_L1TE5", 25, 3.2, 4.9, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j25_p320eta490_L1TE10", 25, 3.2, 4.9, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j30_L1J5", 30, -3.2, 3.2, false));
                           trigger_vec.push_back(new Trigger("HLT_j30_0eta490_L1TE10", 30, -4.9, 4.9, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j35_p320eta490_L1TE10", 35, 3.2, 4.9, false));
                           trigger_vec.push_back(new Trigger("HLT_j40_L1J5", 40, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j40_L1J10", 40, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j45_n200eta320", 45, -3.2, -2, false));
                           trigger_vec.push_back(new Trigger("HLT_j45_p200eta320", 45, 2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j45_p320eta490", 45, 3.2, 4.9, false));
                           trigger_vec.push_back(new Trigger("HLT_j50_L1J10", 50, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j55_n200eta320", 55, -3.2, -2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j55_p200eta320", 55, 2, 3.2, false));
                           trigger_vec.push_back(new Trigger("HLT_j55_p320eta490", 55, 3.2, 4.9, false));
                           trigger_vec.push_back(new Trigger("HLT_j60", 60, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j60_L1J15", 60, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j65_n200eta320", 65, -3.2, -2, false));
                           trigger_vec.push_back(new Trigger("HLT_j65_p200eta320", 65, 2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j65_p320eta490", 65, 3.2, 4.9, false));
                           trigger_vec.push_back(new Trigger("HLT_j75_L1J20", 75, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j75_n200eta320", 75, -3.2, -2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j75_p200eta320", 75, 2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j75_p320eta490", 75, 3.2, 4.9, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j90_L1J20", 90, -3.2, 3.2, false));
                           trigger_vec.push_back(new Trigger("HLT_j100_L1J20", 100, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j110", 110, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j125_L1J30", 125, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j150_L1J30", 150, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j175", 175, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j200", 200, -3.2, 3.2, false));
    if (useLatestTriggers) trigger_vec.push_back(new Trigger("HLT_j250", 250, -3.2, 3.2, false));
    return;
}


/**
 * Initializes triggers complete with momentum and pseudorapidity cutoffs.
 * To add new triggers, simply create a line like the one below with the trigger name, the momentum threshold, and its pseudorapidity interval.
 */
void initialize (int rn=0, bool initTriggerMaps=true, bool skip_irrelevant_triggers=false) {

    if (debugStatements) cout << Form("Initializing trigger system for run %i...", rn);

    /** Store run number as a global variable **/
    runNumber = rn;
    if (rn < 313500) periodA = true;
    else periodA = false;
    if (rn < 313629) useIonTrigs = true;
    else useIonTrigs = false;

    /** Create an array of triggers **/
    if (useDataVersion <= 3) {
        if (!skip_irrelevant_triggers || periodA) add_period_A_triggers();
        if (!skip_irrelevant_triggers || !periodA) add_period_B_triggers();
    } // end v1, v2, v3 trigger generation
    else {
        ifstream triggerListFile ((workPath + "triggerList.txt").c_str());
        if (!triggerListFile.is_open()) {
            cout << "\rtriggerUtil::initialize: Cannot find triggerList.txt!" << endl;
            throw runtime_error("ifstream file not found");
        }

        string currLine;
        int trigBlockLowerRunNumber = 0;
        int trigBlockUpperRunNumber = 0;
        while (getline(triggerListFile, currLine)) { // get each consecutive line
            //parse currline
            if (currLine[0] == 'R') {
                trigBlockLowerRunNumber = trigBlockUpperRunNumber;
                stringstream trigBlockToInt (currLine.substr(4));
                trigBlockToInt >> trigBlockUpperRunNumber;
                continue;
            }
            int space1 = currLine.find(" ", 0);
            int space2 = currLine.find(" ", space1 + 1);
            int space3 = currLine.find(" ", space2 + 1);

            string trigName = currLine.substr(0, space1);

            stringstream trigPtToFloat (currLine.substr(space1 + 1, space2-space1-1));
            float trigPtFloat;
            trigPtToFloat >> trigPtFloat;

            stringstream trigLetaToFloat (currLine.substr(space2 + 1, space3-space2-1));
            float trigLetaFloat;
            trigLetaToFloat >> trigLetaFloat;

            stringstream trigUetaToFloat (currLine.substr(space3 + 1, currLine.length()-space3-1));
            float trigUetaFloat;
            trigUetaToFloat >> trigUetaFloat;

            bool useThisTrigger = (trigBlockLowerRunNumber <= runNumber && runNumber < trigBlockUpperRunNumber) || (runNumber == 0);
            if (!skip_irrelevant_triggers || useThisTrigger) {
                trigger_vec.push_back(new Trigger(trigName, trigPtFloat, trigLetaFloat, trigUetaFloat, trigName.find("ion")!=string::npos, useThisTrigger));
            }
        }
        triggerListFile.close();
    } // end v4+ trigger generation

    cout << Form("\rTriggers initialized for run %i          ", rn) << endl;
    numtrigs = trigger_vec.size();

    if (debugStatements) cout << "Num trigs = " << numtrigs << endl;
    
    string versionString = "v" + to_string(useDataVersion) + "/";
    rootPath = rootPath + versionString;
    dataPath = dataPath + versionString;
    plotPath = plotPath + versionString;
    ptPath = rootPath + "pt_data/";
    trigPath = rootPath + "trig_data/";

    if (debugStatements) cout << "Trigger pt bin path is " << trigPath << endl;    
    if (debugStatements) cout << "Saving plots to " << plotPath << endl;

    /** Assign indices for tree branching. **/

    for (int i = 0; i < numtrigs; i++) {
        trigger_vec[i]->index = i;
    }


    /** Instantiate the pseudorapidity interval trigger vectors if required. **/
    
    if (initTriggerMaps) {

        if (debugStatements) cout << "Starting loop over triggers...";
    
        // Find all trigger analysis files.
        TSystemDirectory dir(trigPath.c_str(), trigPath.c_str());
        TList* files = dir.GetListOfFiles();
        std::vector<TString> filenames;
        if (files) {
            TSystemFile *file;
            TString fname;
            TIter next(files);
            while ((file=(TSystemFile*)next())) {
                fname = file->GetName();
                if (!file->IsDirectory() && fname.EndsWith(".root")) {
                    filenames.push_back(fname);
                    if (rn == 313063) cout << "Found " << fname.Data() << endl;
                }
            }
        }

        // First combine trigger data from all runs into one histogram for each trigger. If the trigger never fired in a run, assume it wasn't on so don't add its luminosity.
        TH1D* thishist;
        int rnIndex = 0;
        const int numruns = filenames.size();
        double* total_lumi_vec = new double[numtrigs*numpbins*numetabins];
        int* numtrigfires = new int[numpbins * numetabins];
        int* best_bin_allruns_vec = new int[numtrigs * numpbins * numetabins];
        TVectorD* numTrigFirings;
        for (int n = 0; n < numpbins*numetabins*numtrigs; n++) {
            total_lumi_vec[n] = 0;
            if (n < numpbins*numetabins) kinematic_lumi_vec[n] = 0;
            if (n < numpbins*numetabins) numtrigfires[n] = 0;
            if (n < numtrigs*numpbins*numetabins) best_bin_allruns_vec[n] = 0;
        }

        for (TString filename : filenames) {
            TFile* thisfile = new TFile(Form("%s%s", trigPath.c_str(), filename.Data()), "READ");
    
            int thisRunNumber = (int)(*((TVectorD*)thisfile->Get("run_vec")))[0];
            if (skipRun(thisRunNumber)) { // Only allow desired runs to be considered
                thisfile->Close();
                continue;
            }
    
            double thisLuminosity = (*((TVectorD*)thisfile->Get("lum_vec")))[0];
            
            assert (numetabins == (int)(*((TVectorD*)thisfile->Get("run_vec")))[1]); // Quick assertions that the trigger data matches what we want.
            assert (numtrigs == (int)(*((TVectorD*)thisfile->Get("run_vec")))[2]);

            numTrigFirings = (TVectorD*)thisfile->Get("trig_fire_vec");
            
            for (Trigger* trig : trigger_vec) {
                int index = trig->index;
                if (trig->iontrigger == thisRunNumber > 313603) continue; // Use the correct triggers for this run

                // Integrate the number of times the trigger fired over eta and pt. If the result is 0, assume that the trigger was effectively inactive, so don't add the luminosity to that bin.
                double integral_deta = 0;
                for (int ebin = 0; ebin < numetabins; ebin++) {
                    thishist = (TH1D*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, index, ebin));
                    integral_deta += thishist->Integral();
                }
    
                for (int ebin = 0; ebin < numetabins; ebin++) {
                    for (int pbin = 0; pbin < numpbins; pbin++) {
                        if (integral_deta > 0 && trig->min_pt[ebin] <= pbins[pbin] && trig->lower_eta <= etabins[ebin] && etabins[ebin+1] <= trig->upper_eta) total_lumi_vec[index + (pbin + ebin*numpbins)*numtrigs] = thisLuminosity;
                        else total_lumi_vec[index + (pbin + ebin*numpbins)*numtrigs] = 0;
                    }
                }
            }
            thisfile->Close();
            delete thisfile;
            // Calculate the best trigger to use for each bin, and be sure to scale by the correct deta, number of events, and luminosity.
            for (int pbin = 0; pbin < numpbins; pbin++) {
                for (int ebin = 0; ebin < numetabins; ebin++) {
                    double maxtrigfirings = 0;
                    int best_bin_allruns = 0;
                    int last_best_bin = 0;
                    for (Trigger* trig : trigger_vec) {
                        if (pbins[pbin] < trig->min_pt[ebin] || (trig->iontrigger == thisRunNumber > 313500) || trig->disabled) continue;
                        int index = trig->index;
                        if ((*numTrigFirings)[index + (pbin + ebin*numpbins)*numtrigs] > maxtrigfirings) {
                            maxtrigfirings = (*numTrigFirings)[index + (pbin + ebin*numpbins)*numtrigs];
                            last_best_bin = best_bin_allruns;
                            best_bin_allruns = index;
                            if (rn == thisRunNumber) best_bins[pbin + ebin*numpbins] = index;
                        }
                    }
                    if (rn == thisRunNumber) numtrigfires[pbin + ebin*numpbins] = maxtrigfirings; 
//                    best_bin_allruns_vec[rnIndex + (pbin + ebin*numpbins)*numruns] = last_best_bin;
//                    kinematic_lumi_vec[pbin + ebin*numpbins] += total_lumi_vec[last_best_bin + (pbin + ebin*numpbins)*numtrigs];
                    best_bin_allruns_vec[rnIndex + (pbin + ebin*numpbins)*numruns] = best_bin_allruns;
                    kinematic_lumi_vec[pbin + ebin*numpbins] += total_lumi_vec[best_bin_allruns + (pbin + ebin*numpbins)*numtrigs];
                }
            }
            rnIndex++;
        }
        delete [] total_lumi_vec;
        delete [] best_bin_allruns_vec;
        delete [] numtrigfires;
/*        rnIndex = 0;
        for (TString filename : filenames) {
            TFile* thisfile = new TFile(Form("%s%s", trigPath.c_str(), filename.Data()), "READ");
    
            int thisRunNumber = (int)(*((TVectorD*)thisfile->Get("run_vec")))[0];
            if (skipRun(thisRunNumber)) { // Only allow desired runs to be considered
                thisfile->Close();
                continue;
            }
            for (int ebin = 0; ebin < numetabins; ebin++) {
                int act_ebin = ebin;
                if (thisRunNumber < 313500) act_ebin = numetabins - ebin - 1;
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    int best_bin_allruns = best_bin_allruns_vec[rnIndex + (pbin + act_ebin*numpbins)*numruns];
                    kinematic_lumi_vec[pbin + ebin*numpbins] += total_lumi_vec[best_bin_allruns + (pbin + act_ebin*numpbins)*numruns];
                }
            }
            rnIndex++; 
        }*/
        if (debugStatements) cout << Form("\rInitialization complete for run number %i", runNumber) << endl;  

        if (runNumber == 313063 && debugStatements) {
            cout << Form("Example trigger assignment (run %i):", runNumber) << endl << endl;
            for (int ebin = 0; ebin < numetabins; ebin++) {
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    cout << Form("ebin=%i,\tpbin=%i, \teta=(%.1f, %.1f),\tp=(%i, %i),\ttrig=%i,\tcounts=%i", ebin, pbin, etabins[ebin], etabins[ebin+1], (int)pbins[pbin], (int)pbins[pbin+1], best_bins[pbin + ebin*numpbins], numtrigfires[pbin + ebin*numpbins]) << endl;
                }
            } 
        }
    }
    return;
}

