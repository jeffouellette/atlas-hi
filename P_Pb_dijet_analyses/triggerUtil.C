#include <vector>

using namespace std;

//=========================================================================================================
/** Variable declarations **/

// Global parameters
int runNumber;
const double dijet_pt_ratio_cutoff = 0.7; // Minimum subleading-to-leading jet ratio for the event to be considered a dijet.

const int global_max_pt = 6000; // Maximum allowed transverse momentum
const double min_eta = -4.9; // Minimum detectable pseudorapidity in hadronic calorimeter
const double max_eta = 4.9; // Maximum detectable pseudorapidity
const int trigthres = 10; // Jet pt threshold for triggers

// Directory information
const string workPath = "/Users/jeffouellette/Research/atlas-hi/P_Pb_dijet_analyses/";
const string rootPath = workPath + "rootFiles/";
string dataPath;
string trigPath;
string plotPath;
const int orig_run_list[27] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314112, 314157, 314170};
const int updated_run_list[14] = {313063, 313107, 313136, 313259, 313603, 313630, 313688, 313695, 313878, 313929, 314014, 314105, 314157, 314170};

// Useful constants
const float Z = 82;   // value of Z for Pb
const float A = 208;  // value of A for Pb
const float sqrt_s_nn = 8160; // Collision energy in CoM frame (GeV)

// Transverse momentum and pseudorapidity binning
const double pbins[42] = {25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 220., 240., 260., 280., 300., 350., 400., 500., 600., 800., 1100., 1500., 2000., 2500., 6000.};
const int numpbins = sizeof(pbins)/sizeof(pbins[0]) - 1;
const double etabins[9] = {-4.9, -3.2, -2., -1., 0, 1, 2., 3.2, 4.9};
//const double etabins[2] = {-4.9, 4.9}; // Used for avoiding eta-binning.
const int numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;
int numtrigs; // Total number of triggers

const bool runOldData = false;

const bool runPeriodA = true;
const bool runPeriodB = false;
bool periodA;

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
    if (!periodA) return prefactor * (jpt0 + TMath::Exp(-jeta0) + jpt1 * TMath::Exp(-jeta1));
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

    int min_pt;
    int max_pt;
    double lower_eta;    
    double upper_eta;
    int index;
    bool enabled;
    bool iontrigger;

    Trigger(string, int, double, double, bool);
    Trigger(const Trigger* t);

};

/**
 * Creates a Trigger object. By default, the maximum momentum and branching index are both 0. It is
 * expected that these values will be nonzero by the time the object is used purposefully.
 */
Trigger::Trigger(string thisname, int thismin_pt, double etal, double etau, bool thisiontrigger) {
    name = thisname;
    min_pt = thismin_pt;
    lower_eta = etal;
    upper_eta = etau;
    max_pt = 0;
    index = 0;
    enabled = true;
    iontrigger = thisiontrigger;
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
    enabled = t->enabled;
    iontrigger = t->iontrigger;
}

/**
 * Determines whether analysis of this run should be skipped.
 */
static bool skipRun (int rn) {
    if (rn < 313500 && !runPeriodA) return true;
    if (rn > 313500 && !runPeriodB) return true;
    if (runOldData) {
        bool contains_rn = false;
        int i = 0;
        while (i < sizeof(orig_run_list)/sizeof(int) && !contains_rn) {
            contains_rn = orig_run_list[i] == rn;
            i++;
        }
        return !contains_rn;

    } else {
        bool contains_rn = false;
        int i = 0;
        while (i < sizeof(updated_run_list)/sizeof(int) && !contains_rn) {
            contains_rn = updated_run_list[i] == rn;
            i++;
        }
        return !contains_rn;
    }
}

/**
 * Returns a pointer to an std:vector<int> of the run numbers currently being processed.
 */ 
static std::vector<int>* getRunNumbers() {
    std::vector<int>* rns = new std::vector<int>(0);
    if (runOldData) {
        for (int i = 0; i < sizeof(orig_run_list)/sizeof(int); i++) {
            rns->push_back(orig_run_list[i]);
        }
    } else {
        for (int i = 0; i < sizeof(updated_run_list)/sizeof(int); i++) {
            rns->push_back(updated_run_list[i]);
        }
    }
    return rns;
}


/** Trigger vectors used for pt and eta binning. **/
std::vector<Trigger*> trigger_vec(0);
int best_bins[numpbins*numetabins];
double kinematic_lumi_vec[numpbins*numetabins];

/**
 * Initializes triggers complete with momentum and pseudorapidity cutoffs.
 * To add new triggers, simply create a line like the one below with the trigger name, the momentum threshold, and its pseudorapidity interval.
 */
void initialize (int rn=0, bool initTriggerMaps = true) {

    cout << Form("Initializing trigger system for run number %i", rn) << endl;

    /** Store run number as a global variable **/
    runNumber = rn;
    if (rn < 313500) periodA = true;
    else periodA = false;

    /** Create an array of triggers **/

    // TODO use ION triggers for period A, but use NON-ION triggers for period B.
    // Ion triggers
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j10_ion_p320eta490_L1MBTS_1_1", 10+trigthres, 3.2, 4.9, true)); // SECONDARY
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j15_ion_p320eta490_L1MBTS_1_1", 15+trigthres, 3.2, 4.9, true));
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j25_ion_p320eta490_L1TE5", 25+trigthres, 3.2, 4.9, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j25_ion_p320eta490_L1TE10", 25+trigthres, 3.2, 4.9, true)); // SECONDARY
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j30_ion_0eta490_L1TE10", 30+trigthres, -4.9, 4.9, true));
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j30_ion_L1J5", 30+trigthres, min_eta, max_eta, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j35_ion_p320eta490_L1TE10", 35+trigthres, 3.2, 4.9, true)); // SECONDARY
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j40_ion_L1J5", 40+trigthres, min_eta, max_eta, true));
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j40_ion_L1J10", 40+trigthres, min_eta, max_eta, true)); // SECONDARY
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j45_ion_p200eta320", 45+trigthres, 2, 3.2, true));
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j45_ion_n200eta320", 45+trigthres, -3.2, -2, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j45_ion_p320eta490", 45+trigthres, 3.2, 4.9, true)); // SECONDARY
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j50_ion_L1J10", 50+trigthres, min_eta, max_eta, true));
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j55_ion_p200eta320", 55+trigthres, 2, 3.2, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j55_ion_n200eta320", 55+trigthres, -3.2, -2, true)); // SECONDARY
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j55_ion_p320eta490", 55+trigthres, 3.2, 4.9, true));
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j60_ion_L1J15", 60+trigthres, min_eta, max_eta, true)); // SECONDARY
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j60_ion_L1J20", 60+trigthres, min_eta, max_eta, true));
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j65_ion_p200eta320", 65+trigthres, 2, 3.2, true));
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j65_ion_n200eta320", 65+trigthres, -3.2, -2, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j65_ion_p320eta490", 65+trigthres, 3.2, 4.9, true)); // SECONDARY
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j75_ion_L1J20", 75+trigthres, min_eta, max_eta, true));
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j75_ion_p200eta320", 75+trigthres, 2, 3.2, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j75_ion_n200eta320", 75+trigthres, -3.2, -2, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j75_ion_p320eta490", 75+trigthres, 3.2, 4.9, true)); // SECONDARY
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j100_ion_L1J20", 100+trigthres, min_eta, max_eta, true));
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j110_ion_L1J30", 110+trigthres, min_eta, max_eta, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j125_ion_L1J30", 125+trigthres, min_eta, max_eta, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j150_ion_L1J30", 150+trigthres, min_eta, max_eta, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j175_ion_L1J50", 175+trigthres, min_eta, max_eta, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j200_ion_L1J50", 200+trigthres, min_eta, max_eta, true)); // SECONDARY
    if (!runOldData) trigger_vec.push_back(new Trigger("HLT_j250_ion_L1J50", 250+trigthres, min_eta, max_eta, true)); // SECONDARY

    // Alternate p-p triggers 
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j15_p320eta490_L1MBTS_1_1", 15+trigthres, 3.2, 4.9, false));
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j30_0eta490_L1TE10", 30+trigthres, -4.9, 4.9, false));
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j40_L1J5", 40+trigthres, min_eta, max_eta, false));
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j45_p200eta320", 45+trigthres, 2, 3.2, false));
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j50_L1J10", 50+trigthres, min_eta, max_eta, false));
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j55_p320eta490", 55+trigthres, 3.2, 4.9, false));
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j60", 60+trigthres, min_eta, max_eta, false));
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j65_p200eta320", 65+trigthres, 2, 3.2, false));
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j75_L1J20", 75+trigthres, min_eta, max_eta, false));
    if (runOldData) trigger_vec.push_back(new Trigger("HLT_j100_L1J20", 100+trigthres, min_eta, max_eta, false));

    numtrigs = trigger_vec.size();

    if (runOldData) {
        dataPath = workPath + "rundata/original/";
        plotPath = workPath + "Plots/original/";
        trigPath = workPath + "rootFiles/pt_data/trig_bin/original/"; 
    }
    else {
        dataPath = workPath + "rundata/updated/";
        plotPath = workPath + "Plots/updated/";
        trigPath = workPath + "rootFiles/pt_data/trig_bin/updated/";
    }
    
    /** Assign indices for tree branching. **/

    for (int i = 0; i < numtrigs; i++) {
        trigger_vec[i]->index = i;
    }


    /** Instantiate the pseudorapidity interval trigger vectors if required. **/
    
    if (initTriggerMaps) {

        cout << "Starting loop over triggers..." << endl;
    
        // Find all trigger analysis files.
        TSystemDirectory dir(trigPath.c_str(), trigPath.c_str());
        TList *files = dir.GetListOfFiles();
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
        double total_lumi_vec[numtrigs*numpbins*numetabins];
        int numtrigfires[numpbins * numetabins] = {};
        TVectorD* numTrigFirings;
        for (int n = 0; n < numpbins*numetabins*numtrigs; n++) {
            total_lumi_vec[n] = 0;
            if (n < numpbins*numetabins) kinematic_lumi_vec[n] = 0;
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
                if (trig->iontrigger == thisRunNumber > 313500) continue; // Use the correct triggers for this run

                // Integrate the number of times the trigger fired over eta and pt. If the result is 0, assume that the trigger was effectively inactive, so don't add the luminosity to that bin.
                double integral_deta = 0;
                for (int ebin = 0; ebin < numetabins; ebin++) {
                    thishist = (TH1D*)thisfile->Get(Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, index, ebin));
                    integral_deta += thishist->Integral();
                }
    
                for (int ebin = 0; ebin < numetabins; ebin++) {
                    for (int pbin = 0; pbin < numpbins; pbin++) {
                        if (integral_deta > 0) total_lumi_vec[index + (pbin + ebin*numpbins)*numtrigs] = thisLuminosity;
                        else total_lumi_vec[index + (pbin + ebin*numpbins)*numtrigs] = 0;
                    }
                }
            }
            thisfile->Close();
            // Calculate the best trigger to use for each bin, and be sure to scale by the correct deta, number of events, and luminosity.
            for (int pbin = 0; pbin < numpbins; pbin++) {
                for (int ebin = 0; ebin < numetabins; ebin++) {
                    double maxtrigfirings = 0;
                    int best_bin_allruns = 0;
                    for (Trigger* trig : trigger_vec) {
                        if (pbins[pbin] < trig->min_pt || (trig->iontrigger == thisRunNumber > 313500)) continue;
                        int index = trig->index;
                        if ((*numTrigFirings)[index + (pbin + ebin*numpbins)*numtrigs] > maxtrigfirings) {
                            maxtrigfirings = (*numTrigFirings)[index + (pbin + ebin*numpbins)*numtrigs];
                            best_bin_allruns = index;
                            if (rn == thisRunNumber) best_bins[pbin + ebin*numpbins] = index;
                        }
                    }
                    if (rn == thisRunNumber) numtrigfires[pbin + ebin*numpbins] = maxtrigfirings;
                    kinematic_lumi_vec[pbin + ebin*numpbins] += total_lumi_vec[best_bin_allruns + (pbin + ebin*numpbins)*numtrigs];
                }
            }
            rnIndex++;
        }

        if (runNumber == 313063) {
            cout << Form("Example trigger assignment (run %i):", runNumber) << endl << endl;
            for (int ebin = 0; ebin < numetabins; ebin++) {
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    cout << Form("ebin=%i,\tpbin=%i, \teta=(%.1f, %.1f),\tp=(%i, %i),\ttrig=%i,\tcounts=%i", ebin, pbin, etabins[ebin], etabins[ebin+1], (int)pbins[pbin], (int)pbins[pbin+1], best_bins[pbin + ebin*numpbins], numtrigfires[pbin + ebin*numpbins]) << endl;
                }
            } 
        }
    }

    cout << Form("Initialization complete for run number %i", runNumber) << endl;  
    return;
}

