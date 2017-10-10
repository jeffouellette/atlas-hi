#include <vector>

using namespace std;

//=========================================================================================================
/** Variable declarations **/

// Global parameters
int runNumber;
const double dijet_pt_frac_cutoff = 0.2; // Maximum fraction of transverse momentum of all jets but the
                                        // first two for the event to be considered a dijet.
const int global_max_pt = 6000; // Maximum allowed transverse momentum
const double min_eta = -4.9; // Minimum detectable pseudorapidity in hadronic calorimeter
const double max_eta = 4.9; // Maximum detectable pseudorapidity
const int trigthres = 10; // Jet pt threshold for triggers
int numtrigs; // Total number of triggers

// Directory information
const char* trig_dir = "/Users/jeffouellette/Research/atlas-hi/trig_data/";

// Useful constants
const float Z = 82;   // value of Z for Pb
const float A = 208;  // value of A for Pb
const float sqrt_s_nn = 8160; // Collision energy in CoM frame (GeV)

// Transverse momentum and pseudorapidity binning
const double pbins[42] = {25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 220., 240., 260., 280., 300., 350., 400., 500., 600., 800., 1100., 1500., 2000., 2500., 6000.};
const int numpbins = sizeof(pbins)/sizeof(pbins[0]) - 1;
const double etabins[9] = {-4.9, -3.2, -2., -1., 0, 1, 2., 3.2, 4.9};
const int numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;


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
        double* arr = linspace(TMath::Log2(hi/(2*num)), loghi, num);
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
    //if (!periodA) return prefactor * (jpt0 * TMath::Exp(jeta0) + jpt1 * TMath::Exp(jeta1));
    return prefactor * (jpt0 * TMath::Exp(-jeta0) + jpt1 * TMath::Exp(-jeta1));
}


/**
 * Returns xa for the event.
 */
double get_xa(double jpt0, double jpt1, double jeta0, double jeta1, bool periodA) {
    double prefactor = TMath::Sqrt(A/Z) / sqrt_s_nn;
    //if (!periodA) return prefactor * (jpt0 + TMath::Exp(-jeta0) + jpt1 * TMath::Exp(-jeta1));
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

    Trigger(string, int, double, double);
    Trigger(const Trigger* t);

};

/**
 * Creates a Trigger object. By default, the maximum momentum and branching index are both 0. It is
 * expected that these values will be nonzero by the time the object is used purposefully.
 */
Trigger::Trigger(string thisname, int thismin_pt, double etal, double etau) {
    name = thisname;
    min_pt = thismin_pt;
    lower_eta = etal;
    upper_eta = etau;
    max_pt = 0;
    index = 0;
    enabled = true;
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
}


/** Trigger vectors used for pt and eta binning. **/
std::vector<TString> trig_file_names(0);
std::vector<Trigger*> trigger_vec(0);
std::vector<std::vector<Trigger*>*>* trigger_pt_eta_bin_map;


/**
 * Initializes triggers complete with momentum and pseudorapidity cutoffs.
 * To add new triggers, simply create a line like the one below with the trigger name, the momentum threshold, and its pseudorapidity interval.
 */
void initialize (int rn, bool initTriggerMaps = true) {

    cout << Form("Initializing trigger system for run number %i", rn) << endl;

    /** Store run number as a global variable **/
    runNumber = rn;

    /** Create an array of triggers **/

    // Ion triggers
    trigger_vec.push_back(new Trigger("HLT_j15_ion_p320eta490_L1MBTS_1_1", 15+trigthres, 3.2, 4.9));
    trigger_vec.push_back(new Trigger("HLT_j30_ion_0eta490_L1TE10", 30+trigthres, -4.9, 4.9));
    trigger_vec.push_back(new Trigger("HLT_j40_ion_L1J5", 40+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j45_ion_p200eta320", 45+trigthres, 2, 3.2));
    trigger_vec.push_back(new Trigger("HLT_j50_ion_L1J10", 50+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j55_ion_p320eta490", 55+trigthres, 3.2, 4.9));
    trigger_vec.push_back(new Trigger("HLT_j60_ion_L1J20", 60+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j65_ion_p200eta320", 65+trigthres, 2, 3.2));
    trigger_vec.push_back(new Trigger("HLT_j75_ion_L1J20", 75+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j100_ion_L1J20", 100+trigthres, min_eta, max_eta));

    // Alternate p-p triggers 
/*  trigger_vec.push_back(new Trigger("HLT_j15_p320eta490_L1MBTS_1_1", 15+trigthres, 3.2, 4.9);
    trigger_vec.push_back(new Trigger("HLT_j30_0eta490_L1TE10", 30+trigthres, -4.9, 4.9));
    trigger_vec.push_back(new Trigger("HLT_j40_L1J5", 40+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j45_p200eta320", 45+trigthres, 2, 3.2));
    trigger_vec.push_back(new Trigger("HLT_j50_L1J10", 50+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j55_p320eta490", 55+trigthres, 3.2, 4.9));
    trigger_vec.push_back(new Trigger("HLT_j60", 60+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j65_p200eta320", 65+trigthres, 2, 3.2));
    trigger_vec.push_back(new Trigger("HLT_j75_L1J20", 75+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j100_L1J20", 100+trigthres, min_eta, max_eta));*/

    numtrigs = trigger_vec.size();
    
    /** Assign indices for tree branching. **/

    for (int i = 0; i < numtrigs; i++) {
        trigger_vec[i]->index = i;
    }

    Trigger* newtrig;
    Trigger* maxtrig;
    /** Instantiate the pseudorapidity interval trigger vectors if required. **/
    if (initTriggerMaps) {

        TFile* runFile = new TFile(Form("./trig_data/run_%i.root", runNumber), "READ");
        TH3D* h_meta = (TH3D*)runFile->Get("eta_pt_trig");

        trigger_pt_eta_bin_map = new std::vector<std::vector<Trigger*>*>(numetabins, new std::vector<Trigger*>(numpbins, 0));
        for (int ebin = 0; ebin < numetabins; ebin++) {
            for (int pbin = 0; pbin < numpbins; pbin++) {
                maxtrig = trigger_vec[0];
                for (Trigger* trig : trigger_vec) {
                    if (h_meta->GetBinContent((trig->index)+1, pbin+1, ebin+1) > h_meta->GetBinContent((maxtrig->index)+1, pbin+1, ebin+1)) {
                        maxtrig = trig;
                    }
                }
                //newtrig = new Trigger(*maxtrig);
                (*((*trigger_pt_eta_bin_map)[ebin]))[pbin] = new Trigger(*maxtrig);
                if (h_meta->GetBinContent((maxtrig->index)+1, pbin+1, ebin+1) == 0.) {
                    (*((*trigger_pt_eta_bin_map)[ebin]))[pbin]->enabled = false;
                }
                /*if (runNumber == 313259) {
                    cout << Form("ebin=%i, pbin=%i, trig=%s with value %f, enabled=%i", ebin, pbin, maxtrig->name.c_str(), h_meta->GetBinContent((maxtrig->index)+1, pbin+1, ebin+1), (*((*trigger_pt_eta_bin_map)[ebin]))[pbin]->enabled) << endl;
                }*/
            }
        }
        
        // Assign upper pt cuts. Triggers are already arranged properly by construction of the mapping.
        for (int ebin = 0; ebin < numetabins; ebin++) {
            for (int pbin = 0; pbin < numpbins-1; pbin++) {
                (*((*trigger_pt_eta_bin_map)[ebin]))[pbin]->max_pt = (*((*trigger_pt_eta_bin_map)[ebin]))[pbin+1]->min_pt;
            }
            (*((*trigger_pt_eta_bin_map)[ebin]))[numpbins-1]->max_pt = pbins[numpbins];
        }
        runFile->Close();

        if (rn == 313063) {
            cout << Form("Example trigger assignment (run %i):", runNumber) << endl << endl;
            for (int ebin = 0; ebin < numetabins; ebin++) {
                for (int pbin = 0; pbin < numpbins; pbin++) {
                    cout << Form("ebin=%i,\tpbin=%i, \teta=(%.1f, %.1f),\tp=(%i, %i),\ttrig=%s,\tenabled=%i", ebin, pbin, etabins[ebin], etabins[ebin+1], (int)pbins[pbin], (int)pbins[pbin+1], (*((*trigger_pt_eta_bin_map)[ebin]))[pbin]->name.c_str(), (*((*trigger_pt_eta_bin_map)[ebin]))[pbin]->enabled) << endl;
                }
            } 
        }

    }

    cout << Form("Initialization complete for run number %i", runNumber) << endl;  
    return;
}
