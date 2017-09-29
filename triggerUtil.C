#include <vector>

using namespace std;

//=========================================================================================================

// Global parameters

const float dijet_pt_frac_cutoff = 0.2; // Maximum fraction of transverse momentum of all jets but the
                                        // first two for the event to be considered a dijet.
const int global_max_pt = 2000; // Maximum allowed transverse momentum

const double min_eta = -1e3; // Minimum detectable pseudorapidity
const double max_eta = 1e3; // Maximum detectable pseudorapidity

const int trigthres = 10; // Jet pt threshold for triggers
int numtrigs; // Total number of triggers

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
    double loglo;
    if (lo > 0) loglo = TMath::Log2(lo);
    else loglo = TMath::Log2(0.01*(hi-lo));
    double loghi = TMath::Log2(hi);
    double* arr = linspace(loglo, loghi, num);
    for (int i = 0; i <= num; i++) {
        arr[i] = TMath::Power(2, arr[i]);
    }
    return arr;
}

/**
 * Returns xp for the event.
 */
double get_xp(double jpt0, double jpt1, double jeta0, double jeta1, bool periodA) {
    double prefactor = TMath::Sqrt(Z/A) / sqrt_s_nn;
    if (!periodA) return prefactor * (jpt0 * TMath::Exp(jeta0) + jpt1 * TMath::Exp(jeta1));
    else return prefactor * (jpt0 * TMath::Exp(-jeta0) + jpt1 * TMath::Exp(-jeta1));
}

/**
 * Returns xa for the event.
 */
double get_xa(double jpt0, double jpt1, double jeta0, double jeta1, bool periodA) {
    double prefactor = TMath::Sqrt(A/Z) / sqrt_s_nn;
    if (!periodA) return prefactor * (jpt0 + TMath::Exp(-jeta0) + jpt1 * TMath::Exp(-jeta1));
    else return prefactor * (jpt0 * TMath::Exp(jeta0) + jpt1 * TMath::Exp(jeta1));
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


//=========================================================================================================
// Store trigger names as an array of strings and loop over them to set the branch addresses --- DEPRECATED
/*static const int trigLength = 24;
const char* m_trig_string[trigLength] = {
    "HLT_j15_ion_p320eta490_L1MBTS_1_1",
    "HLT_j15_p320eta490_L1MBTS_1_1",
    "HLT_j30_ion_0eta490_L1TE10",
    "HLT_j30_0eta490_L1TE10",
    "HLT_j40_ion_L1J5",
    "HLT_j40_L1J5",
    "HLT_j45_ion_p200eta320",
    "HLT_j45_p200eta320",
    "HLT_j50_ion_L1J10",
    "HLT_j50_L1J10",
    "HLT_j55_ion_p320eta490",
    "HLT_j55_p320eta490",
    "HLT_j60_ion_L1J20",
    "HLT_j60",
    "HLT_j65_ion_p200eta320",
    "HLT_j65_p200eta320",
    "HLT_j75_ion_L1J20",
    "HLT_j75_L1J20",
    "HLT_j100_ion_L1J20",
    "HLT_j100_L1J20",
    "HLT_2j10_ion_p320eta490_L1TE10",
    "HLT_2j10_p320eta490_L1TE10",
    "HLT_2j30_ion_p320eta490",
    "HLT_2j30_p320eta490"
};*/
//=========================================================================================================


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

    Trigger(string, int, double, double);
    Trigger(const Trigger &t);

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
}

Trigger::Trigger(const Trigger &t) {
    name = t.name;
    min_pt = t.min_pt;
    max_pt = t.max_pt;
    lower_eta = t.lower_eta;
    upper_eta = t.upper_eta;
    index = t.index;
}


std::vector<Trigger*> trigger_vec(0);
std::vector<Trigger*> triggers_n490eta320(0);
std::vector<Trigger*> triggers_n320eta200(0);
std::vector<Trigger*> triggers_n200eta100(0);
std::vector<Trigger*> triggers_n100eta0(0);
std::vector<Trigger*> triggers_p0eta100(0);
std::vector<Trigger*> triggers_p100eta200(0);
std::vector<Trigger*> triggers_p200eta320(0);
std::vector<Trigger*> triggers_p320eta490(0);

std::vector<Trigger*> get_trigs_eta_range(double etal, double etau) {
    std::vector<Trigger*> triglist;

    // First grab all the triggers in the eta range specified
    for (Trigger* trig : trigger_vec) {
        if (trig->lower_eta <= etal && etau <= trig->upper_eta) {
            triglist.push_back(new Trigger(*trig));
        }
    }
    // Then sort the triggers by jet momentum (pt) cutoff so that we can assign the upper pt cutoff
    int ntrigs = triglist.size();
    for (int n1 = 0; n1 < ntrigs; n1++) {
        int nextminjetind = 0;
        for (int n2 = n1; n2 < ntrigs; n2++) {
            if (triglist[n2]->min_pt < triglist[nextminjetind]->min_pt) nextminjetind = n2;
        }
        swap(triglist[nextminjetind], triglist[n1]);
    }
    // Now that the triggers have been sorted, eassign upper pt cuts.
    for (int n = 0; n < ntrigs-1; n++) {
        triglist[n]->max_pt = triglist[n+1]->min_pt;
    }
    triglist[ntrigs-1]->max_pt = global_max_pt;

    return triglist;
}


/**
 * Initializes triggers complete with momentum and pseudorapidity cutoffs.
 * To add new triggers, simply create a line like the one below with the trigger name, the momentum threshold, and its pseudorapidity interval.
 */
void initialize () {

    // Create a vector of triggers.
    trigger_vec.push_back(new Trigger("HLT_j15_p320eta490_L1MBTS_1_1", 15+trigthres, 3.2, 4.9));
    trigger_vec.push_back(new Trigger("HLT_j30_0eta490_L1TE10", 30+trigthres, -4.9, 4.9));
    trigger_vec.push_back(new Trigger("HLT_j40_L1J5", 40+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j45_p200eta320", 45+trigthres, 2, 3.2));
    trigger_vec.push_back(new Trigger("HLT_j50_L1J10", 50+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j55_p320eta490", 55+trigthres, 3.2, 4.9));
    trigger_vec.push_back(new Trigger("HLT_j60", 60+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j65_p200eta320", 65+trigthres, 2, 3.2));
    trigger_vec.push_back(new Trigger("HLT_j75_L1J20", 75+trigthres, min_eta, max_eta));
    trigger_vec.push_back(new Trigger("HLT_j100_L1J20", 100+trigthres, min_eta, max_eta));

    numtrigs = trigger_vec.size();

    // Assign indices for tree branching.
    for (int i = 0; i < numtrigs; i++) {
        trigger_vec[i]->index = i;
    }

    // Instantiate the pseudorapidity interval trigger vectors.
    triggers_n490eta320 = get_trigs_eta_range(-4.9, -3.2);
    triggers_n320eta200 = get_trigs_eta_range(-3.2, -2.0);
    triggers_n200eta100 = get_trigs_eta_range(-2.0, -1.0);
    triggers_n100eta0 = get_trigs_eta_range(-1.0, 0);
    triggers_p0eta100 = get_trigs_eta_range(0, 1.0);
    triggers_p100eta200 = get_trigs_eta_range(1.0, 2.0);
    triggers_p200eta320 = get_trigs_eta_range(2.0, 3.2);
    triggers_p320eta490 = get_trigs_eta_range(3.2, 4.9);
  
}
