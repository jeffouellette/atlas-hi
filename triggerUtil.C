
const float dijet_pt_frac_cutoff = 0.1; // Maximum fraction of transverse momentum of all jets but the first two for the event to be considered a dijet

const float Z = 82;   // value of Z for Pb
const float A = 208;  // value of A for Pb
const float sqrt_s_nn = 8160; // p-Pb collision energy in CoM frame (GeV)
// Store trigger names as an array of strings and loop over them to set the branch addresses
static const int trigLength = 24;
const char* m_trig_string[trigLength] = {
    "HLT_j15_p320eta490_L1MBTS_1_1",
    "HLT_j15_ion_p320eta490_L1MBTS_1_1",
    "HLT_j30_0eta490_L1TE10",
    "HLT_j30_ion_0eta490_L1TE10",
    "HLT_j40_L1J5",
    "HLT_j40_ion_L1J5",
    "HLT_j45_p200eta320",
    "HLT_j45_ion_p200eta320",
    "HLT_j50_L1J10",
    "HLT_j50_ion_L1J10",
    "HLT_j55_p320eta490",
    "HLT_j55_ion_p320eta490",
    "HLT_j60",
    "HLT_j60_ion_L1J20",
    "HLT_j65_p200eta320",
    "HLT_j65_ion_p200eta320",
    "HLT_j75_L1J20",
    "HLT_j75_ion_L1J20",
    "HLT_j100_L1J20",
    "HLT_j100_ion_L1J20",
    "HLT_2j10_p320eta490_L1TE10",
    "HLT_2j10_ion_p320eta490_L1TE10",
    "HLT_2j30_p320eta490",
    "HLT_2j30_ion_p320eta490"
};

//TTree* tree;
                               
double get_xp(double jpt0, double jpt1, double jeta0, double jeta1) {
        return (TMath::Sqrt(Z/A) / sqrt_s_nn) * (jpt0* TMath::Exp(jeta0)+jpt1*TMath::Exp(jeta1));
}

double get_xa(double jpt0, double jpt1, double jeta0, double jeta1) {
        return (TMath::Sqrt(A/Z) / sqrt_s_nn) * (jpt0*TMath::Exp(-jeta0)+jpt1*TMath::Exp(-jeta1));
}

double get_q2(double xp, double je, double jpt){
    return (double)TMath::Sqrt(TMath::Sqrt(A/Z)*sqrt_s_nn*xp*(je-TMath::Sqrt(je*je-jpt*jpt)));
}

double get_mjj(TLorentzVector jet0, TLorentzVector jet1) {
    return (jet0+jet1).Mag();
}

// Initialize maps of trigger numbers (as ordered above) to the appropriate jet cutoffs for a specific eta range. This is done to obtain continuous coverage over the p_t spectrum.

std::map<int, int> get_trig_lower_n200eta490 () {
    std::map<int, int> trig_lower_n200eta490;
    trig_lower_n200eta490[3] = 40;
    trig_lower_n200eta490[5] = 50;
    trig_lower_n200eta490[9] = 60;
    trig_lower_n200eta490[13] = 70;
    trig_lower_n200eta490[17] = 85;
    trig_lower_n200eta490[19] = 110;
    return trig_lower_n200eta490; 
}


std::map<int, int> get_trig_upper_n200eta490 () {
    std::map<int, int> trig_upper_n200eta490;
    trig_upper_n200eta490[3] = 50;
    trig_upper_n200eta490[5] = 60;
    trig_upper_n200eta490[9] = 70;
    trig_upper_n200eta490[13] = 85;
    trig_upper_n200eta490[17] = 110;
    trig_upper_n200eta490[19] = 2000;
    return trig_upper_n200eta490;
}

std::map<int, int> get_trig_lower_0eta200 () {
    std::map<int, int> trig_lower_0eta200; // Eta range -2 < eta < 2
    trig_lower_0eta200[3] = 40;
    trig_lower_0eta200[5] = 50;
    trig_lower_0eta200[9] = 60;
    trig_lower_0eta200[13] = 70;
    trig_lower_0eta200[17] = 85;
    trig_lower_0eta200[19] = 110;
    return trig_lower_0eta200;
}

std::map<int, int> get_trig_upper_0eta200 () {
    std::map<int, int> trig_upper_0eta200;
    trig_upper_0eta200[3] = 50;
    trig_upper_0eta200[5] = 60;
    trig_upper_0eta200[9] = 70;
    trig_upper_0eta200[13] = 85;
    trig_upper_0eta200[17] = 110;
    trig_upper_0eta200[19] = 2000;
    return trig_upper_0eta200;
}

std::map<int, int> get_trig_lower_p200eta320 () {
    std::map<int, int> trig_lower_p200eta320; // Eta range: 2 <= eta < 3.2
    trig_lower_p200eta320[3] = 40;
    trig_lower_p200eta320[5] = 50;
    trig_lower_p200eta320[7] = 55;
    trig_lower_p200eta320[9] = 60;
    trig_lower_p200eta320[13] = 70;
    trig_lower_p200eta320[15] = 75;
    trig_lower_p200eta320[17] = 85;
    trig_lower_p200eta320[19] = 110;
    return trig_lower_p200eta320;
}

std::map<int, int> get_trig_upper_p200eta320 () {
    std::map<int, int> trig_upper_p200eta320;
    trig_upper_p200eta320[3] = 50;
    trig_upper_p200eta320[5] = 55;
    trig_upper_p200eta320[7] = 60;
    trig_upper_p200eta320[9] = 70;
    trig_upper_p200eta320[13] = 75;
    trig_upper_p200eta320[15] = 85;
    trig_upper_p200eta320[17] = 110;
    trig_upper_p200eta320[19] = 2000;
    return trig_upper_p200eta320;
}

std::map<int, int> get_trig_lower_p320eta490 () {
    std::map<int, int> trig_lower_p320eta490; // Eta range: 3.2 <= eta < 4.9
    trig_lower_p320eta490[1] = 25;
    trig_lower_p320eta490[3] = 40;
    trig_lower_p320eta490[5] = 50;
    trig_lower_p320eta490[9] = 60;
    trig_lower_p320eta490[11] = 65;
    trig_lower_p320eta490[13] = 70;
    trig_lower_p320eta490[17] = 85;
    trig_lower_p320eta490[19] = 110;
    return trig_lower_p320eta490;        
}

std::map<int, int> get_trig_upper_p320eta490 () {
    std::map<int, int> trig_upper_p320eta490;
    trig_upper_p320eta490[1] = 40;
    trig_upper_p320eta490[3] = 50;
    trig_upper_p320eta490[5] = 60;
    trig_upper_p320eta490[9] = 65;
    trig_upper_p320eta490[11] = 70;
    trig_upper_p320eta490[13] = 85;
    trig_upper_p320eta490[17] = 110;
    trig_upper_p320eta490[19] = 2000;
    return trig_upper_p320eta490;
}
// End trigger-jet cutoff map initialization


