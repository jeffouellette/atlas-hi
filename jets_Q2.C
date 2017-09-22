#include "triggerUtil.C"

void jets_Q2(const int runNumber, // Run number identifier.
             float luminosity, // Integrated luminosity for this run. Presumed constant over the run period.
             const bool PbP = false) // Stores if this run was a Pb-p run (opposite direction). This makes the lab frame boost in the opposite direction. 
{

    float eta_lab;
    if (PbP) eta_lab = -0.465;
    else eta_lab = 0.465;

    const int numhists = 8;
    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    // Initialize maps of trigger numbers (as ordered above) to the appropriate jet cutoffs for a specific eta range. This is done to obtain continuous coverage over the p_t spectrum.

    std::map<int, int> trig_lower_n200eta490 = get_trig_lower_n200eta490(); // Eta range : -4.9 < eta <= -3.2
    std::map<int, int> trig_upper_n200eta490 = get_trig_upper_n200eta490();
    std::map<int, int> trig_lower_0eta200 = get_trig_lower_0eta200();       // Eta range -2 < eta < 2
    std::map<int, int> trig_upper_0eta200 = get_trig_upper_0eta200();
    std::map<int, int> trig_lower_p200eta320 = get_trig_lower_p200eta320(); // Eta range: 2 <= eta < 3.2
    std::map<int, int> trig_upper_p200eta320 = get_trig_upper_p200eta320();
    std::map<int, int> trig_lower_p320eta490 = get_trig_lower_p320eta490(); // Eta range: 3.2 <= eta < 4.9
    std::map<int, int> trig_upper_p320eta490 = get_trig_upper_p320eta490();

    // End trigger-jet cutoff map initialization


    TTree* tree = (TTree*)(new TFile(Form("./rundata/run_%i_raw.root", runNumber)))->Get("tree");

    //Useful arrays
    const double xbins[17] = {25, 30, 40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};
    const float jet_cuts_n200eta490[7] = {40, 50, 60, 70, 85, 110, 6000};
    const float jet_cuts_0eta200[7] = {40, 50, 60, 70, 85, 110, 6000};
    const float jet_cuts_p200eta320[9] = {40, 50, 55, 60, 70, 75, 85, 110, 6000};
    const float jet_cuts_p320eta490[9] = {25, 40, 50, 60, 65, 70, 85, 110, 6000};

    const double d_eta[8] = {1.7, 1.2, 1, 1, 1, 1, 1.2, 1.7};
    const float eta_cuts[9] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
    const double harr_scales[8] = {0.005, 0.03, 0.1, 0.5, 1, 0.3, 0.05, 0.01};   // rescaling factors so the histograms don't overlap
    //const double harr_scales[8] = {1, 1, 1, 1, 1, 1, 1, 1};

    // Create an array of 6 histograms, one for each rapidity region.   
    TH1D* harr[numhists];
    for (int i = 0; i < numhists; i++) {
        harr[i] = new TH1D(Form("%ieta%i", runNumber, i), Form("%1.1f < #eta < %1.1f (#times %1.3f);#it{Q}^{dijet}_{12} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{Q}_{12}dy #left[pb (GeV/#it{c})^{-1}#right]", eta_cuts[i], eta_cuts[i+1], harr_scales[i]), sizeof(xbins)/sizeof(xbins[0])-1, xbins);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }

    // Create branching addresses:  
    // Create arrays to store trigger values for each event
    bool m_trig_bool[trigLength];   // stores whether trigger was triggered
    float m_trig_prescale[trigLength];      // stores the prescaling factor for the trigger
    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    float j_e[60] = {};
    int njet = 0;

    // Set branch addresses
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("j_e", j_e);
    tree->SetBranchAddress("njet", &njet);
    for (int i = 0; i < trigLength; i++) {
        tree->SetBranchAddress(m_trig_string[i], &m_trig_bool[i]); 
        tree->SetBranchAddress(Form("%s_prescale", m_trig_string[i]), &m_trig_prescale[i]);
    }


    // Iterate over each event
    const int numentries = tree->GetEntries();

    // |eta| < 2 selections add to histograms 2, 3, 4, & 5.
    // 2 <= |eta| < 3.2 selections add to histograms 1 & 6.
    // 3.2 <= |eta| < 4.9 selections add to histograms 0 & 7.
    const std::vector<int> trig_n200eta490 = {3, 5, 9, 13, 17, 19};
    const std::vector<int> trig_0eta200 = {3, 5, 9, 13, 17, 19};
    const std::vector<int> trig_p200eta320 = {3, 5, 7, 9, 13, 15, 17, 19};
    const std::vector<int> trig_p320eta490 = {1, 3, 5, 9, 11, 13, 17, 19};
    double jeta, jpt, q2low, q2high, extra_jpt_sum;
    bool takeEvent;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses

        extra_jpt_sum = 0;
        for (int j = 2; j < njet; j++) {
            extra_jpt_sum += j_pt[j];
        }
        takeEvent = extra_jpt_sum / (j_pt[0] + j_pt[1] + extra_jpt_sum) <= dijet_pt_frac_cutoff;
        if (takeEvent) { // select 2 jet events or events with 2 dominating jets

            jpt = (double)j_pt[0];
            jeta = (double)j_eta[0];

            if (TMath::Abs(jeta) < 2 && TMath::Abs(jeta) >= 0) {
                for (int trig_num : trig_0eta200) { // iterate over each trigger
                    if (m_trig_bool[trig_num] && jpt >= trig_lower_0eta200[trig_num] && jpt < trig_upper_0eta200[trig_num]) { // if triggered, check whether the jet momentum falls in the correct range
                        for (int k = 2; k < 6; k++) {
                            if (jeta >= eta_cuts[k] && jeta < eta_cuts[k+1]) {
                                float xp = (TMath::Sqrt(Z/A) / sqrt_s_nn) * (j_pt[0]*TMath::Exp(j_eta[0]+eta_lab)+j_pt[1]*TMath::Exp(j_eta[1]+eta_lab)); 
                                q2low = get_q2(xp, (double)j_e[0], (double)j_pt[0]); 
                                q2high = get_q2(xp, (double)j_e[1], (double)j_pt[1]);
                                harr[k]->Fill(0.5*(q2high+q2low), m_trig_prescale[trig_num]);
                                break;
                            }
                        }
                        break; // any jet should only be plotted once.
                    }
                }
            }
            else if (jeta < 3.2 && jeta >= 2) {
                for (int trig_num : trig_p200eta320) {
                    if (m_trig_bool[trig_num] && jpt >= trig_lower_p200eta320[trig_num] && jpt < trig_upper_p200eta320[trig_num]) {
                        float xp = (TMath::Sqrt(Z/A) / sqrt_s_nn) * (j_pt[0]*TMath::Exp(j_eta[0]+eta_lab)+j_pt[1]*TMath::Exp(j_eta[1]+eta_lab)); 
                        q2low = get_q2(xp, (double)j_e[0], (double)j_pt[0]); 
                        q2high = get_q2(xp, (double)j_e[1], (double)j_pt[1]);
                        harr[6]->Fill(0.5*(q2high+q2low), m_trig_prescale[trig_num]);
                        break;
                    }
                }
            }
            else if (jeta < 4.9 && jeta >= 3.2) {
                for (int trig_num : trig_p320eta490) {
                    if (m_trig_bool[trig_num] && jpt >= trig_lower_p320eta490[trig_num] && jpt < trig_upper_p320eta490[trig_num]) {
                        float xp = (TMath::Sqrt(Z/A) / sqrt_s_nn) * (j_pt[0]*TMath::Exp(j_eta[0]+eta_lab)+j_pt[1]*TMath::Exp(j_eta[1]+eta_lab)); 
                        q2low = get_q2(xp, (double)j_e[0], (double)j_pt[0]); 
                        q2high = get_q2(xp, (double)j_e[1], (double)j_pt[1]);
                        harr[7]->Fill(0.5*(q2high+q2low), m_trig_prescale[trig_num]);
                        break;
                    }
                }
            }
            else if (jeta <= -2 && jeta > -3.2) {
                for (int trig_num : trig_n200eta490) {
                    if (m_trig_bool[trig_num] && jpt >= trig_lower_n200eta490[trig_num] && jpt < trig_upper_n200eta490[trig_num]) {
                        float xp = (TMath::Sqrt(Z/A) / sqrt_s_nn) * (j_pt[0]*TMath::Exp(j_eta[0]+eta_lab)+j_pt[1]*TMath::Exp(j_eta[1]+eta_lab)); 
                        q2low = get_q2(xp, (double)j_e[0], (double)j_pt[0]); 
                        q2high = get_q2(xp, (double)j_e[1], (double)j_pt[1]);
                        harr[1]->Fill(0.5*(q2high+q2low), m_trig_prescale[trig_num]);
                        break;
                    }
                }
            }
            else if (jeta <= -3.2 && jeta > -4.9) {
                for (int trig_num : trig_n200eta490) {
                    if (m_trig_bool[trig_num] && jpt >= trig_lower_n200eta490[trig_num] && jpt < trig_upper_n200eta490[trig_num]) {
                        float xp = (TMath::Sqrt(Z/A) / sqrt_s_nn) * (j_pt[0]*TMath::Exp(j_eta[0]+eta_lab)+j_pt[1]*TMath::Exp(j_eta[1]+eta_lab)); 
                        q2low = get_q2(xp, (double)j_e[0], (double)j_pt[0]); 
                        q2high = get_q2(xp, (double)j_e[1], (double)j_pt[1]);
                        harr[0]->Fill(0.5*(q2high+q2low), m_trig_prescale[trig_num]);
                        break;
                    }
                }
            }
        }
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("./Q2_data/run_%i.root", runNumber), "RECREATE");
    for (int i = 0; i< numhists; i++) {
        harr[i]->Scale((harr_scales[i]) / (A * luminosity * d_eta[i]), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        harr[i]->Write();
    }
    output->Close();
}
