#include "triggerUtil.C"

void jets_xa_xp(int runNumber, // Run number identifier.
                float luminosity, // Integrated luminosity for this run. Presumed constant over the run period.
                bool periodA)
{

    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    // Initialize maps of trigger numbers (as ordered above) to the appropriate jet cutoffs for a specific eta range. This is done to obtain continuous coverage over the p_t spectrum.
    std::map<int, int> trig_lower_n200eta490 = get_trig_lower_n200eta490(); // Eta range : -4.9 < eta <= -3.2
    std::map<int, int> trig_upper_n200eta490 = get_trig_upper_n200eta490();
    std::map<int, int> trig_lower_0eta200 = get_trig_lower_0eta200(); // Eta range -2 < eta < 2
    std::map<int, int> trig_upper_0eta200 = get_trig_upper_0eta200();
    std::map<int, int> trig_lower_p200eta320 = get_trig_lower_p200eta320(); // Eta range: 2 <= eta < 3.2
    std::map<int, int> trig_upper_p200eta320 = get_trig_upper_p200eta320();
    std::map<int, int> trig_lower_p320eta490 = get_trig_lower_p320eta490(); // Eta range: 3.2 <= eta < 4.9
    std::map<int, int> trig_upper_p320eta490 = get_trig_upper_p320eta490();

    // End trigger-jet cutoff map initialization

    TTree* tree = (TTree*)(new TFile(Form("./rundata/run_%i_raw.root", runNumber)))->Get("tree");

    const float jet_cuts_n200eta490[7] = {40, 50, 60, 70, 85, 110, 2000};
    const float jet_cuts_0eta200[7] = {40, 50, 60, 70, 85, 110, 2000};
    const float jet_cuts_p200eta320[9] = {40, 50, 55, 60, 70, 75, 85, 110, 2000};
    const float jet_cuts_p320eta490[9] = {25, 40, 50, 60, 65, 70, 85, 110, 2000};
    
    const float d_eta[8] = {1.7, 1.2, 1, 1, 1, 1, 1.2, 1.7};
    const float eta_cuts[9] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
    const float harr_scales[8] = {1, 1, 1, 1, 1, 1, 1, 1};   // rescaling factors so the histograms don't overlap

    const float xbins[34] = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.80, 0.84, 0.88, 0.92, 0.96, 1.00, 1.08, 1.16, 1.24, 1.32, 1.40, 1.48, 1.56, 1.64};
    const int nbins = sizeof(xbins)/sizeof(xbins[0]) - 1; 
    // Create an array of 16 histograms, one for each rapidity region and one for x_p, x_a. 
    const int numhists = 16;
    TH1F* harr[numhists];
    for (int i = 0; i < numhists/2; i++) {
        harr[i] = new TH1F(Form("%ieta%i", runNumber, i), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d#sigma^{2}/d#it{x}_{p} dy #left[pb#right]", eta_cuts[i], eta_cuts[i+1]), nbins, xbins);
        harr[i]->Sumw2();
    }
    for (int i = numhists/2; i < numhists; i++) {
        harr[i] = new TH1F(Form("%ieta%i", runNumber, i), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d#sigma^{2}/d#it{x}_{a} dy #left[pb#right]", eta_cuts[i%(numhists/2)], eta_cuts[(i%(numhists/2))+1]), nbins, xbins);
        harr[i]->Sumw2();
    }

    // Create arrays to store trigger values for each event
    bool m_trig_bool[trigLength];
    float m_trig_prescale[trigLength];

    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    float j_phi[60] = {};
    float j_e[60] = {};
    int njet = 0;
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("j_e", j_e);
    tree->SetBranchAddress("njet", &njet);
    tree->SetBranchAddress("j_phi", j_phi);

    // Set branch addresses
    for (int i = 0; i < trigLength; i++) {
        tree->SetBranchAddress(m_trig_string[i], &m_trig_bool[i]); 
        tree->SetBranchAddress(Form("%s_prescale", m_trig_string[i]), &m_trig_prescale[i]);
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();
        
    const std::vector<int> trig_n200eta490 = {3, 5, 9, 13, 17, 19};
    const std::vector<int> trig_0eta200 = {3, 5, 9, 13, 17, 19};
    const std::vector<int> trig_p200eta320 = {3, 5, 7, 9, 13, 15, 17, 19};
    const std::vector<int> trig_p320eta490 = {1, 3, 5, 9, 11, 13, 17, 19};
    double jeta0, jpt0, jeta1, jpt1, xp, xa, extra_jpt_sum;
    bool takeEvent;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses

        extra_jpt_sum = 0;
        for (int j = 2; j < njet; j++) {
            extra_jpt_sum += j_pt[j];
        }
        takeEvent = extra_jpt_sum / (j_pt[0] + j_pt[1] + extra_jpt_sum) <= dijet_pt_frac_cutoff;
        if (takeEvent) {

            jpt0 = (double)j_pt[0];
            jeta0 = (double)j_eta[0];
            jpt1 = (double)j_pt[1];
            jeta1 = (double)j_eta[1];

            xp = get_xp(jpt0, jpt1, jeta0, jeta1, periodA);
            xa = get_xa(jpt0, jpt1, jeta0, jeta1, periodA);

            if (TMath::Abs(jeta0) < 2 && TMath::Abs(jeta0) >= 0) {
                for (int trig_num : trig_0eta200) { // iterate over each trigger
                    if (m_trig_bool[trig_num] && jpt0 >= trig_lower_0eta200[trig_num] && jpt0 < trig_upper_0eta200[trig_num]) { // if triggered, check whether the jet momentum falls in the correct range
                        for (int k = 2; k < 6; k++) {
                            if (jeta0 >= eta_cuts[k] && jeta0 < eta_cuts[k+1]) {
                                harr[k]->Fill(xp, m_trig_prescale[trig_num]);
                                harr[k+8]->Fill(xa, m_trig_prescale[trig_num]);
                                break;
                            }
                        }
                        break; // any jet should only be plotted once.
                    }
                }
            }
            else if (jeta0 < 3.2 && jeta0 >= 2) {
                for (int trig_num : trig_p200eta320) {
                    if (m_trig_bool[trig_num] && jpt0 >= trig_lower_p200eta320[trig_num] && jpt0 < trig_upper_p200eta320[trig_num]) {
                        harr[6]->Fill(xp, m_trig_prescale[trig_num]);
                        harr[14]->Fill(xa, m_trig_prescale[trig_num]);
                        break;
                    }
                }
            }
            else if (jeta0 < 4.9 && jeta0 >= 3.2) {
                for (int trig_num : trig_p320eta490) {
                    if (m_trig_bool[trig_num] && jpt0 >= trig_lower_p320eta490[trig_num] && jpt0 < trig_upper_p320eta490[trig_num]) {
                        harr[7]->Fill(xp, m_trig_prescale[trig_num]);
                        harr[15]->Fill(xa, m_trig_prescale[trig_num]);
                        break;
                    }
                }
            }
            else if (jeta0 <= -2 && jeta0 > -3.2) {
                for (int trig_num : trig_n200eta490) {
                    if (m_trig_bool[trig_num] && jpt0 >= trig_lower_n200eta490[trig_num] && jpt0 < trig_upper_n200eta490[trig_num]) {
                        harr[1]->Fill(xp, m_trig_prescale[trig_num]);
                        harr[9]->Fill(xa, m_trig_prescale[trig_num]);
                        break;
                    }
                }
            }
            else if (jeta0 <= -3.2 && jeta0 > -4.9) {
                for (int trig_num : trig_n200eta490) {
                    if (m_trig_bool[trig_num] && jpt0 >= trig_lower_n200eta490[trig_num] && jpt0 < trig_upper_n200eta490[trig_num]) {
                        harr[0]->Fill(xp, m_trig_prescale[trig_num]);
                        harr[8]->Fill(xa, m_trig_prescale[trig_num]);
                        break;
                    }
                }
            }
        }
    }
    // Save to root file
    TFile* output = new TFile(Form("./xdata/run_%i.root", runNumber), "RECREATE");
    for (int i = 0; i < numhists; i++) {
        harr[i]->Scale(harr_scales[i%(numhists/2)]/(A * luminosity * d_eta[i%(numhists/2)]), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        harr[i]->Write();
    }
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    output->Close();
}
