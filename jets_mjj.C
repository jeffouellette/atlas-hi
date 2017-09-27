#include "triggerUtil.C"

void jets_mjj(const int runNumber, // Run number identifier.
        double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{

    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    const int numhists = 5;

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

    //Useful arrays

    const double xbins [17] = {25, 30, 40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};
    const int len_xbins = sizeof(xbins)/sizeof(xbins[0]);

    const double harr_scales[5] = {1, 1, 1, 1, 1};   // rescaling factors so the histograms don't overlap

    const float etastarcuts[numhists+1] = {0, 0.5, 1, 1.5, 2, 3};
    double d_eta[numhists];
    for (int i = 0; i < numhists; i++) {
        d_eta[i] = etastarcuts[i+1] - etastarcuts[i];
    }

    // Create an array of 6 histograms, one for each rapidity region.   
    TH1D* harr[numhists];
    for (int i = 0; i < numhists; i++) {
        harr[i] = new TH1D(Form("%ieta%i", runNumber, i), Form("%g #leq #left|#eta*#right| #leq %g; #it{M}_{JJ} #left[GeV/#it{c}^{2}#right];d^{2}#sigma/d#it{M}_{JJ}d#left|#eta*#right| #left[pb (GeV/#it{c}^{#it{2}})^{-1}#right]", etastarcuts[i], etastarcuts[i+1]), len_xbins-1, xbins);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }

    // Create branching addresses:  
    // Create arrays to store trigger values for each event
    bool m_trig_bool[trigLength];   // stores whether trigger was triggered
    float m_trig_prescale[trigLength];      // stores the prescaling factor for the trigger
    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    float j_phi[60] = {};
    float j_e[60] = {};
    int njet = 0;

    // Set branch addresses
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("j_phi", j_phi);
    tree->SetBranchAddress("j_e", j_e);
    tree->SetBranchAddress("njet", &njet);
    for (int i = 0; i < trigLength; i++) {
        tree->SetBranchAddress(m_trig_string[i], &m_trig_bool[i]); 
        tree->SetBranchAddress(Form("%s_prescale", m_trig_string[i]), &m_trig_prescale[i]);
    }


    // Iterate over each event
    const int numentries = tree->GetEntries();

    // Dijet cone radius for event selection
    //const double r_squared = 1;

    // |eta| < 2 selections add to histograms 2, 3, 4, & 5.
    // 2 <= |eta| < 3.2 selections add to histograms 1 & 6.
    // 3.2 <= |eta| < 4.9 selections add to histograms 0 & 7.
    const std::vector<int> trig_n200eta490 = {3, 5, 9, 13, 17, 19};
    const std::vector<int> trig_0eta200 = {3, 5, 9, 13, 17, 19};
    const std::vector<int> trig_p200eta320 = {3, 5, 7, 9, 13, 15, 17, 19};
    const std::vector<int> trig_p320eta490 = {1, 3, 5, 9, 11, 13, 17, 19};
    double jeta0, jeta1, jphi0, jphi1, jpt0, jpt1, je0, je1, etastar, mjj, jpz0, jpz1, extra_jpt_sum;
    TLorentzVector jet0, jet1;
    bool takeEvent;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses
        // Convert to doubles to obtain higher precision.

        jpt0 = (double)j_pt[0];
        jphi0 = (double)j_phi[0];
        jeta0 = (double)j_eta[0];
        jpt1 = (double)j_pt[1];
        jphi1 = (double)j_phi[1];
        jeta1 = (double)j_eta[1];
        
        extra_jpt_sum = 0;
        for (int j = 2; j < njet; j++) {
            extra_jpt_sum += j_pt[j];
        }
        takeEvent = (extra_jpt_sum / (jpt0 + jpt1 + extra_jpt_sum) <= dijet_pt_frac_cutoff);
        if (takeEvent) {    // select 2 jet events or 3 jet events with 2 dominating jets

            je0 = (double)j_e[0];
            je1 = (double)j_e[1];

            etastar = TMath::Abs(jeta0-jeta1)/2;

            //mjj = 2 * TMath::Abs(0.5*(jpt0+jpt1)) * TMath::CosH(etastar);
            jet0.SetPtEtaPhiE(jpt0, jeta0, jphi0, je0);
            jet1.SetPtEtaPhiE(jpt1, jeta1, jphi1, je1);
            mjj = get_mjj(jet0, jet1); 

            if (TMath::Abs(jeta0) < 2 && TMath::Abs(jeta0) >= 0) {
                for (int trig_num : trig_0eta200) { // iterate over each trigger
                    if (m_trig_bool[trig_num] && jpt0 >= trig_lower_0eta200[trig_num] && jpt0 < trig_upper_0eta200[trig_num]) { // if triggered, check the jet momentum
                        for (int i = 0; i < numhists; i++) {
                            if (etastarcuts[i] <= etastar && etastar < etastarcuts[i+1]) {
                                    harr[i]->Fill(mjj, (double)m_trig_prescale[trig_num]);
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
                        for (int i = 0; i < numhists; i++) {
                             if (etastarcuts[i] <= etastar && etastar < etastarcuts[i+1]) {
                                 harr[i]->Fill(mjj, (double)m_trig_prescale[trig_num]);
                                 break;
                             }
                        }
                        break;
                    }
                }
            }
            else if (jeta0 < 4.9 && jeta0 >= 3.2) {
                for (int trig_num : trig_p320eta490) {
                    if (m_trig_bool[trig_num] && jpt0 >= trig_lower_p320eta490[trig_num] && jpt0 < trig_upper_p320eta490[trig_num]) {
                        for (int i = 0; i < numhists; i++) {
                            if (etastarcuts[i] <= etastar && etastar < etastarcuts[i+1]) {
                                harr[i]->Fill(mjj, (double)m_trig_prescale[trig_num]);
                                break;
                            }
                        }
                        break;
                    }
                }
            }
            else if (jeta0 <= -2 && jeta0 > -3.2) {
                for (int trig_num : trig_n200eta490) {
                    if (m_trig_bool[trig_num] && jpt0 >= trig_lower_n200eta490[trig_num] && jpt0 < trig_upper_n200eta490[trig_num]) {
                        for (int i = 0; i < numhists; i++) {
                            if (etastarcuts[i] <= etastar && etastar < etastarcuts[i+1]) {
                                harr[i]->Fill(mjj, (double)m_trig_prescale[trig_num]);
                                break;
                            }
                        }
                        break;
                    }
                }
            }
            else if (jeta0 <= -3.2 && jeta0 > -4.9) {
                for (int trig_num : trig_n200eta490) {
                    if (m_trig_bool[trig_num] && jpt0 >= trig_lower_n200eta490[trig_num] && jpt0 < trig_upper_n200eta490[trig_num]) {
                        for (int i = 0; i < numhists; i++) {
                            if (etastarcuts[i] <= etastar && etastar < etastarcuts[i+1]) {
                                harr[i]->Fill(mjj, (double)m_trig_prescale[trig_num]);
                                break;
                            }
                        }
                        break;
                    }
                }
            }
        }
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("./mjj_data/run_%i.root", runNumber), "RECREATE");
    for (int i = 0; i < numhists; i++) {
        harr[i]->Scale((harr_scales[i]) / (A * luminosity * d_eta[i]), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        harr[i]->Write();
    }

    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    output->Close();
}
