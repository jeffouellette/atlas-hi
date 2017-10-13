#include "../triggerUtil.C"

void jets_Q2(const int thisRunNumber, // Run number identifier.
             double luminosity, // Integrated luminosity for this run. Presumed constant over the run period.
             bool periodA)
{


    const int numhists = 8;
    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    initialize(thisRunNumber);

    TTree* tree = (TTree*)(new TFile(Form("./rundata/run_%i_raw.root", thisRunNumber)))->Get("tree");

    //Useful arrays
    const double xbins[42] = {25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 220., 240., 260., 280., 300., 350., 400., 500., 600., 800., 1100., 1500., 2000., 2500., 6000.};
    const int numbins = sizeof(xbins)/sizeof(xbins[0]) - 1;
    const double d_eta[8] = {1.7, 1.2, 1, 1, 1, 1, 1.2, 1.7};
    const double eta_cuts[9] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
    const double harr_scales[8] = {0.005, 0.03, 0.1, 0.5, 1, 0.3, 0.05, 0.01};   // rescaling factors so the histograms don't overlap

    // Create an array of 6 histograms, one for each rapidity region.   
    TH1D* harr[numhists];
    for (int i = 0; i < numhists; i++) {
        harr[i] = new TH1D(Form("%ieta%i", thisRunNumber, i), Form("%g < #eta < %g (#times %g);#it{Q}^{avg}_{JJ} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{Q}_{JJ}dy #left[pb (GeV/#it{c})^{-1}#right]", eta_cuts[i], eta_cuts[i+1], harr_scales[i]), numbins, xbins);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }

    // Create branching addresses:  
    // Create arrays to store trigger values for each event
    bool m_trig_bool[numtrigs];   // stores whether trigger was triggered
    float m_trig_prescale[numtrigs];      // stores the prescaling factor for the trigger
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
    for (Trigger* trig : trigger_vec) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &m_trig_bool[trig->index]);
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &m_trig_prescale[trig->index]);
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();

    double jpt0, jpt1, jeta0, jeta1, je0, je1, xp, q2, extra_jpt_sum;
    bool takeEvent;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses

        extra_jpt_sum = 0;
        for (int j = 2; j < njet; j++) {
            extra_jpt_sum += j_pt[j];
        }
        takeEvent = extra_jpt_sum / (j_pt[0] + j_pt[1] + extra_jpt_sum) <= dijet_pt_frac_cutoff;
        if (takeEvent) { // select 2 jet events or events with 2 dominating jets

            jpt0 = (double)j_pt[0];
            jpt1 = (double)j_pt[1];
            jeta0 = (double)j_eta[0];
            jeta1 = (double)j_eta[1];
            je0 = (double)j_e[0];
            je1 = (double)j_e[1];
            xp = get_xp(jpt0, jpt1, jeta0, jeta1, periodA); 
            q2 = 0.5 * (get_q2(xp, je0, jpt0) + get_q2(xp, je1, jpt1));


            if (0 < jeta0 && jeta0 < 1) {
                for (Trigger* trig : triggers_p0eta100) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0; // Only take the event if triggered and trigger was not disabled
                    if (takeEvent && trig->min_pt <= jpt0 && jpt0 < trig->max_pt) {
                        harr[4]->Fill(q2, m_trig_prescale[trig->index]);
                        break; // Break to ensure that only one trigger allows the event to be recorded
                    }
                }
            }
            else if (-1 < jeta0 && jeta0 < 0) {
                for (Trigger* trig : triggers_n100eta0) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt0 && jpt0 < trig->max_pt) {
                        harr[3]->Fill(q2, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (1 < jeta0 && jeta0 < 2) {
                for (Trigger* trig : triggers_p100eta200) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt0 && jpt0 < trig->max_pt) {
                        harr[5]->Fill(q2, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (-2 < jeta0 && jeta0 < -1) {
                for (Trigger* trig : triggers_n200eta100) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt0 && jpt0 < trig->max_pt) {
                        harr[2]->Fill(q2, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (2 < jeta0 && jeta0 < 3.2) {
                for (Trigger* trig : triggers_p200eta320) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt0 && jpt0 < trig->max_pt) {
                        harr[6]->Fill(q2, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (-3.2 < jeta0 && jeta0 < -2) {
                for (Trigger* trig : triggers_n320eta200) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt0 && jpt0 < trig->max_pt) {
                        harr[1]->Fill(q2, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (3.2 < jeta0 && jeta0 < 4.9) {
                for (Trigger* trig : triggers_p320eta490) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt0 && jpt0 < trig->max_pt) {
                        harr[7]->Fill(q2, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (-4.9 < jeta0 && jeta0 < -3.2) {
                for (Trigger* trig : triggers_n490eta320) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt0 && jpt0 < trig->max_pt) {
                        harr[0]->Fill(q2, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
        }
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("./Q2_data/run_%i.root", thisRunNumber), "RECREATE");
    for (int i = 0; i< numhists; i++) {
        harr[i]->Scale((harr_scales[i]) / (A * luminosity * d_eta[i]), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        harr[i]->Write();
    }
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    output->Close();
}
