#include "../triggerUtil.C"

void jets_xa_xp(int thisRunNumber, // Run number identifier.
                double luminosity, // Integrated luminosity for this run. Presumed constant over the run period.
                bool periodA)
{

    const int numbins = 40;
    const int numhists = 16;
    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    initialize(thisRunNumber);

    TTree* tree = (TTree*)(new TFile(Form("./rundata/run_%i_raw.root", thisRunNumber)))->Get("tree");

    const double d_eta[8] = {1.7, 1.2, 1, 1, 1, 1, 1.2, 1.7};
    const double eta_cuts[9] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
    const double harr_scales[8] = {1, 1, 1, 1, 1, 1, 1, 1};   // rescaling factors so the histograms don't overlap

    const double* xbins = logspace(0, 1.6, numbins);
    // Create an array of 16 histograms, one for each rapidity region and one for x_p, x_a. 
    TH1D* harr[numhists];
    for (int i = 0; i < numhists/2; i++) {
        harr[i] = new TH1D(Form("%ieta%i", thisRunNumber, i), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d#sigma^{2}/d#it{x}_{p} dy #left[pb#right]", eta_cuts[i], eta_cuts[i+1]), numbins, xbins);
        harr[i]->Sumw2();
    }
    for (int i = numhists/2; i < numhists; i++) {
        harr[i] = new TH1D(Form("%ieta%i", thisRunNumber, i), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d#sigma^{2}/d#it{x}_{a} dy #left[pb#right]", eta_cuts[i%(numhists/2)], eta_cuts[(i%(numhists/2))+1]), numbins, xbins);
        harr[i]->Sumw2();
    }

//    cout << numtrigs << endl;

    // Create arrays to store trigger values for each event
    bool m_trig_bool[numtrigs];
    float m_trig_prescale[numtrigs];

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
    for (Trigger* trig : trigger_vec) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &m_trig_bool[trig->index]);
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &m_trig_prescale[trig->index]);
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();
        
    double leading_jpt, leading_jeta, jeta0, jpt0, jeta1, jpt1, xp, xa, extra_jpt_sum;
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
            jpt1 = (double)j_pt[1];
            
            jeta0 = (double)j_eta[0];
            jeta1 = (double)j_eta[1];

            leading_jpt = std::max(jpt0, jpt1);
            if (jpt0 == leading_jpt) leading_jeta = jeta0;
            else leading_jeta = jeta1;

            xp = get_xp(jpt0, jpt1, jeta0, jeta1, periodA);
            xa = get_xa(jpt0, jpt1, jeta0, jeta1, periodA);

            if (0 < leading_jeta && leading_jeta < 1) {
                for (Trigger* trig : triggers_p0eta100) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0; // Only take the event if triggered and trigger was not disabled
                    if (takeEvent && trig->min_pt <= leading_jpt && leading_jpt < trig->max_pt) {
                        if (!periodA) {
                            harr[(numhists/4)]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)]->Fill(xa, m_trig_prescale[trig->index]);
                        } else {
                            harr[(numhists/4)-1]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)-1]->Fill(xa, m_trig_prescale[trig->index]);
                        }
                        //harr[4]->Fill(xp, m_trig_prescale[trig->index]);
                        //harr[12]->Fill(xa, m_trig_prescale[trig->index]);
                        break; // Break to ensure that only one trigger allows the event to be recorded
                    }
                }
            }
            else if (-1 < leading_jeta && leading_jeta < 0) {
                for (Trigger* trig : triggers_n100eta0) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= leading_jpt && leading_jpt < trig->max_pt) {
                        if (!periodA) {
                            harr[(numhists/4)-1]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)-1]->Fill(xa, m_trig_prescale[trig->index]);
                        } else {
                            harr[(numhists/4)]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)]->Fill(xa, m_trig_prescale[trig->index]);
                        }
                        //harr[3]->Fill(xp, m_trig_prescale[trig->index]);
                        //harr[11]->Fill(xa, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (1 < leading_jeta && leading_jeta < 2) {
                for (Trigger* trig : triggers_p100eta200) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= leading_jpt && leading_jpt < trig->max_pt) {
                        if (!periodA) {
                            harr[(numhists/4)+1]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)+1]->Fill(xa, m_trig_prescale[trig->index]);
                        } else {
                            harr[(numhists/4)-2]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)-2]->Fill(xa, m_trig_prescale[trig->index]);
                        }
                        //harr[5]->Fill(xp, m_trig_prescale[trig->index]);
                        //harr[13]->Fill(xa, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (-2 < leading_jeta && leading_jeta < -1) {
                for (Trigger* trig : triggers_n200eta100) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= leading_jpt && leading_jpt < trig->max_pt) {
                        if (!periodA) {
                            harr[(numhists/4)-2]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)-2]->Fill(xa, m_trig_prescale[trig->index]);
                        } else {
                            harr[(numhists/4)+1]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)+1]->Fill(xa, m_trig_prescale[trig->index]);
                        }
                        //harr[2]->Fill(xp, m_trig_prescale[trig->index]);
                        //harr[10]->Fill(xa, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (2 < leading_jeta && leading_jeta < 3.2) {
                for (Trigger* trig : triggers_p200eta320) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= leading_jpt && leading_jpt < trig->max_pt) {
                        if (!periodA) {
                            harr[(numhists/4)+2]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)+2]->Fill(xa, m_trig_prescale[trig->index]);
                        } else {
                            harr[(numhists/4)-3]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)-3]->Fill(xa, m_trig_prescale[trig->index]);
                        }
                        //harr[6]->Fill(xp, m_trig_prescale[trig->index]);
                        //harr[14]->Fill(xa, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (-3.2 < leading_jeta && leading_jeta < -2) {
                for (Trigger* trig : triggers_n320eta200) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= leading_jpt && leading_jpt < trig->max_pt) {
                        if (!periodA) {
                            harr[(numhists/4)-3]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)-3]->Fill(xa, m_trig_prescale[trig->index]);
                        } else {
                            harr[(numhists/4)+2]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)+2]->Fill(xa, m_trig_prescale[trig->index]);
                        }
                        //harr[1]->Fill(xp, m_trig_prescale[trig->index]);
                        //harr[9]->Fill(xa, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (3.2 < leading_jeta && leading_jeta < 4.9) {
                for (Trigger* trig : triggers_p320eta490) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= leading_jpt && leading_jpt < trig->max_pt) {
                        if (!periodA) {
                            harr[(numhists/4)+3]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)+3]->Fill(xa, m_trig_prescale[trig->index]);
                        } else {
                            harr[(numhists/4)-4]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)-4]->Fill(xa, m_trig_prescale[trig->index]);
                        }
                        //harr[7]->Fill(xp, m_trig_prescale[trig->index]);
                        //harr[15]->Fill(xa, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (-4.9 < leading_jeta && leading_jeta < -3.2) {
                for (Trigger* trig : triggers_n490eta320) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= leading_jpt && leading_jpt < trig->max_pt) {
                        if (!periodA) {
                            harr[(numhists/4)-4]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)-4]->Fill(xa, m_trig_prescale[trig->index]);
                        } else {
                            harr[(numhists/4)+3]->Fill(xp, m_trig_prescale[trig->index]);
                            harr[(3*numhists/4)+3]->Fill(xa, m_trig_prescale[trig->index]);
                        }
                        //harr[0]->Fill(xp, m_trig_prescale[trig->index]);
                        //harr[8]->Fill(xa, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
        }
    }
    // Save to root file
    TFile* output = new TFile(Form("./xdata/run_%i.root", thisRunNumber), "RECREATE");
    for (int i = 0; i < numhists; i++) {
        harr[i]->Scale(harr_scales[i%(numhists/2)]/(A * luminosity * d_eta[i%(numhists/2)]), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        harr[i]->Write();
    }
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    output->Close();
}
