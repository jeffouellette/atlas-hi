#include "../triggerUtil.C"

// IDEA: plot the number of times each trigger fired for a particular run number.
// For overlapping triggers, this will improve statistics by choosing the trigger that fires more often.

void triggers(const int thisRunNumber, // Run number identifier.
             double luminosity, // Integrated luminosity for this run. Presumed constant over the run period.
             bool periodA)
{

    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    initialize(thisRunNumber, false);

    TTree* tree = (TTree*)(new TFile(Form("./rundata/run_%i_raw.root", thisRunNumber)))->Get("tree");

    const double* trigbins = linspace(-0.5, numtrigs+0.5, numtrigs+1);
    const double* ybins = linspace(-0.5, numpbins+0.5, numpbins+1);
    const double* ebins = linspace(-0.5, numetabins+0.5, numetabins+1);

    if (thisRunNumber == 313063) {
        cout << Form("Numtrigs: %i", numtrigs) << endl;
        for (Trigger* trig : trigger_vec) {
            int i = trig->index;
            cout << Form("Trigger %s is bin %i centered at %.3f", trig->name.c_str(), i, 0.5*(trigbins[i+1]+trigbins[i])) << endl;
        }
        for (int pbin = 0; pbin < numpbins; pbin++) {
            cout << Form("Momentum range %.0f...%.0f is bin %i centered at %.0f", pbins[pbin], pbins[pbin+1], pbin, 0.5*(ybins[pbin+1]+ybins[pbin])) << endl;
        }
    }


    TH1D* harr[numtrigs];
    for (Trigger* trig : trigger_vec) {
        int i = trig->index;
        harr[i] = new TH1D(Form("trig%i", i), ";Trigger;Number of times triggered / integrated luminosity #left[pb^{-1}#right]", numtrigs, trigbins);
        harr[i]->Sumw2();
    }

    TH2D* h2_meta = new TH2D("pt_trig", ";Trigger bin;p^{jet}_{T} bin #left[GeV #it{c}^{-1}#right];Counts / dy Luminosity", numtrigs, trigbins, numpbins, ybins);

    // Create an additional data storage histogram for binning in eta as well.
    TH3D* h_meta = new TH3D("eta_pt_trig", "", numtrigs, trigbins, numpbins, ybins, numetabins, ebins);

    // Create branching addresses:  
    // Create arrays to store trigger values for each event
    bool m_trig_bool[numtrigs];   // stores whether trigger was triggered
    float m_trig_prescale[numtrigs];      // stores the prescaling factor for the trigger
    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    int njet = 0;

    // Set branch addresses
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("njet", &njet);
    for (Trigger* trig : trigger_vec) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &m_trig_bool[trig->index]);
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &m_trig_prescale[trig->index]);
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();
    Trigger* trig;
    bool takeEvent;
    int numticks = 0;
    int index, pbin, ebin;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses

        for (Trigger* trig : trigger_vec) { // Just find if the trigger was fired
            index = trig->index;
            if (!m_trig_bool[index] || m_trig_prescale <= 0) continue;
        
            numticks++;
            harr[index]->Fill(index);
            
            for (int j = 0; j < njet; j++) {
                pbin = 0;
                while (pbins[pbin] <= j_pt[j]) pbin++;
                pbin--;
                ebin = 0;
                while (etabins[ebin] <= j_eta[j]) ebin++;
                ebin--;

                if (pbin == -1 || ebin == -1 || pbin >= numpbins || ebin >= numetabins) continue;

                if (trig->lower_eta <= etabins[ebin] && etabins[ebin+1] <= trig->upper_eta && pbins[pbin] >= trig->min_pt) {
                    h2_meta->Fill((double)index, (double)pbin);
                    h_meta->Fill((double)index, (double)pbin, (double)ebin);
                    if (thisRunNumber == 313063 && j_pt[j] > 400) cout << Form("Found high-pt jet with pt=%.1f. Trigger %i=%i with prescale %.1f (filling h_meta files, content of h_meta is now %g)", j_pt[j], index, m_trig_bool[index], m_trig_prescale[index], h_meta->GetBinContent(index+1, pbin+1, ebin+1)) << endl;
                }
            }
        }
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("./trig_data/run_%i.root", thisRunNumber), "RECREATE");
    for (Trigger* trig : trigger_vec) {
        int i = trig->index;
        harr[i]->Scale(1/((trig->upper_eta-trig->lower_eta)*luminosity));
        harr[i]->Write();
    }
    
    h2_meta->Write();

    h_meta->Write();

    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");

    TVector run_info(3);
    run_info[0] = thisRunNumber;
    run_info[1] = numtrigs;
    run_info[2] = numentries;
    run_info.Write("run_info");

    output->Close();

    cout << Form("Finished run %i with %i triggers fired", thisRunNumber, numticks) << endl;
}
