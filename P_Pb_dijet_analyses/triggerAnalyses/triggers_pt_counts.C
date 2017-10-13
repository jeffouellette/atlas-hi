#include "triggerUtil.C"

void triggers_pt_counts(const int thisRunNumber, // Run number identifier.
                       double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{

    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    initialize(thisRunNumber, false);
    const int numhists = numtrigs;

    TTree* tree = (TTree*)(new TFile(Form("./rundata/run_%i_raw.root", thisRunNumber)))->Get("tree");

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

    TH1D* harr[numhists];
    int pbin, ebin, index;
    for (Trigger* trig : trigger_vec) {
        index = trig->index;
        TString histname = Form("trig_pt_counts_run%i_trig%i", thisRunNumber, index);
        harr[index] = new TH1D(histname, ";#it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
        harr[index]->Sumw2(); // instruct each histogram to propagate errors
    }


    // Iterate over each event
    const int numentries = tree->GetEntries();

    double jpt, jeta;

    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses

        for (Trigger* trig : trigger_vec) {
            index = trig->index;
            if (!m_trig_bool[index] || m_trig_prescale[index] <= 0) continue;

            for (int j = 0; j < njet; j++) {
                jpt = (double)j_pt[j];
                jeta = (double)j_eta[j];
                if (jpt >= trig->min_pt && trig->lower_eta <= jeta && jeta < trig->upper_eta) harr[index]->Fill(jpt, m_trig_prescale[index]);
            } 

        }       
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("./pt_data/trig_bin/run_%i.root", thisRunNumber), "RECREATE");
    for (Trigger* trig : trigger_vec) {
        index = trig->index;
        harr[index]->Scale(1/(A * (trig->upper_eta - trig->lower_eta)), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        harr[index]->Write();
    }
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    output->Close();
}


