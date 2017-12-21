#include "../triggerUtil.C"

void triggers_pt_counts(const int thisRunNumber, // Run number identifier.
                       double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;

    initialize(thisRunNumber, false);
    vector<Trigger*> triggerSubList(0);
    for (Trigger* trig : trigger_vec) {
        if (useIonTrigs != trig->iontrigger) continue;
        triggerSubList.push_back(trig);
    }

    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)
    const int numhists = numtrigs * numetabins;

    TTree* tree = (TTree*)(new TFile(Form("%srun_%i_raw.root", dataPath.c_str(), thisRunNumber)))->Get("tree");

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
    for (Trigger* trig : triggerSubList) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &m_trig_bool[trig->index]);
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &m_trig_prescale[trig->index]);
    }

    TH1F* harr[numhists];
    int pbin, ebin, index;
    for (Trigger* trig : trigger_vec) {
        index = trig->index;
        for (ebin = 0; ebin < numetabins; ebin++) {
            TString histname = Form("trig_pt_counts_run%i_trig%i_ebin%i", thisRunNumber, index, ebin);
            harr[index + ebin*numtrigs] = new TH1F(histname, ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
            harr[index + ebin*numtrigs]->Sumw2(); // instruct each histogram to propagate errors
        }
    }


    // Iterate over each event
    const int numentries = tree->GetEntries();

    double jpt, jeta;
    TVectorD numtrigfirings(numtrigs * numpbins * numetabins);
    for (int n = 0; n < numtrigs*numpbins*numetabins; n++) {
        numtrigfirings[n] = 0;
    }

    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses

        for (Trigger* trig : triggerSubList) {
//            if (useIonTrigs != trig->iontrigger) continue; // Only consider ion triggers in p. A and non-ion triggers in p. B
            index = trig->index;
            if (!m_trig_bool[index] || m_trig_prescale[index] <= 0) continue; // if the trigger wasn't fired (or was disabled in some way) just continue.
        
            for (pbin = 0; pbin < numpbins; pbin++) {
                for (ebin = 0; ebin < numetabins; ebin++) {
                    for (int j = 0; j < njet; j++) {
                        jpt = (double)j_pt[j];
                        jeta = (double)j_eta[j];
                        if (etabins[ebin] <= jeta && jeta < etabins[ebin+1] && pbins[pbin] <= jpt && jpt < pbins[pbin+1] && trig->lower_eta <= jeta && jeta < trig->upper_eta && trig->min_pt[ebin] <= jpt) {
                            numtrigfirings[index + (pbin + ebin*numpbins)*numtrigs]++;
                            break;
                        }
                    }
                }
            }

            for (int j = 0; j < njet; j++) {
                jpt = (double)j_pt[j];
                jeta = (double)j_eta[j];

                ebin = 0;
                while (etabins[ebin] <= jeta) ebin++;
                ebin--;
                if (ebin == -1 || ebin >= numetabins) continue;

                if (trig->min_pt[ebin] <= jpt && trig->lower_eta <= jeta && jeta < trig->upper_eta) {
                    harr[index + ebin*numtrigs]->Fill(jpt, m_trig_prescale[index]);
                }
            } 
        }       
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("%srun_%i.root", ptPath.c_str(), thisRunNumber), "RECREATE");
    for (Trigger* trig : trigger_vec) {
        index = trig->index;
        for (ebin = 0; ebin < numetabins; ebin++) {
            harr[index + ebin*numtrigs]->Scale(1/A); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
            harr[index + ebin*numtrigs]->Write();
        }
    }
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");

    TVectorD run_vec(3);
    run_vec[0] = thisRunNumber;
    run_vec[1] = numetabins;
    run_vec[2] = numtrigs;
    run_vec.Write("run_vec");

    numtrigfirings.Write("trig_fire_vec");
        
    output->Close();
}


