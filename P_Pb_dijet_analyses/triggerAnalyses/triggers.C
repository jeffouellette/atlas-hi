#include "../triggerUtil.C"

// IDEA: plot the number of times each trigger fired for a particular run number.
// For overlapping triggers, this will improve statistics by choosing the trigger that fires more often.

void triggers(const int thisRunNumber, // Run number identifier.
             double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;

    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    initialize(thisRunNumber, false, true);

    //TTree* tree = (TTree*)(new TFile(Form("%srun_%i_raw.root", dataPath.c_str(), thisRunNumber)))->Get("tree");
    TTree* tree = NULL;
    TSystemDirectory dir(dataPath.c_str(), dataPath.c_str());
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                TFile* thisfile = new TFile(dataPath+fname, "READ");
                TTree* thistree = (TTree*)thisfile->Get("tree");
                int rn;
                thistree->SetBranchAddress("runNumber", &rn);
                thistree->GetEvent(0);
                if (rn == thisRunNumber) {
                    tree = thistree;
                    break;
                }
                thisfile->Close();
                thisfile->Delete();
            }
        }
    }
    if (tree == NULL) {
        cout << "TTree not obtained for given run number. Quitting." << endl;
        return;
    }

    const double* trigbins = linspace(0, numtrigs+1, numtrigs+1);
    const double* ybins_pt = linspace(0, numpbins+1, numpbins+1);
    const double* ybins_eta = linspace(0, numetabins+1, numetabins+1);

    if (thisRunNumber == 313063) {
        cout << Form("Numtrigs: %i", numtrigs) << endl;
        for (Trigger* trig : trigger_vec) {
            int i = trig->index;
            cout << Form("Trigger %s is bin %i centered at %.1f", trig->name.c_str(), i, 0.5*(trigbins[i+1]+trigbins[i])) << endl;
        }
        for (int pbin = 0; pbin < numpbins; pbin++) {
            cout << Form("Momentum range %.0f...%.0f is bin %i centered at %.1f", pbins[pbin], pbins[pbin+1], pbin, 0.5*(ybins_pt[pbin+1]+ybins_pt[pbin])) << endl;
        }
    }


    TH1D* harr[numtrigs];
    for (Trigger* trig : trigger_vec) {
        int i = trig->index;
        harr[i] = new TH1D(Form("trig%i", i), ";Trigger;N^{trig}_{counts} / #it{L}_{int} #left[pb#right]", numtrigs, trigbins);
        harr[i]->Sumw2();
    }

    TH2D* hist_pt = new TH2D("pt_trig", ";Trigger;#it{p}^{jet}_{T} #left[GeV/#it{c}#right];dN^{trig}_{counts} / #it{L}d#it{p}_{T} #left[GeV^{-1} pb#right]", numtrigs, trigbins, numpbins, ybins_pt);
    TH2D* hist_eta = new TH2D("eta_trig", ";Trigger;#it{#eta}^{jet};dN^{trig}_{counts} / #it{L}d#it{#eta} #left[pb#right]", numtrigs, trigbins, numetabins, ybins_eta);

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
        if (useIonTrigs != trig->iontrigger) continue;
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &m_trig_bool[trig->index]);
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &m_trig_prescale[trig->index]);
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();
    Trigger* trig;
    bool takeEvent;
    int numticks = 0;
    int index, pbin, ebin;
    double jpt, jeta;
    for (long long i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses

        for (Trigger* trig : trigger_vec) { // Just find if the trigger was fired
            if (periodA != trig->iontrigger) continue;
            index = trig->index;
            if (!m_trig_bool[index] || m_trig_prescale <= 0) continue;
        
            numticks++;
            harr[index]->Fill(((double)index)+0.5);
            
            for (int j = 0; j < njet; j++) {
                jpt = (double)j_pt[j];
                jeta = (double)j_eta[j];

                pbin = 0;
                while (pbins[pbin] <= jpt) pbin++;
                pbin--;
                ebin = 0;
                while (etabins[ebin] <= jeta) ebin++;
                ebin--;

                if (pbin == -1 || pbin >= numpbins || ebin == -1 || ebin >= numetabins) continue;

                if (trig->lower_eta <= jeta && jeta < trig->upper_eta && trig->min_pt[ebin] <= pbins[pbin]) {
                    hist_pt->Fill(((double)index)+0.5, ((double)pbin)+0.5);
                    hist_eta->Fill(((double)index)+0.5, ((double)ebin)+0.5);
                }
            }
        }
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("%srun_%i.root", trigPath.c_str(), thisRunNumber), "RECREATE");
    for (Trigger* trig : trigger_vec) {
        int i = trig->index;
        harr[i]->Write();
    }
    
    hist_pt->Write();
    hist_eta->Write();

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
    return;
}
