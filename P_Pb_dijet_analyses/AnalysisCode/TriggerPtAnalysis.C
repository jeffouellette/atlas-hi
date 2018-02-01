#include "../triggerUtil.C"

// IDEA: plot the number of times each trigger fired for a particular run number.
// For overlapping triggers, this will improve statistics by choosing the trigger that fires more often.

void TriggerPtAnalysis(const int thisRunNumber, // Run number identifier.
             double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;

    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    initialize(thisRunNumber, false, false);


    /**** Generate list of physics triggers ****/
    vector<Trigger*> triggerSubList(0);
    for (Trigger* trig : triggerVec) {
        if (trig->lowerRunNumber <= thisRunNumber && thisRunNumber < trig->upperRunNumber && trig->name != minbiasTriggerName) triggerSubList.push_back(trig);
    }
    if (debugStatements) {
        cout << "Status: In triggers.C (22): Processing run " << thisRunNumber << " with triggers:" << endl;
        for (Trigger* trig : triggerSubList) {
            cout << "\t" << trig->name << endl;
        }
    }
    double* triggerEfficiencies = getTriggerEfficiencies();
    const int numtrigs_sublist = triggerSubList.size();
    const double* trigbins = linspace(0, numtrigs_sublist+1, numtrigs_sublist+1);
    /**** End generate list of physics triggers ****/


    /**** Find the relevant TTree for this run ****/
    TFile* thisFile = NULL;
    TTree* tree = NULL;
    {
        TSystemDirectory dir(dataPath.c_str(), dataPath.c_str());
        TList* sysfiles = dir.GetListOfFiles();
        if (sysfiles) {
            TSystemFile* sysfile;
            TString fname;
            TIter next(sysfiles);

            while ((sysfile = (TSystemFile*)next())) {
                fname = sysfile->GetName();
                if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
                    if (debugStatements) cout << "Status: In triggers_pt_counts.C (47): Found " << fname.Data() << endl; 
                    if (fname.Contains(to_string(thisRunNumber))) {
                        thisFile = new TFile(dataPath+fname, "READ");
                        tree = (TTree*)thisFile->Get("tree");
                        break;
                    }
                }
            }
        }
    }
    if (tree == NULL) {
        cout << "Error: In triggers_pt_counts.C (58): TTree not obtained for given run number. Quitting." << endl;
        return;
    }
    /**** End find TTree ****/

    /**** Disable loading of unimportant branch values - speeds up entry retrieval ****/
    {
        vector<string> interestingBranchNames = {"njet", "j_pt", "j_eta"};
        TObjArray* branches = (TObjArray*)(tree->GetListOfBranches());
        bool interestingBranch;
        for (TObject* obj : *branches) {
            TString branchName = (TString)obj->GetName();
            if (debugStatements) cout << "Status: In triggers_pt_counts.C (63): Tree contains branch \"" << branchName.Data() << "\"" << endl;
            interestingBranch = false;
            for (string s : interestingBranchNames) {
                interestingBranch = interestingBranch || (branchName.Data() == s);
            }
            if (!interestingBranch) {
                for (Trigger* trig : triggerSubList) {
                    if (branchName == trig->name || branchName == trig->name + "_prescale") {
                        interestingBranch = true;
                        break;
                    }
                }
            }
            if (!interestingBranch) {
                tree->SetBranchStatus(branchName, 0);
            }
        }
    }
    /**** End disable unimportant branches ****/

    /*if (thisRunNumber == 313063) {
        cout << Form("Numtrigs: %i", numtrigs) << endl;
        for (Trigger* trig : triggerVec) {
            int i = trig->index;
            cout << Form("Trigger %s is bin %i centered at %.1f", trig->name.c_str(), i, 0.5*(trigbins[i+1]+trigbins[i])) << endl;
        }
        for (int pbin = 0; pbin < numpbins; pbin++) {
            cout << Form("Momentum range %.0f...%.0f is bin %i centered at %.1f", pbins[pbin], pbins[pbin+1], pbin, 0.5*(ybins_pt[pbin+1]+ybins_pt[pbin])) << endl;
        }
    }*/


    TH1D* integratedCountsHistArr[numtrigs_sublist];
    TH1D* countsHistArr[numtrigs_sublist*numetabins];
    for (int t = 0; t < numtrigs_sublist; t++) {
        Trigger* trig = triggerSubList[t];
        // integrated counts plot for each trigger -- I don't think I actually need to fill these since I can just integrate the pt spectrum defined below instead
        TString histName = Form("integratedCounts_%s_run%i", trig->name.c_str(), thisRunNumber);
        integratedCountsHistArr[t] = new TH1D(histName, ";Trigger;Counts / L_{int} #left[nb#right]", numtrigs_sublist, trigbins);
        // pt spectrum for each trigger
        for (int etabin = 0; etabin < numetabins; etabin++) {
            histName = Form("counts_%s_etabin_%i_run%i", trig->name.c_str(), etabin, thisRunNumber);
            countsHistArr[t + etabin*numtrigs_sublist] = new TH1D(histName, ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T}dy #left[nb GeV^{-1}#right]", numpbins, pbins);
            countsHistArr[t + etabin*numtrigs_sublist]->Sumw2();
        }
    }

    // Create branching addresses:  
    // Create arrays to store trigger values for each event
    //bool m_trig_bool[numtrigs];   // stores whether trigger was triggered
    //float m_trig_prescale[numtrigs];      // stores the prescaling factor for the trigger
    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    int njet = 0;
    // Set branch addresses
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("njet", &njet);
    for (Trigger* trig : triggerSubList) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->m_trig_bool));
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &(trig->m_trig_prescale));
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();
    Trigger* trig;
    int index, pbin, etabin;
    double jpt, jeta, eff;
    for (long long event = 0; event < numentries; event++) {
        tree->GetEntry(event); // stores trigger values and data in the designated branch addresses

        for (int t = 0; t < numtrigs_sublist; t++) {
            trig = triggerSubList[t];
            index = trig->index;
            if (!(trig->m_trig_bool) || trig->m_trig_prescale <= 0) continue;
            
            for (int j = 0; j < njet; j++) {
                jpt = (double)j_pt[j];
                jeta = (double)j_eta[j];
                if (jpt < pbins[0] || jeta < etabins[0] || jeta > etabins[numetabins]) continue;

                pbin = 0;
                while (pbins[pbin] <= jpt) pbin++;
                pbin--;
                etabin = 0;
                while (etabins[etabin] <= jeta) etabin++;
                etabin--;

                if (trig->lower_eta <= jeta && jeta < trig->upper_eta && trig->min_pt <= jpt) {
                    eff = triggerEfficiencies[index + pbin*numtrigs];
                    if (eff > 0.) countsHistArr[t + etabin*numtrigs_sublist]->Fill(jpt, 1./eff);
                }
            }
        }
    }
    for (int t = 0; t < numtrigs_sublist; t++) {
        double integral = 0;
        for (int etabin = 0; etabin < numetabins; etabin++) integral += countsHistArr[t + etabin*numtrigs_sublist]->Integral();
        integratedCountsHistArr[t]->SetBinContent(t+1, integral);
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("%srun_%i.root", trigPath.c_str(), thisRunNumber), "RECREATE");
    for (int t = 0; t < numtrigs_sublist; t++) {
        Trigger* trig = triggerSubList[t];
        for (int etabin = 0; etabin < numetabins; etabin++) countsHistArr[t + etabin*numtrigs_sublist]->Write();
        integratedCountsHistArr[t]->Write();
    }
    
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write(Form("lum_vec_%i", thisRunNumber));

    TVector run_info(3);
    run_info[0] = thisRunNumber;
    run_info[1] = numtrigs;
    run_info[2] = numentries;
    run_info.Write(Form("run_info_%i", thisRunNumber));

    output->Close();

    if (debugStatements) cout << "Status: In triggers.C (193): Finished run " << thisRunNumber << endl;

    delete[] triggerEfficiencies;
    return;
}
