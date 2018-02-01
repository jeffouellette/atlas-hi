#include "../triggerUtil.C"

void IdealPtAnalysis(const int thisRunNumber, // Run number identifier.
                       double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;

    initialize(thisRunNumber, true, true);
    
    /**** Generate list of physics triggers ****/
    vector<Trigger*> triggerSubList(0);
    for (Trigger* trig : triggerVec) {
        if (trig->lowerRunNumber <= thisRunNumber && thisRunNumber < trig->upperRunNumber && trig->name != minbiasTriggerName) triggerSubList.push_back(trig);
    }
    if (debugStatements) {
        cout << "Status: In triggers_pt_counts.C (16): Processing run " << thisRunNumber << " with triggers:" << endl;
        for (Trigger* trig : triggerSubList) {
            cout << "\t" << trig->name << endl;
        }
    }
    /**** End generate list of physics triggers ****/

    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)
    const int numhists = numtrigs * numetabins;

    /**** Find the relevant TTree for this run ****/
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
                    if (debugStatements) cout << "Status: In triggers_pt_counts.C (39): Found " << fname.Data() << endl; 
                    if (fname.Contains(to_string(thisRunNumber))) {
                        tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
                        break;
                    }
                }
            }
        }
    }
    if (tree == NULL) {
        cout << "Error: In triggers_pt_counts.C (49): TTree not obtained for given run number. Quitting." << endl;
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
            if (debugStatements) cout << "Status: In triggers_pt_counts.C (62): Tree contains branch \"" << branchName.Data() << "\"" << endl;
            interestingBranch = false;
            for (string s : interestingBranchNames) {
                interestingBranch = interestingBranch || (branchName.Data() == s);
            }
            if (!interestingBranch) {
                for (Trigger* trig : triggerSubList) {
                    if (branchName == trig->name) {
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


    /**** Set branching addresses ****/
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
    //    tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &(trig->m_trig_prescale));
    }
    /**** End set branch addresses ****/


    /**** Histogram initialization ****/
    TH1F* histArr[numhists];
    int pbin, etabin, index;
    for (etabin = 0; etabin < numetabins; etabin++) {
        TString histName = Form("trig_pt_counts_run%i_etabin%i", thisRunNumber, etabin);
        histArr[etabin] = new TH1F(histName, ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T}dy #left[nb (GeV/#it{c})^{-1}#right]", numpbins, pbins);
        histArr[etabin]->Sumw2(); // instruct each histogram to propagate errors
    }


    /**** Iterate over each event ****/
    const int numentries = tree->GetEntries();

    double jpt, jeta, eff;
    Trigger* bestTrigger = NULL;
    for (long long entry = 0; entry < numentries; entry++) {
        tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses

        for (int j = 0; j < njet; j++) {
            jpt = (double)j_pt[j];
            jeta = (double)j_eta[j];

            etabin = getEtabin(jeta);
            pbin = getPbin(jpt);

            if (pbin < 0 || etabin < 0 || pbin > numpbins || etabin > numetabins) continue; // this checks that the jets fall within the pt, eta bins

            bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins];
            eff = kinematicEfficiencyVec[pbin + etabin*numpbins]; 
            if (bestTrigger == NULL || eff == 0) continue; // make sure we're not trying to look at a null trigger
            if (bestTrigger->m_trig_bool) histArr[etabin]->Fill(jpt, 1./eff);
        }
    }
    /**** End event iteration ****/


    /**** Write output histograms to a root file ****/
    TFile* output = new TFile(Form("%srun_%i.root", ptPath.c_str(), thisRunNumber), "RECREATE");
    for (etabin = 0; etabin < numetabins; etabin++) {
        histArr[etabin]->Scale(1/A); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        histArr[etabin]->Write();
    }
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write(Form("lum_vec_%i", thisRunNumber));

    TVectorD run_vec(4);
    run_vec[0] = thisRunNumber;
    run_vec[1] = numetabins;
    run_vec[2] = numtrigs;
    run_vec[3] = numpbins;
    run_vec.Write(Form("run_vec_%i", thisRunNumber));

    output->Close();
    /**** End write output ****/

    if (debugStatements) cout << "Status: In triggers_pt_counts.C (163): Finished calculating pt spectrum for run " << thisRunNumber << endl;
    return;
}
