#include "../triggerUtil.C"

void bootstrap(Trigger* trig, TH1F** harr) {
    if (trig->referenceTrigger->index == trig->index) {
        for (int pbin = 0; pbin < numpbins; pbin++) {
            harr[trig->referenceTrigger->index]->SetBinContent(pbin+1, 1);
            harr[trig->referenceTrigger->index]->SetBinError(pbin+1, 0);
        }
        trig->isBootstrapped = true;
    }
    if (trig->isBootstrapped) return;

    int index = trig->index;
    int referenceIndex = trig->referenceTrigger->index;
    if (!(trig->referenceTrigger->isBootstrapped)) {
        if (debugStatements) cout << "Bootstrapping trigger " << trig->referenceTrigger->name << " for " << trig->name << endl;
        bootstrap(trig->referenceTrigger, harr);
    }
    harr[index]->Multiply(harr[referenceIndex]);
    trig->isBootstrapped = true;
    return;
}


void trigger_efficiencies(const int thisRunNumber, // Run number identifier.
                       double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;

    initialize(thisRunNumber, false, false);

    vector<Trigger*> triggerSubList(0);
    for (Trigger* trig : trigger_vec) {
        if (!trig->disabled) triggerSubList.push_back(trig);
    }
    if (debugStatements) {
        cout << "Processing run " << thisRunNumber << " with triggers:" << endl;
        for (Trigger* trig : trigger_vec) {
            cout << trig->name << endl;
        }
    }

    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)
    const int numhists = numtrigs;

    TTree* tree = NULL;
    TSystemDirectory dir(dataPath.c_str(), dataPath.c_str());
    TList* sysfiles = dir.GetListOfFiles();
    if (sysfiles) {
        TSystemFile *sysfile;
        TString fname;
        TIter next(sysfiles);

        int* rn;
        while ((sysfile=(TSystemFile*)next())) {
            fname = sysfile->GetName();
            if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
                if (debugStatements) cout << "Found " << fname.Data() << endl; 
                if (fname.Contains(to_string(thisRunNumber))) {
                    tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
                    break;
                }
            }
        }
    }
    if (tree == NULL) {
        cout << "TTree not obtained for given run number. Quitting." << endl;
        return;
    }

    // Disable loading of unimportant branch values
    vector<string> interestingBranchNames = {"njet", "j_pt", "j_eta", "hlt_ion_j_pt", "hlt_ion_j_eta", "hlt_ion_njet", "hlt_j_pt", "hlt_j_eta", "hlt_njet"};
    TObjArray* branches = (TObjArray*)(tree->GetListOfBranches());
    for (TObject* obj : *branches) {
        TString branchName = (TString)obj->GetName();
        if (debugStatements) cout << "Tree contains branch \"" << branchName.Data() << "\"" << endl;
        bool interestingBranch;
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

    // Create branching addresses:  
    // Create arrays to store trigger values for each event
    //bool m_trig_bool[numtrigs];   // stores whether trigger was triggered
    //float m_trig_prescale[numtrigs];      // stores the prescaling factor for the trigger
    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    int njet = 0;

    float hlt_j_pt[60] = {};
    float hlt_j_eta[60] = {};
    int hlt_njet = 0;

    // Set branch addresses
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("njet", &njet);
    string hlt_j_pt_name, hlt_j_eta_name, hlt_njet_name;
    if (thisRunNumber <= 313603) {
/*        hlt_j_pt_name = "hlt_ion_j_pt";
        hlt_j_eta_name = "hlt_ion_j_eta";
        hlt_njet_name = "hlt_ion_njet";
    } else {
        hlt_j_pt_name = "hlt_j_pt";
        hlt_j_eta_name = "hlt_j_eta";
        hlt_njet_name = "hlt_njet";
    }
    tree->SetBranchAddress(hlt_j_pt_name.c_str(), hlt_j_pt);
    tree->SetBranchAddress(hlt_j_eta_name.c_str(), hlt_j_eta);
    tree->SetBranchAddress(hlt_njet_name.c_str(), &hlt_njet);*/

        tree->SetBranchAddress("hlt_ion_j_pt", hlt_j_pt);
        tree->SetBranchAddress("hlt_ion_j_eta", hlt_j_eta);
        tree->SetBranchAddress("hlt_ion_njet", &hlt_njet);
    } else {
        tree->SetBranchAddress("hlt_j_pt", hlt_j_pt);
        tree->SetBranchAddress("hlt_j_eta", hlt_j_eta);
        tree->SetBranchAddress("hlt_njet", &hlt_njet);
    }
    for (Trigger* trig : triggerSubList) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->m_trig_bool));
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &(trig->m_trig_prescale));
    }

    int pbin, ebin, index, referenceIndex;
    TH1F* harr[2*numtrigs];
    for (Trigger* trig : triggerSubList) {
        index = trig->index;
        // standard trigger firings
        TString histname = Form("%s_efficiency_run%i", trig->name.c_str(), thisRunNumber);
        harr[index] = new TH1F(histname, ";#it{p}_{T}^{jet} #left[GeV#right];#epsilon", numpbins, pbins);
        harr[index]->Sumw2(); // instruct each histogram to propagate errors
        // reference trigger firings
        histname = Form("%s_reference_run%i", trig->name.c_str(), thisRunNumber);
        harr[index+numtrigs] = new TH1F(histname, ";#it{p}_{T}^{jet} #left[GeV#right];#epsilon", numpbins, pbins);
        harr[index+numtrigs]->Sumw2();
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();
    double max_j_pt, max_hlt_j_pt, p1, p2, p_adj;
    if (debugStatements) cout << "Looping over " << numentries << " events in run " << thisRunNumber << endl;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses

        /**          BOOTSTRAPPING METHOD          **/
        /** Method works by calculating lowest pt  **/
        /** jet trigger fire rate relative to min  **/
        /** bias trigger, then 2nd lowest pt       **/
        /** relative to lowest pt multiplied by    **/
        /** lowest trigger efficiency, ad          **/
        /** infinitum.                             **/

        for (Trigger* trig : triggerSubList) {
            index = trig->index;
            referenceIndex = trig->referenceTrigger->index;
            max_j_pt = 0; // find leading (reconstructed) jet
            for (int j = 0; j < njet; j++) {
                if (j_pt[j] > max_j_pt && trig->lower_eta <= j_eta[j] && j_eta[j] <= trig->upper_eta) {
                    max_j_pt = (double)j_pt[j];
                }
            }
            max_hlt_j_pt = 0; // find leading jet as seen by the HLT
            for (int j = 0; j < hlt_njet; j++) {
                if (hlt_j_pt[j] > max_hlt_j_pt && trig->lower_eta <= hlt_j_eta[j] && hlt_j_eta[j] <= trig->upper_eta) {
                    max_hlt_j_pt = (double)hlt_j_pt[j];
                }
            }
            if (trig->referenceTrigger->m_trig_bool) { // only consider events where the reference trigger fired
                //p1 = (double)(m_trig_prescale[index]);
                //p2 = (double)(m_trig_prescale[referenceIndex]);
                //p_adj = 1./(1./p1 + 1./p2 - 1./(p1*p2));
                harr[index+numtrigs]->Fill(max_j_pt); // fill everytime that the reference trigger fired
     //           if (m_trig_bool[index]) harr[index]->Fill(max_j_pt); // fill if the trigger DID fire
                if (max_hlt_j_pt >= trig->min_pt[0]) harr[index]->Fill(max_j_pt); // fill if the trigger fired too
            }

            /*if (!m_trig_bool[referenceIndex]) continue;
            //p2 = (double)(m_trig_prescale[referenceIndex]);
            p2 = 1.;
            if (!m_trig_bool[index]) {
                harr[index+numtrigs]->Fill(jpt, p2);
            }
            else {
                p1 = (double)(m_trig_prescale[index]);
                p_adj = 1./(1./p1 + 1./p2 - 1./(p1*p2));
                p1 = p_adj;
                p2 = p_adj;
                harr[index]->Fill(jpt, p1);
                harr[index+numtrigs]->Fill(jpt, p2);
            }*/
        }       
    }
    if (debugStatements) cout << "Finished event loop for run " << thisRunNumber << endl;

    // Write histograms to a root file
    TFile* output = new TFile(Form("%srun_%i.root", effPath.c_str(), thisRunNumber), "RECREATE");
    for (Trigger* trig : triggerSubList) {
        index = trig->index;
        harr[index]->Write();
        harr[index+numtrigs]->Write();
    }
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    TVectorD run_vec(3);
    run_vec[0] = thisRunNumber;
    run_vec[1] = numetabins;
    run_vec[2] = numtrigs;
    run_vec.Write("run_vec");

    output->Close();

    return;
}
