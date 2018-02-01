#include "../triggerUtil.C"

void bootstrap(Trigger* trig, TH1F** histArr) {
    if (trig->referenceTrigger->index == trig->index) {
        for (int pbin = 0; pbin < numpbins; pbin++) {
            histArr[trig->referenceTrigger->index]->SetBinContent(pbin+1, 1);
            histArr[trig->referenceTrigger->index]->SetBinError(pbin+1, 0);
        }
        trig->isBootstrapped = true;
    }
    if (trig->isBootstrapped) return;

    int index = trig->index;
    int referenceIndex = trig->referenceTrigger->index;
    if (!(trig->referenceTrigger->isBootstrapped)) {
        if (debugStatements) cout << "Status: In trigger_efficiencies.C (16): Bootstrapping trigger " << trig->referenceTrigger->name << " for " << trig->name << endl;
        bootstrap(trig->referenceTrigger, histArr);
    }
    histArr[index]->Multiply(histArr[referenceIndex]);
    trig->isBootstrapped = true;
    return;
}


void TriggerEfficiencyAnalysis(const int thisRunNumber, // Run number identifier.
                               double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;

    initialize(thisRunNumber, false, false);

    vector<Trigger*> triggerSubList(0);
    for (Trigger* trig : triggerVec) {
        if (!trig->disabled) triggerSubList.push_back(trig);
    }
    if (debugStatements) {
        cout << "Status: In trigger_efficiencies.C (37): Processing run " << thisRunNumber << " with triggers:" << endl;
        for (Trigger* trig : triggerVec) {
            cout << "\t" << trig->name << endl;
        }
    }

    luminosity *= 1e-3; // convert from nb^(-1) to pb^(-1)
    const int numhists = numtrigs;

    TTree* tree = NULL;
    {
        TSystemDirectory dir(dataPath.c_str(), dataPath.c_str());
        TList* sysfiles = dir.GetListOfFiles();
        if (sysfiles) {
            TSystemFile *sysfile;
            TString fname;
            TIter next(sysfiles);

            while ((sysfile=(TSystemFile*)next())) {
                fname = sysfile->GetName();
                if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
                    if (debugStatements) cout << "Status: In trigger_efficiencies.C (58): Found " << fname.Data() << endl; 
                    if (fname.Contains(to_string(thisRunNumber))) {
                        tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
                        break;
                    }
                }
            }
        }
    }
    if (tree == NULL) {
        cout << "Error: In trigger_efficiencies.C (68): TTree not obtained for given run number. Quitting." << endl;
        return;
    }

    // Disable loading of unimportant branch values
    vector<string> interestingBranchNames = {"njet", "j_pt", "j_eta", "hlt_ion_j_pt", "hlt_ion_j_eta", "hlt_ion_njet", "hlt_j_pt", "hlt_j_eta", "hlt_njet"};
    TObjArray* branches = (TObjArray*)(tree->GetListOfBranches());
    for (TObject* obj : *branches) {
        TString branchName = (TString)obj->GetName();
        if (debugStatements) cout << "Status: In trigger_efficiencies.C (77): Tree contains branch \"" << branchName.Data() << "\"" << endl;
        bool interestingBranch;
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

    // Create branching addresses:  
    // Create arrays to store trigger values for each event
    bool m_trig_bool[numtrigs];   // stores whether trigger was triggered
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
    if (thisRunNumber <= 313603) {
        tree->SetBranchAddress("hlt_ion_j_pt", hlt_j_pt);
        tree->SetBranchAddress("hlt_ion_j_eta", hlt_j_eta);
        tree->SetBranchAddress("hlt_ion_njet", &hlt_njet);
    } else {
        tree->SetBranchAddress("hlt_j_pt", hlt_j_pt);
        tree->SetBranchAddress("hlt_j_eta", hlt_j_eta);
        tree->SetBranchAddress("hlt_njet", &hlt_njet);
    }
    for (Trigger* trig : triggerSubList) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &m_trig_bool[trig->index]);
    }

    int pbin, ebin, index, referenceIndex;
    TH1F* histArr[2*numtrigs];
    for (Trigger* trig : triggerSubList) {
        index = trig->index;
        // standard trigger firings
        TString histname = Form("%s_efficiency_run%i", trig->name.c_str(), thisRunNumber);
        histArr[index] = new TH1F(histname, ";#it{p}_{T}^{jet} #left[GeV#right];Efficiency #epsilon", numpbins, pbins);
        histArr[index]->Sumw2(); // instruct each histogram to propagate errors
        // reference trigger firings
        histname = Form("%s_reference_run%i", trig->name.c_str(), thisRunNumber);
        histArr[index+numtrigs] = new TH1F(histname, ";#it{p}_{T}^{jet} #left[GeV#right];Efficiency #epsilon", numpbins, pbins);
        histArr[index+numtrigs]->Sumw2();
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();
    double max_j_pt, max_hlt_j_pt, p1, p2, p_adj;
    if (debugStatements) cout << "Status: In trigger_efficiencies.C (140): Looping over " << numentries << " events in run " << thisRunNumber << endl;
    for (int entry = 0; entry < numentries; entry++) {
        tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses

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
            if (m_trig_bool[referenceIndex]) { // only consider events where the reference trigger fired
                histArr[index+numtrigs]->Fill(max_j_pt); // fill everytime that the reference trigger fired
                if (max_hlt_j_pt >= trig->min_pt) histArr[index]->Fill(max_j_pt); // fill if the trigger fired too
            }
        }       
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("%srun_%i.root", effPath.c_str(), thisRunNumber), "RECREATE");
    for (Trigger* trig : triggerSubList) {
        index = trig->index;
        histArr[index]->Write();
        histArr[index+numtrigs]->Write();
    }
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write(Form("lum_vec_%i", thisRunNumber));
    TVectorD run_vec(3);
    run_vec[0] = thisRunNumber;
    run_vec[1] = numetabins;
    run_vec[2] = numtrigs;
    run_vec.Write(Form("run_vec_%i", thisRunNumber));

    output->Close();

    if (debugStatements) cout << "Status: In triggerEfficiencies.C (192):Finished event loop for run " << thisRunNumber << endl;
    return;
}
