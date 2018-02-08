#include "../triggerUtil.C"

void TriggerEfficiencyAnalysis(const int thisRunNumber, // Run number identifier.
                               double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;

    initialize(thisRunNumber, false);

    vector<Trigger*> triggerSubList(0);
    for (Trigger* trig : triggerVec) {
        if ((trig->lowerRunNumber <= thisRunNumber && thisRunNumber < trig->upperRunNumber) || trig->referenceTrigger == trig) triggerSubList.push_back(trig);
    }
    if (debugStatements) {
        cout << "Status: In TriggerEfficiencyAnalysis.C (37): Processing run " << thisRunNumber << " with triggers:" << endl;
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
                    if (debugStatements) cout << "Status: In TriggerEfficiencyAnalysis.C (58): Found " << fname.Data() << endl; 
                    if (fname.Contains(to_string(thisRunNumber))) {
                        tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
                        break;
                    }
                }
            }
        }
    }
    if (tree == NULL) {
        cout << "Error: In TriggerEfficiencyAnalysis.C (68): TTree not obtained for given run number. Quitting." << endl;
        return;
    }

    // Disable loading of unimportant branch values
    vector<string> interestingBranchNames = {"njet", "j_pt", "j_eta", "hlt_ion_j_pt", "hlt_ion_j_eta", "hlt_ion_njet", "hlt_j_pt", "hlt_j_eta", "hlt_njet", "nvert", "vert_type"};
    TObjArray* branches = (TObjArray*)(tree->GetListOfBranches());
    for (TObject* obj : *branches) {
        TString branchName = (TString)obj->GetName();
        if (debugStatements) cout << "Status: In TriggerEfficiencyAnalysis.C (77): Tree contains branch \"" << branchName.Data() << "\"" << endl;
        bool interestingBranch = false;
        for (string s : interestingBranchNames) {
            interestingBranch = interestingBranch || (branchName.Data() == s);
        }
        if (!interestingBranch) {
            for (Trigger* trig : triggerSubList) {
                interestingBranch = interestingBranch || (branchName == trig->name);
            }
        }
        if (!interestingBranch) {
            tree->SetBranchStatus(branchName, 0);
        }
    }

    // Create branching addresses:  
    float j_pt[60] = {};
    float j_eta[60] = {};
    int njet = 0;
    float hlt_j_pt[60] = {};
    float hlt_j_eta[60] = {};
    int hlt_njet = 0;
    int nvert = 0;
    int vert_type[60] = {};

    // Set branch addresses
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("njet", &njet);
    tree->SetBranchAddress("nvert", &nvert);
    tree->SetBranchAddress("vert_type", vert_type);
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
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->m_trig_bool));
    }

    int pbin, ebin, index, referenceIndex;
//    TH1F* histArr[2*numtrigs];
    TEfficiency* effArr[numtrigs];
    for (Trigger* trig : triggerSubList) {
        index = trig->index;
        effArr[index] = new TEfficiency(Form("%s_t_efficiency_run%i", trig->name.c_str(), thisRunNumber), "", numpbins, pbins);
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();
    double max_j_pt, max_hlt_j_pt;
    if (debugStatements) cout << "Status: In TriggerEfficiencyAnalysis.C (140): Looping over " << numentries << " events in run " << thisRunNumber << endl;
    for (long long entry = 0; entry < numentries; entry++) {
        tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses

        if ((nvert == 0) || (nvert > 0 && vert_type[0] != 1)) continue; // Basic event selection: require a primary vertex

        /**          BOOTSTRAPPING METHOD          **/
        /** Method works by calculating lowest pt  **/
        /** jet trigger fire rate relative to min  **/
        /** bias trigger, then 2nd lowest pt       **/
        /** relative to lowest pt multiplied by    **/
        /** lowest trigger efficiency, ad          **/
        /** infinitum.                             **/

        for (Trigger* trig : triggerSubList) {
            if (!(trig->referenceTrigger->m_trig_bool)) continue; // only consider events where the reference trigger fired
            index = trig->index;

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

            effArr[index]->Fill((max_hlt_j_pt >= trig->threshold_pt), max_j_pt);
        }       
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("%srun_%i.root", effPath.c_str(), thisRunNumber), "RECREATE");
    for (Trigger* trig : triggerSubList) {
        index = trig->index;
        effArr[index]->Write();
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
