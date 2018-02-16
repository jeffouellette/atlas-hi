#include "../triggerUtil.C"

void IdealPtAnalysis(const int thisRunNumber, // Run number identifier.
                       double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;

    initialize(thisRunNumber, true);
    vector<TF1*>* triggerEfficiencyFunctions = getTriggerEfficiencyFunctions();
    
    /**** Generate list of physics triggers ****/
    vector<Trigger*>* triggerSubvector = getTriggerSubvector(thisRunNumber);
    if (debugStatements) {
        cout << "Status: In IdealPtAnalysis.C (breakpoint A): Processing run " << thisRunNumber << " with triggers:" << endl;
        for (Trigger* trig : (*triggerSubvector)) {
            cout << "\t" << trig->name << endl;
        }
    }
    /**** End generate list of physics triggers ****/


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
                    if (debugStatements) cout << "Status: In IdealPtAnalysis.C (breakpoint B): Found " << fname.Data() << endl; 
                    if (fname.Contains(to_string(thisRunNumber))) {
                        tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
                        break;
                    }
                }
            }
        }
    }
    if (tree == NULL) {
        cout << "Error: In IdealPtAnalysis.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
        return;
    }
    /**** End find TTree ****/


    /**** Disable loading of unimportant branch values - speeds up entry retrieval ****/
    {
        vector<string> interestingBranchNames = {"njet", "j_pt", "j_eta", "j_phi"};
        TObjArray* branches = (TObjArray*)(tree->GetListOfBranches());
        bool interestingBranch;
        for (TObject* obj : *branches) {
            TString branchName = (TString)obj->GetName();
            if (debugStatements) cout << "Status: In IdealPtAnalysis.C (breakpoint D): Tree contains branch \"" << branchName.Data() << "\"" << endl;
            interestingBranch = false;
            for (string s : interestingBranchNames) {
                interestingBranch = interestingBranch || (branchName.Data() == s);
            }
            if (!interestingBranch) {
                for (Trigger* trig : (*triggerSubvector)) {
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
    float j_pt[60] = {};
    float j_eta[60] = {};
    float j_phi[60] = {};
    int njet = 0;

    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("j_phi", j_phi);
    tree->SetBranchAddress("njet", &njet);
    for (Trigger* trig : (*triggerSubvector)) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->m_trig_bool));
    }
    /**** End set branch addresses ****/


    /**** Histogram initialization ****/
    const int numhists = numtrigs * numetabins;
    TH1F* histArr[numhists];
    TH2F* etaPhiHist;
    int pbin, etabin, index;
    for (etabin = 0; etabin < numetabins; etabin++) {
        TString histName = Form("trig_pt_counts_run%i_etabin%i", thisRunNumber, etabin);
        histArr[etabin] = new TH1F(histName, ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/Ad#it{p}_{T}d#eta #left[nb GeV^{-1}#right]", numpbins, pbins);
        histArr[etabin]->Sumw2(); // instruct each histogram to propagate errors
    }
    {
        TString histName = Form("etaPhiHist_run%i", thisRunNumber);
        etaPhiHist = new TH2F(histName, ";#eta;#phi;", 98, -4.9, 4.9, 100, 0, 2*TMath::Pi());
    }

    /**** Iterate over each event ****/
    const int numentries = tree->GetEntries();

    double jpt, jeta, jphi, eff;
    Trigger* bestTrigger = NULL;
    for (long long entry = 0; entry < numentries; entry++) {
        tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses

        for (int j = 0; j < njet; j++) {
            jpt = (double)j_pt[j];
            jeta = (double)j_eta[j];
            jphi = (double)j_phi[j];
            while (jphi < 0) jphi += 2.*pi;

            etabin = getEtabin(jeta);
            pbin = getPbin(jpt);

            if (pbin < 0 || etabin < 0 || pbin > numpbins || etabin > numetabins) continue; // this checks that the jets fall within the pt, eta bins

            bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins];
            //eff = kinematicEfficiencyVec[pbin + etabin*numpbins]; 
            if (bestTrigger == NULL) continue; // make sure we're not trying to look at a null trigger
            eff = (*triggerEfficiencyFunctions)[bestTrigger->index]->Eval(jpt);
            if (eff == 0) continue; // avoid dividing by 0 
            if (bestTrigger->m_trig_bool) {
                if (jphi <= lowerPhiCut || jphi >= upperPhiCut) histArr[etabin]->Fill(jpt, ((2*pi)/(2*pi- (upperPhiCut-lowerPhiCut)))*(1./eff));
                //if (thisRunNumber < 313500) etaPhiHist->Fill(-jeta, jphi);
                //else etaPhiHist->Fill(jeta, jphi);
                etaPhiHist->Fill(jeta, jphi);
            }
        }
    }
    /**** End event iteration ****/


    /**** Write output histograms to a root file ****/
    TFile* output = new TFile(Form("%srun_%i.root", ptPath.c_str(), thisRunNumber), "RECREATE");
    for (etabin = 0; etabin < numetabins; etabin++) {
        histArr[etabin]->Scale(1/A); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        histArr[etabin]->Write();
    }
    etaPhiHist->Write();

    TVectorD run_vec(5);
    run_vec[0] = thisRunNumber;
    run_vec[1] = numetabins;
    run_vec[2] = numtrigs;
    run_vec[3] = numpbins;
    run_vec[4] = luminosity;
    run_vec.Write(Form("run_vec_%i", thisRunNumber));

    output->Close();
    /**** End write output ****/

    if (debugStatements) cout << "Status: In IdealPtAnalysis.C (breakpoint E): Finished calculating pt spectrum for run " << thisRunNumber << endl;
    return;
}
