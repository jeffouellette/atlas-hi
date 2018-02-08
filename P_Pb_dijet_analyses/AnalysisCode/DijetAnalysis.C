#include "../triggerUtil.C"

void DijetAnalysis(int thisRunNumber, // Run number identifier.
                double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;
    initialize(thisRunNumber, true);
    vector<TF1*>* triggerEfficiencyFunctions = getTriggerEfficiencyFunctions();

    const int numxbins = 100;
    const int numqbins = 8;
    const int numhists = 2*numetabins;
    const bool periodA = thisRunNumber < 313500;
    const double histArrScales[8] = {1, 1, 1, 1, 1, 1, 1, 1};   // rescaling factors so the histograms don't overlap- all 1 right now during debugging
    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)


    /**** Generate list of physics triggers ****/
    vector<Trigger*> triggerSubList(0);
    for (Trigger* trig : triggerVec) {
        if (trig->lowerRunNumber <= thisRunNumber && thisRunNumber < trig->upperRunNumber && trig->name != minbiasTriggerName) triggerSubList.push_back(trig);
    }
    if (debugStatements) {
        cout << "Status: In IdealPtAnalysis.C (16): Processing run " << thisRunNumber << " with triggers:" << endl;
        for (Trigger* trig : triggerSubList) {
            cout << "\t" << trig->name << endl;
        }
    }
    /**** End generate list of physics triggers ****/


    /**** Find the relevant TTree for this run ****/
    //TTree* tree = (TTree*)(new TFile(Form("%srun_%i_raw.root", dataPath.c_str(), thisRunNumber)))->Get("tree");
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
                    if (debugStatements) cout << "Status: In IdealPtAnalysis.C (39): Found " << fname.Data() << endl; 
                    if (fname.Contains(to_string(thisRunNumber))) {
                        tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
                        break;
                    }
                }
            }
        }
    }
    if (tree == NULL) {
        cout << "Error: In IdealPtAnalysis.C (49): TTree not obtained for given run number. Quitting." << endl;
        return;
    }
    /**** End find TTree ****/

    /**** Disable loading of unimportant branch values - speeds up entry retrieval ****/
    {
        vector<string> interestingBranchNames = {"njet", "j_pt", "j_eta", "j_phi", "j_e", "eventNumber", "fcalA_et", "fcalC_et", "nvert", "vert_type"};
        TObjArray* branches = (TObjArray*)(tree->GetListOfBranches());
        bool interestingBranch;
        for (TObject* obj : *branches) {
            TString branchName = (TString)obj->GetName();
            if (debugStatements) cout << "Status: In IdealPtAnalysis.C (62): Tree contains branch \"" << branchName.Data() << "\"" << endl;
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


    const double* xbins = logspace(8e-5, 1.6, numxbins);
    const double* qbins = logspace(20, 1000, numqbins);
    // Create an array of 16 histograms, one for each rapidity region and one for x_p, x_a. 
    TH1D* histArr[numhists];
    for (int i = 0; i < numhists/2; i++) {
        histArr[i] = new TH1D(Form("%ieta%i", thisRunNumber, i), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d^{2}#sigma/d#it{x}_{p} dy #left[pb#right]", etabins[i], etabins[i+1]), numxbins, xbins);
        histArr[i]->Sumw2();
    }
    for (int i = numhists/2; i < numhists; i++) {
        histArr[i] = new TH1D(Form("%ieta%i", thisRunNumber, i), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d^{2}#sigma/d#it{x}_{a} dy #left[pb#right]", etabins[i%(numhists/2)], etabins[(i%(numhists/2))+1]), numxbins, xbins);
        histArr[i]->Sumw2();
    }
    TH2D* xaxpcorr = new TH2D(Form("xaxpcorr_run%i", thisRunNumber), ";#it{x}_{a};#it{x}_{p};d^{2}#sigma/d#it{x}_{p}d#it{x}_{a}", numxbins, xbins, numxbins, xbins);

    const int numfcalbins = 60;
    const double* fcalbins = logspace(10, 500, numfcalbins);
    TH2D* fcalhist = new TH2D(Form("fcalhist_run%i", thisRunNumber), ";#it{x}_{p};FCAL energy deposited;", numxbins, xbins, numfcalbins, fcalbins);
//    cout << numtrigs << endl;

    // Create arrays to store trigger values for each event
//    bool m_trig_bool[numtrigs];
//    float m_trig_prescale[numtrigs];

    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    float j_phi[60] = {};
    float j_e[60] = {};
    int vert_type[60] = {};
    int nvert = 0;
    int eventNumber = 0;
    int njet = 0;
    float fcal_et = 0;
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("j_e", j_e);
    tree->SetBranchAddress("njet", &njet);
    tree->SetBranchAddress("j_phi", j_phi);
    tree->SetBranchAddress("nvert", &nvert);
    tree->SetBranchAddress("vert_type", vert_type);
    if (!periodA) tree->SetBranchAddress("fcalC_et", &fcal_et);
    else tree->SetBranchAddress("fcalA_et", &fcal_et);
    tree->SetBranchAddress("eventNumber", &eventNumber);

    // Set branch addresses
    for (Trigger* trig : triggerSubList) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->m_trig_bool));
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &(trig->m_trig_prescale));
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();

    double leadingjpt, leadingjeta, subleadingjpt, subleadingjeta, leadingjphi, subleadingjphi, leadingje, subleadingje, xp, xa, eff, lumi, scale, hardness_q;
    int leadingj, subleadingj, etabin, pbin, index, qbin;
    bool takeEvent;
    Trigger* bestTrigger;
    for (long long entry = 0; entry < numentries; entry++) {
        tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses
        if ((nvert == 0) || (nvert > 0 && vert_type[0] != 1)) continue; // Basic event selection: require a primary vertex

        leadingj = 0;
        subleadingj = 1;
        for (int j = 1; j < njet; j++) {
            if (j_pt[j] > (double)j_pt[leadingj]) {
                subleadingj = leadingj;
                leadingj = j;
            }
            else if (j_pt[j] > (double)j_pt[subleadingj]) {
                subleadingj = j;
            }
        }
        
        leadingjeta = (double)j_eta[leadingj];
        subleadingjeta = (double)j_eta[subleadingj];
        leadingjpt = (double)j_pt[leadingj];
        subleadingjpt = (double)j_pt[subleadingj];
        leadingjphi = (double)j_phi[leadingj];
        subleadingjphi = (double)j_phi[subleadingj];
        leadingje = (double)j_e[leadingj];
        subleadingje = (double)j_e[subleadingj];

        if (leadingjpt > 2000) cout << Form("High pt (%.0f GeV) jet detected in run %i, event %i!", leadingjpt, thisRunNumber, eventNumber) << endl;

        // Consider events where the dijet ratio meets the cutoff
        takeEvent = subleadingjpt/leadingjpt >= dijet_pt_ratio_cutoff;
        if (!takeEvent) continue;

        etabin = getEtabin(leadingjeta);
        pbin = getPbin(leadingjpt);
        if (pbin == -1 || pbin > numpbins || etabin == -1 || etabin > numetabins) continue; // make sure we are in a valid kinematic bin

        bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins];
        if (bestTrigger == NULL) continue; // make sure its not a null trigger

        index = bestTrigger->index;
        lumi = kinematicLumiVec[pbin + etabin*numpbins];
        takeEvent = bestTrigger->m_trig_bool && bestTrigger->m_trig_prescale > 0;
        if (!takeEvent) continue; // make sure the trigger fired and was not prescaled

        xp = get_xp(leadingjpt, subleadingjpt, leadingjeta, subleadingjeta, periodA);
        xa = get_xa(leadingjpt, subleadingjpt, leadingjeta, subleadingjeta, periodA);
        // Equivalent to:
        // xp = get_xp(leadingjpt, subleadingjpt, -leadingjeta, -subleadingjeta, false);
        // xa = get_xa(leadingjpt, subleadingjpt, -leadingjeta, -subleadingjeta, false);

        //hardness_q = 0.5*(get_q(xp, leadingje, leadingjpt) + get_q(xp, subleadingje, subleadingjpt));
        /*qbin = 0;
        while (qbin < numqbins+1 && qbins[qbin++] < hardness_q);
        qbin -= 2;
        if (qbin == -1 || qbin > numqbins) continue;*/

        lumi = kinematicLumiVec[pbin + etabin*numpbins];
        eff = (*triggerEfficiencyFunctions)[index]->Eval(leadingjpt);

        if (lumi == 0 || eff == 0) continue; // avoid dividing by 0 

        scale = 1e3/(eff*lumi);
        histArr[etabin]->Fill(xp, scale);
        histArr[etabin+numetabins]->Fill(xa, scale);
        if ((!periodA && leadingjeta <= 3.2) || (periodA && leadingjeta >= -3.2)) {
            xaxpcorr->Fill(xa, xp, scale);
            fcalhist->Fill(xp, fcal_et, scale);
        }
    }
    // Save to root file
    string output_name = xPath;
    if (runPeriodA && !runPeriodB) output_name = output_name + "periodA/";
    else if (!runPeriodA && runPeriodB) output_name = output_name + "periodB/";
    else output_name = output_name + "periodAB/";
    output_name = Form("%srun_%i.root", output_name.c_str(), thisRunNumber);

    TFile* output = new TFile(output_name.c_str(), "RECREATE");
    for (int i = 0; i < numhists; i++) {
        histArr[i]->Scale(histArrScales[i%(numetabins)]/(etabins[(i%numetabins)+1] - etabins[(i%numetabins)]), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        histArr[i]->Write();
    }
    xaxpcorr->Write();
    fcalhist->Write();
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    output->Close();
}
