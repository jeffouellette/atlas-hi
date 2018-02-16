#include "../triggerUtil.C"

const int numxbins = 50;
const int numqbins = 8;
const int numq2bins = 110;
const int numq2xbins = 100;
const int nummbins = 50;
const int numfcalbins = 60;

const double* xbins = logspace(1.6e-4, 1.6, numxbins);
const double* qbins = logspace(20, 1200, numqbins);
const double* q2bins = logspace(1, 500000, numq2bins);
const double* q2xbins = logspace(1.6e-4, 1.6, numq2xbins);
const double* mbins = logspace(20, 1200, nummbins);
const double* fcalbins = logspace(10, 500, numfcalbins);

void DijetAnalysis(int thisRunNumber, // Run number identifier.
                   double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;
    initialize(thisRunNumber, true);
    vector<TF1*>* triggerEfficiencyFunctions = getTriggerEfficiencyFunctions();

    const bool periodA = thisRunNumber < 313500;

    /**** Generate list of physics triggers ****/
    vector<Trigger*>* triggerSubvector = getTriggerSubvector(thisRunNumber);
    if (debugStatements) {
        cout << "Status: In DijetAnalysis.C (breakpoint A): Processing run " << thisRunNumber << " with triggers:" << endl;
        for (Trigger* trig : (*triggerSubvector)) {
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
                    if (debugStatements) cout << "Status: In DijetAnalysis.C (breakpoint B): Found " << fname.Data() << endl; 
                    if (fname.Contains(to_string(thisRunNumber))) {
                        tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
                        break;
                    }
                }
            }
        }
    }
    if (tree == NULL) {
        cout << "Error: In DijetAnalysis.C (breakpoint C): TTree not obtained for given run number. Quitting." << endl;
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
            if (debugStatements) cout << "Status: In DijetAnalysis.C (breakpoint D): Tree contains branch \"" << branchName.Data() << "\"" << endl;
            interestingBranch = false;
            for (string s : interestingBranchNames) {
                interestingBranch = interestingBranch || (branchName.Data() == s);
            }
            if (!interestingBranch) {
                for (Trigger* trig : (*triggerSubvector)) {
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


    // Create an array of 16 histograms, one for each rapidity region and one for x_p, x_a. 
    TH1D* xHistArr[2*numetabins];
    TH1D* xqHistArr[2*numqbins];
    TH1D* mHistArr[numetabins];

    for (int etabin = 0; etabin < numetabins; etabin++) {
        xHistArr[etabin] = new TH1D(Form("%ieta%i", thisRunNumber, etabin), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d^{2}N/L_{int}d#it{x}_{p}d#eta #left[nb#right]", etabins[etabin], etabins[etabin+1]), numxbins, xbins);
        xHistArr[etabin]->Sumw2();
    }
    for (int etabin = numetabins; etabin < 2*numetabins; etabin++) {
        xHistArr[etabin] = new TH1D(Form("%ieta%i", thisRunNumber, etabin), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d^{2}N/L_{int}d#it{x}_{a}d#eta #left[nb#right]", etabins[etabin%numetabins], etabins[(etabin%numetabins)+1]), numxbins, xbins);
        xHistArr[etabin]->Sumw2();
    }

    for (int qbin = 0; qbin < numqbins; qbin++) {
        xqHistArr[qbin] = new TH1D(Form("%iq%i", thisRunNumber, qbin), Form("%g < #it{Q} < %g;#it{x}_{p};d^{2}N/L_{int}d#it{x}_{p}d#it{Q} #left[nb GeV^{-1}#right]", qbins[qbin], qbins[qbin+1]), numxbins, xbins);
        xqHistArr[qbin]->Sumw2();
    }
    for (int qbin = numqbins; qbin < 2*numqbins; qbin++) {
        xqHistArr[qbin] = new TH1D(Form("%iq%i", thisRunNumber, qbin), Form("%g < #it{Q} < %g;#it{x}_{a};d^{2}N/L_{int}d#it{x}_{a}d#it{Q} #left[nb GeV^{-1}#right]", qbins[qbin%numqbins], qbins[(qbin%numqbins)+1]), numxbins, xbins);
        xqHistArr[qbin]->Sumw2();
    }

    for (int etabin = 0; etabin < numetabins; etabin++) {
        mHistArr[etabin] = new TH1D(Form("mjj_%ieta%i", thisRunNumber, etabin), Form("%1.1f < #eta < %1.1f;#it{m}_{JJ} #left[GeV#right];d^{2}N/L_{int}d#it{m}_{JJ}d#eta #left[nb GeV^{-1}#right]", etabins[etabin], etabins[etabin+1]), nummbins, mbins);
        mHistArr[etabin]->Sumw2();
    }

    TH2D* qxcorr = new TH2D(Form("xqcorr_run%i", thisRunNumber), ";#it{x}_{a};#it{#bar{Q}}^{2} #left[GeV^{2}#right];d^{2}N/L_{int}d{x}_{a}d#it{#bar{Q}}^{2} #left[nb GeV^{-2}#right]", numq2xbins, q2xbins, numq2bins, q2bins);
    TH2D* xaxpcorr = new TH2D(Form("xaxpcorr_run%i", thisRunNumber), ";#it{x}_{a};#it{x}_{p};d^{2}N/L_{int}d#it{x}_{p}d#it{x}_{a}", numxbins, xbins, numxbins, xbins);
    TH2D* fcalhist = new TH2D(Form("fcalhist_run%i", thisRunNumber), ";#it{x}_{p};FCAL energy #left[GeV#right];", numxbins, xbins, numfcalbins, fcalbins);

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
    tree->SetBranchAddress("j_phi", j_phi);
    tree->SetBranchAddress("j_e", j_e);
    tree->SetBranchAddress("njet", &njet);
    tree->SetBranchAddress("nvert", &nvert);
    tree->SetBranchAddress("vert_type", vert_type);
    if (!periodA) tree->SetBranchAddress("fcalC_et", &fcal_et);
    else tree->SetBranchAddress("fcalA_et", &fcal_et);
    tree->SetBranchAddress("eventNumber", &eventNumber);

    // Set branch addresses
    for (Trigger* trig : (*triggerSubvector)) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->m_trig_bool));
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &(trig->m_trig_prescale));
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();

    double leadingjpt, leadingjeta, subleadingjpt, subleadingjeta, leadingjphi, subleadingjphi, leadingje, subleadingje, subsubleadingjpt, deltaphi; // jet parameters
    double xp, xa, eff, lumi, scale, q_avg, mjj, etajj;
    int leadingj, subleadingj, subsubleadingj, etabin, actetabin, pbin, index, qbin;
    bool takeEvent;
    Trigger* bestTrigger;
    TLorentzVector leadingj_tlv;
    TLorentzVector subleadingj_tlv;
    TLorentzVector dijet_tlv;
    int numGoodEvents = 0;
    for (long long entry = 0; entry < numentries; entry++) {
        tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses
        if ((nvert == 0) || (nvert > 0 && vert_type[0] != 1) || njet < 2) continue; // Basic event selection: require a primary vertex and there to be at least 2 jets

        /** Find the leading dijet pair **/
        leadingj = 0;
        for (int j = 1; j < njet; j++) {
            if (j_pt[j] > j_pt[leadingj]) leadingj = j;
        }
        subleadingj = 0;
        for (int j = 1; j < njet; j++) {
            if (j == leadingj) continue; // ensures subleadingj != leadingj
            // if the subleading jet candidate IS the leading jet candidate ~OR~ if the next jet has a higher pt AND is not the leading jet, THEN jet j is a better subleading jet candidate.
            if (subleadingj == leadingj || j_pt[j] > j_pt[subleadingj]) subleadingj = j;
        }
        subsubleadingj = 0;
        for (int j = 1; j < njet; j++) {
            if (j == leadingj || j == subleadingj) continue; // ensures subsubleadingj != leadingj AND subsubleadingj != subleadingj
            if (subsubleadingj == leadingj || subsubleadingj == subleadingj || j_pt[j] > j_pt[subsubleadingj]) subsubleadingj = j;
        }

        /** Stores parameters of the leading dijet pair **/
        leadingjpt = (double)j_pt[leadingj];
        subleadingjpt = (double)j_pt[subleadingj];
        if (njet != 2) subsubleadingjpt = (double)j_pt[subsubleadingj];
        else subsubleadingjpt = 0;

        leadingjphi = (double)j_phi[leadingj];
        subleadingjphi = (double)j_phi[subleadingj];
        leadingjeta = (double)j_eta[leadingj];
        subleadingjeta = (double)j_eta[subleadingj];

        leadingje = (double)j_e[leadingj];
        subleadingje = (double)j_e[subleadingj];

        // make sure phi variables are in the range 0 to 2pi
        while (leadingjphi < 0) leadingjphi += 2*pi;
        while (subleadingjphi < 0) subleadingjphi += 2*pi;

        deltaphi = TMath::Abs(leadingjphi - subleadingjphi);
        if (deltaphi > pi) deltaphi = 2*pi - deltaphi;
        /** End find leading dijets **/

        /** Event selection **/
        if ((lowerPhiCut < leadingjphi && leadingjphi < upperPhiCut) || (lowerPhiCut < subleadingjphi && subleadingjphi < upperPhiCut)) continue; // select outside disabled HEC region
        if (deltaphi < 7.*pi/8.) continue; //require deltaPhi gap of 7pi/8
        if (leadingjpt < dijetMinimumPt || subleadingjpt < dijetMinimumPt) continue; // minimum pt cut
        if (subsubleadingjpt > thirdJetMaximumPt) continue;
        //if (subleadingjpt/leadingjpt < dijet_pt_ratio_cutoff) continue; // select on pt balance
        numGoodEvents++;
        /** End event selection **/


        if (leadingjpt > 1200) cout << Form("High pt (%.0f GeV) jet detected in run %i, event %i!", leadingjpt, thisRunNumber, eventNumber) << endl;

        etabin = getEtabin(leadingjeta);
        pbin = getPbin(leadingjpt);
        if (pbin == -1 || pbin > numpbins || etabin == -1 || etabin > numetabins) continue; // make sure we are in a valid kinematic bin

        bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins]; // kinematicTriggerVec is created per run, so we do not need to flip etas
        if (bestTrigger == NULL) continue; // make sure its not a null trigger

        if (periodA) actetabin = numetabins - etabin - 1;
        else actetabin = etabin;

        index = bestTrigger->index;
        lumi = kinematicLumiVec[pbin + actetabin*numpbins]; // kinematicLumiVec is created for all runs, so we DO need to flip etas!!!!
        eff = (*triggerEfficiencyFunctions)[index]->Eval(leadingjpt); // same with trigger efficiency, but there is no eta dependence
        takeEvent = bestTrigger->m_trig_bool && bestTrigger->m_trig_prescale > 0 && bestTrigger->min_pt <= leadingjpt && lumi != 0 && eff != 0.;
        if (!takeEvent) continue; // make sure the trigger fired, and that we're not going to be dividing by 0

        scale = 1e3/(eff*lumi); // set the scale to convert counts -> "efficiency corrected cross-section"-esque measurement
        //scale = 1;

        if (periodA) {
            leadingjeta *= -1;
            subleadingjeta *= -1;
        }

        xp = get_xp(leadingjpt, subleadingjpt, leadingjeta, subleadingjeta, false);
        xa = get_xa(leadingjpt, subleadingjpt, leadingjeta, subleadingjeta, false);

        // Fill hardness (Q^2) plots
        q_avg = TMath::Sqrt(0.5*(get_q2(xp, leadingje, leadingjpt) + get_q2(xp, subleadingje, subleadingjpt)));
        {
            qbin = 0;
            while (qbin < numqbins+1 && qbins[qbin++] < q_avg);
            qbin -= 2;
            if (!(qbin == -1 || qbin > numqbins)) {
                xqHistArr[qbin]->Fill(xp, scale);
                xqHistArr[qbin+numqbins]->Fill(xa, scale);
            }
        }
 
        // Fill q2 xa correlation plot
        qxcorr->Fill(xa, q_avg*q_avg, scale);

        // Fill xa xp correlation plot
        xaxpcorr->Fill(xa, xp, scale);

        // Fill Pb-going FCAL distribution plot
        fcalhist->Fill(xp, fcal_et, scale);

        leadingj_tlv.SetPtEtaPhiE(leadingjpt, leadingjeta, leadingjphi, leadingje);
        subleadingj_tlv.SetPtEtaPhiE(subleadingjpt, subleadingjeta, subleadingjphi, subleadingje);
        dijet_tlv = leadingj_tlv + subleadingj_tlv;

        mjj = dijet_tlv.Mag();
        etajj = dijet_tlv.Eta();

        etabin = getEtabin(etajj);
        if (etabin == -1 || etabin > numetabins) continue;
        
        mHistArr[etabin]->Fill(mjj, scale);
        // Fill xa xp distribution plots
        xHistArr[etabin]->Fill(xp, scale); // by filling with the period B etabin, we avoid having to flip the addition later on
        xHistArr[etabin+numetabins]->Fill(xa, scale);

    }
    // Save to root file
    string output_name = xPath;
    if (runPeriodA && !runPeriodB) output_name = output_name + "periodA/";
    else if (!runPeriodA && runPeriodB) output_name = output_name + "periodB/";
    else output_name = output_name + "periodAB/";
    output_name = Form("%srun_%i.root", output_name.c_str(), thisRunNumber);

    TFile* output = new TFile(output_name.c_str(), "RECREATE");
    for (int etabin = 0; etabin < 2*numetabins; etabin++) {
        scale = 1. / (etabins[(etabin%numetabins)+1] - etabins[etabin%numetabins]);
        xHistArr[etabin]->Scale(scale, "width");
        xHistArr[etabin]->Write();
    }
    for (int qbin = 0; qbin < 2*numqbins; qbin++) {
        scale = 1. / (qbins[(qbin%numqbins)+1] - qbins[qbin%numqbins]);
        xqHistArr[qbin]->Scale(scale, "width");
        xqHistArr[qbin]->Write();
    }

    for (int etabin = 0; etabin < numetabins; etabin++) {
        scale = 1. / (etabins[etabin+1] - etabins[etabin]);
        mHistArr[etabin]->Scale(scale, "width");
        mHistArr[etabin]->Write();
    }
    qxcorr->Scale(1., "width");
    xaxpcorr->Scale(1., "width");
    fcalhist->Scale(1., "width");

    qxcorr->Write();
    xaxpcorr->Write();
    fcalhist->Write();

    TVectorD run_vec(2);
    run_vec[0] = luminosity;
    run_vec[1] = numGoodEvents;
    run_vec.Write("run_vec");
    output->Close();
}
