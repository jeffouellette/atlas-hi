#include "../triggerUtil.C"

const double etaCoM = -0.465; // boost into CoM frame in period B kinematics
const int numppPtbins = 36;
const int numppYbins = 12;
const int numpPbYbins = 18;
double ppPtbins[(numppPtbins+1)*numppYbins] = {}; //= new double[(numppPtbins+1)*numppYbins];
const double ppYbins[numppYbins+1] = {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3};
double pPbYbins[numpPbYbins+1] = {-3.465, -3.2, -2.965, -2.465, -2, -1.965, -1.465, -1, -0.965, -0.465, 0, 0.035, 0.535, 1, 1.035, 1.535, 2, 2.035, 2.535};

int getppPbin (double pt) {
    int bin = 0;
    while (bin < numppPtbins+1 && ppPtbins[bin] < pt) bin++;
    return bin - 1;
}


int getppYbin (double eta) {
    int bin = 0;
    while (bin < numppYbins+1 && ppYbins[bin] < eta) bin++;
    return bin - 1;
}


int getpPbYbin (double eta) {
    int bin = 0;
    while (bin < numpPbYbins+1 && pPbYbins[bin] < eta) bin++;
    return bin - 1;
}

/**
 * Loads binning information from the 8TeV pp reference data.
 */
TH1D** setupPPConfiguration() {

    double* ppJetSpectrum = new double[numppPtbins*numppYbins];
    double* ppJetSpectrumStatError = new double[numppPtbins*numppYbins];

    for (int ppPtbin = 0; ppPtbin < numppPtbins; ppPtbin++) {
        for (int ppYbin = 0; ppYbin < numppYbins; ppYbin++) {
            ppPtbins[ppPtbin + ppYbin*(numppPtbins+1)] = 0;
            ppJetSpectrum[ppPtbin + ppYbin*numppPtbins] = 0;
            ppJetSpectrumStatError[ppPtbin + ppYbin*numppPtbins] = 0;
        }
    }

    string path = workPath + "HepData (#28325609)/";
    string dummyline;
    double pt_center, pt_left, pt_right, dsigma, stat_err;
    int ppPtbin, ppYbin, ppYbinComplement;
    vector<string> files = {"atlas_2012_jet_antiktr04_incljetpt_eta6_highmu_data.txt", "atlas_2012_jet_antiktr04_incljetpt_eta5_highmu_data.txt", "atlas_2012_jet_antiktr04_incljetpt_eta4_highmu_data.txt", "atlas_2012_jet_antiktr04_incljetpt_eta3_highmu_data.txt", "atlas_2012_jet_antiktr04_incljetpt_eta2_highmu_data.txt", "atlas_2012_jet_antiktr04_incljetpt_eta1_highmu_data.txt"};

    ppYbin = 0;
    ppYbinComplement = 11;
    for (string file : files) {
        ifstream dataStream;
        dataStream.open(path + file);
        if (!dataStream.is_open()) {
            cout << "Error: In IdealRpPbAnalysis.C (breakpoint A): Cannot find input file!" << endl;
            throw runtime_error("File not found");
        }
        for (int i = 0; i < 641; i++) getline(dataStream, dummyline); // first 640 lines are systematic errors so skip them, 641st line is a table layout

        // Creates a set of padding bins
        ppPtbin = 0;
        ppPtbins[ppPtbin + ppYbin*(numppPtbins+1)] = 60;
        ppPtbins[ppPtbin + ppYbinComplement*(numppPtbins+1)] = 60;
        ppPtbins[ppPtbin+1 + ppYbin*(numppPtbins+1)] = 60;
        ppPtbins[ppPtbin+1 + ppYbinComplement*(numppPtbins+1)] = 60;
        ppPtbin++;

        // Loads the pp spectrum pt bins
        while (dataStream >> pt_center) {
            if (ppPtbin > numppPtbins && debugStatements) cout << "Warning: In IdealRpPbAnalysis.C (breakpoint B): pt bin exceeding range of allocated values. Binning may be overwritten." << endl;
            dataStream >> pt_left;
            dataStream >> pt_right;
            dataStream >> dsigma;
            dataStream >> stat_err;
            ppPtbins[ppPtbin + ppYbin*(numppPtbins+1)] = pt_left;
            ppPtbins[ppPtbin + ppYbinComplement*(numppPtbins+1)] = pt_left;
            ppPtbins[ppPtbin+1 + ppYbin*(numppPtbins+1)] = pt_right;
            ppPtbins[ppPtbin+1 + ppYbinComplement*(numppPtbins+1)] = pt_right;
            ppPtbin++;
        }

        // Creates a set of padding bins
        while(ppPtbin < numppPtbins) {
            ppPtbins[ppPtbin + ppYbin*(numppPtbins+1)] = pt_right;
            ppPtbins[ppPtbin + ppYbinComplement*(numppPtbins+1)] = pt_right;
            ppPtbins[ppPtbin+1 + ppYbin*(numppPtbins+1)] = 2600;
            ppPtbins[ppPtbin+1 + ppYbinComplement*(numppPtbins+1)] = 2600;
            ppPtbin++;
        }

        ppYbin++;
        ppYbinComplement--;

        dataStream.close();
    }

    TH1D** ppHistArr = new TH1D*[numppYbins];
    for (ppYbin = 0; ppYbin < numppYbins; ppYbin++) {
        TString histName = Form("pp_spectrum_ppYbin%i", ppYbin);
        ppHistArr[ppYbin] = new TH1D (histName, ";#it{p}_{T}^{jet} #left[GeV#right];d^{2}#sigma/d#it{p}_{T}d#it{y} #left[pb GeV^{-1}#right]", numppPtbins, &(ppPtbins[ppYbin*(numppPtbins+1)]));
    }

    ppYbin = 0;
    ppYbinComplement = 11;
    for (string file : files) {
        ifstream dataStream;
        dataStream.open(path + file);
        if (!dataStream.is_open()) {
            cout << "Error: In IdealRpPbAnalysis.C (breakpoint C): Cannot find input file!" << endl;
            throw runtime_error("File not found");
        }
        for (int i = 0; i < 641; i++) getline(dataStream, dummyline); // first 640 lines are systematic errors so skip them, 641st line is a table layout
        ppPtbin = 1;
        while (dataStream >> pt_center) {
            if (ppPtbin > numppPtbins && debugStatements) cout << "Warning: In IdealRpPbAnalysis.C (breakpoint D): pt bin exceeding range of allocated values. Binning may be overwritten." << endl;
            dataStream >> pt_left;
            dataStream >> pt_right;
            dataStream >> dsigma;
            dataStream >> stat_err;

            // Rescale the cross section to be accurate in this bin
            double dy = 2*(ppYbins[ppYbin+1]-ppYbins[ppYbin]);
            dsigma = dsigma * dy * 0.5; // cross section in negative+positive rapidity (y) bin is half of the integrated cross section over that bin
            dy = ppYbins[ppYbin+1] - ppYbins[ppYbin]; // always = 0.5
            dsigma = dsigma / dy; // now divide out by the new negative OR positive rapidity bin width
            // given the current binning structure, this won't actually do anything.
            
            ppHistArr[ppYbin]->SetBinContent(ppPtbin+1, dsigma);
            ppHistArr[ppYbinComplement]->SetBinContent(ppPtbin+1, dsigma);
            ppHistArr[ppYbin]->SetBinError(ppPtbin+1, stat_err*dsigma/100.);
            ppHistArr[ppYbinComplement]->SetBinError(ppPtbin+1, stat_err*dsigma/100.);
            ppPtbin++;
        }
        ppYbin++;
        ppYbinComplement--;

        dataStream.close();
    }

    return ppHistArr; 
}


void IdealRpPbAnalysis(const int thisRunNumber, // Run number identifier.
                       double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;

    TH1D** ppHistArr = setupPPConfiguration();
    for (int ppYbin = 0; ppYbin < numppYbins; ppYbin++) delete ppHistArr[ppYbin]; // don't really need these right now so just delete them

    initialize(thisRunNumber, true);
    vector<TF1*>* triggerEfficiencyFunctions = getTriggerEfficiencyFunctions();
    
    /**** Generate list of physics triggers ****/
    vector<Trigger*>* triggerSubvector = getTriggerSubvector(thisRunNumber);
    if (debugStatements) {
        cout << "Status: In IdealRpPbAnalysis.C (breakpoint E): Processing run " << thisRunNumber << " with triggers:" << endl;
        for (Trigger* trig : (*triggerSubvector)) {
            cout << "\t" << trig->name << endl;
        }
    }
    /**** End generate list of physics triggers ****/

    const int numhists = numtrigs * numppYbins;

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
                    if (debugStatements) cout << "Status: In IdealRpPbAnalysis.C (breakpoint F): Found " << fname.Data() << endl; 
                    if (fname.Contains(to_string(thisRunNumber))) {
                        tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
                        break;
                    }
                }
            }
        }
    }
    if (tree == NULL) {
        cout << "Error: In IdealRpPbAnalysis.C (breakpoint G): TTree not obtained for given run number. Quitting." << endl;
        return;
    }
    /**** End find TTree ****/


    /**** Disable loading of unimportant branch values - speeds up entry retrieval ****/
    {
        vector<string> interestingBranchNames = {"njet", "jet_pt", "jet_eta", "jet_phi", "jet_e", "vert_type", "nvert"};
        TObjArray* branches = (TObjArray*)(tree->GetListOfBranches());
        bool interestingBranch;
        for (TObject* obj : *branches) {
            TString branchName = (TString)obj->GetName();
            if (debugStatements) cout << "Status: In IdealRpPbAnalysis.C (breakpoint H): Tree contains branch \"" << branchName.Data() << "\"" << endl;
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


    /**** Load eta-phi histogram for rescaling ****/
    double etaPhiScaleFactors[numpPbYbins];
    {
        TFile* etaPhiFile = new TFile((rootPath + "etaPhiHist_highPtJetsOnly.root").c_str(), "READ");
        TH2D* etaPhiHist = (TH2D*)etaPhiFile->Get("etaPhiHist");
        const int nbins_x = etaPhiHist->GetNbinsX();
        const int nbins_y = etaPhiHist->GetNbinsY();
        for (int pPbYbin = 0; pPbYbin < numpPbYbins; pPbYbin++) {
            double numerator = 0;
            double denominator = 0;
            double x, y, content;
            for (int bin_x = 0; bin_x < nbins_x; bin_x++) {
                x = etaPhiHist->GetXaxis()->GetBinCenter(bin_x+1);
                //if (!(pPbYbins[pPbYbin] < x && x < pPbYbins[pPbYbin+1])) continue; // if we do not overlap with these bins then just continue

                for (int bin_y = 0; bin_y < nbins_y; bin_y++) {
                    y = etaPhiHist->GetYaxis()->GetBinCenter(bin_y+1);
                    content = etaPhiHist->GetBinContent(bin_x+1, bin_y+1);
                    if (!(lowerPhiCut < y && y < upperPhiCut && lowerEtaCut < x && x < upperEtaCut)) denominator += content; // if we're outside the HEC, add to the denominator
                    numerator += content; // always add to the numerator
                }
            }
            if (denominator == 0. && debugStatements) cout << "Warning: In IdealRpPbAnalysis.C (breakpoint I): No jets found outside HEC region!" << endl;
            else if (denominator == 0.) etaPhiScaleFactors[pPbYbin] = 0;
            else etaPhiScaleFactors[pPbYbin] = numerator/denominator;
        }
        etaPhiFile->Close();
    }
    /**** End load eta-phi histogram ****/


    /**** Set branching addresses ****/
    // Create branching addresses:  
    // Create arrays to store trigger values for each event
    //bool m_trig_bool[numtrigs];   // stores whether trigger was triggered
    //float m_trig_prescale[numtrigs];      // stores the prescaling factor for the trigger
    // Create arrays to store jet data for each event
    int njet = 0;
    vector<float>* jet_pt = NULL;
    vector<float>* jet_eta = NULL;
    vector<float>* jet_phi = NULL;
    vector<float>* jet_e = NULL;
    /*float jet_pt[60] = {};
    float jet_eta[60] = {};
    float jet_phi[60] = {};
    float jet_e[60] = {};*/
    int nvert = 0;
    vector<float>* vert_type = NULL;
    //int vert_type[60] = {};

    // Set branch addresses
    tree->SetBranchAddress("njet", &njet);
    tree->SetBranchAddress("jet_pt", &jet_pt);
    tree->SetBranchAddress("jet_eta", &jet_eta);
    tree->SetBranchAddress("jet_phi", &jet_phi);
    tree->SetBranchAddress("jet_e", &jet_e);
    tree->SetBranchAddress("nvert", &nvert);
    tree->SetBranchAddress("vert_type", &vert_type);
    for (Trigger* trig : (*triggerSubvector)) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->m_trig_bool));
    }
    /**** End set branch addresses ****/


    /**** Histogram initialization ****/
    TH1D* pPbCoMHistArr[numhists];
    int pPbYbin, index;
    for (pPbYbin = 0; pPbYbin < numpPbYbins; pPbYbin++) {
        TString histName = Form("pPb_spectrum_run%i_pPbYbin%i", thisRunNumber, pPbYbin);
        int ppYbin_equiv = getppYbin(0.5*(pPbYbins[pPbYbin]+pPbYbins[pPbYbin+1]) - etaCoM);
        pPbCoMHistArr[pPbYbin] = new TH1D(histName, ";#it{p}_{T}^{jet} #left[GeV#right];d#sigma/Ad#it{p}_{T} #left[nb GeV^{-1}#right]", numppPtbins, &(ppPtbins[ppYbin_equiv*(numppPtbins+1)]));
        pPbCoMHistArr[pPbYbin]->Sumw2(); // instruct each histogram to propagate errors
    }


    /**** Iterate over each event ****/
    const int numentries = tree->GetEntries();
    const bool periodA = (thisRunNumber < 313500);

    double jpt, jeta, jphi, je, jy, eff, lumi, scale, numerator, denominator;
    int pbin, etabin, actetabin;
    TLorentzVector tlv;
    Trigger* bestTrigger = NULL;
    for (long long entry = 0; entry < numentries; entry++) {
        tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses
        
        if ((nvert == 0) || (nvert > 0 && vert_type->at(0) != 1)) continue; // Basic event selection: require a primary vertex

        for (int j = 0; j < njet; j++) {
            jpt = (double)jet_pt->at(j);
            jeta = (double)jet_eta->at(j);
            jphi = (double)jet_phi->at(j);
            je = (double)jet_e->at(j);

            // skip jets in the HEC region.
            if (lowerPhiCut < jphi && jphi < upperPhiCut && lowerEtaCut < jeta && jeta < upperEtaCut) continue;

            etabin = getEtabin(jeta);
            pbin = getPbin(jpt);
            if (pbin < 0 || etabin < 0 || pbin > numpbins || etabin > numetabins) continue; // this checks that the jets fall within the pt, eta bins

            bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins];
            if (bestTrigger == NULL || !bestTrigger->m_trig_bool) continue; // make sure we're not trying to look at a null trigger and that it actually fired, and that the jet met the minimum pt cut

            // calculate rescaling factors to obtain a cross-section
            eff = (*triggerEfficiencyFunctions)[bestTrigger->index]->Eval(jpt);
            if (periodA) actetabin = numetabins - etabin - 1;
            else actetabin = etabin; 
            lumi = kinematicLumiVec[pbin + actetabin*numpbins];
            if (eff == 0. || lumi == 0.) continue; // make sure we're not dividing by 0 for some reason

            scale = 1./(eff*lumi);
            scale *= etaPhiScaleFactors[getpPbYbin(jeta)];

            if (periodA) jeta *= -1.; // flips period A pseudorapidities into period B kinematics, so that when we shift into CoM frame we get the correct kinematics
            tlv.SetPtEtaPhiE(jpt, jeta, jphi, je);
            jy = tlv.Rapidity();

            pPbYbin = getpPbYbin(jy);
            if (pPbYbin < 0 || pPbYbin >= numpPbYbins) continue; // this checks that the jets fall within the pp eta bins, which cover a smaller range (so this is an important check)

            pPbCoMHistArr[pPbYbin]->Fill(jpt, scale);
        }
    }
    /**** End event iteration ****/


    /**** Write output histograms to a root file ****/
    TFile* output = new TFile(Form("%srun_%i.root", RpPbPath.c_str(), thisRunNumber), "RECREATE");
    for (pPbYbin = 0; pPbYbin < numpPbYbins; pPbYbin++) {
        pPbCoMHistArr[pPbYbin]->Scale(1e3/A); // divide by A to normalize the histogram to pp results, and scale by 1e3 to convert ub to nb
        pPbCoMHistArr[pPbYbin]->Write();
    }
//    TVectorD lum_vec(1);
//    lum_vec[0] = luminosity;
//    lum_vec.Write(Form("lum_vec_%i", thisRunNumber));

    TVectorD run_vec(6);
    run_vec[0] = thisRunNumber;
    run_vec[1] = numppYbins;
    run_vec[2] = numpPbYbins;
    run_vec[3] = numtrigs;
    run_vec[4] = numppPtbins;
    run_vec[5] = luminosity;
    run_vec.Write(Form("run_vec_%i", thisRunNumber));

    output->Close();
    /**** End write output ****/

    if (debugStatements) cout << "Status: In IdealRpPbAnalysis.C (breakpoint J): Finished calculating pt spectrum for run " << thisRunNumber << endl;
    return;
}
