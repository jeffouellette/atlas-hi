#include "../triggerUtil.C"

const double etaCoM = -0.465; // boost into CoM frame in period B kinematics
const int numppPtbins = 34;
const int numppEtabins = 6;
const int numpPbCoMEtabins = 18;
double ppPtbins[(numppPtbins+1)] = {}; //= new double[(numppPtbins+1)*numppEtabins];
double ppEtabins[numppEtabins+1] = {0, 0.5, 1, 1.5, 2, 2.5, 3};
int numppPtEtabins[numppEtabins] = {};
double pPbCoMEtabins[numpPbCoMEtabins+1] = {-3.465, -3.2, -2.965, -2.465, -2, -1.965, -1.465, -1, -0.965, -0.465, 0, 0.035, 0.535, 1, 1.035, 1.535, 2, 2.035, 2.535};
double* ppJetSpectrum = new double[numppPtbins*numppEtabins];
double* ppJetSpectrumStatError = new double[numppPtbins*numppEtabins];

int getppPbin (double pt) {
    int bin = 0;
    while (bin < numppPtbins+1 && ppPtbins[bin++] < pt);
    return bin - 2;
}


int getppEtabin (double eta) {
    int bin = 0;
    while (bin < numppEtabins+1 && ppEtabins[bin++] < eta);
    return bin - 2;
}


int getpPbCoMEtabin (double eta) {
    int bin = 0;
    while (bin < numpPbCoMEtabins+1 && pPbCoMEtabins[bin++] < eta);
    return bin - 2;
}

/**
 * Loads binning information from the 8TeV pp reference data.
 */
TH1D** setupPPConfiguration() {

    for (int ppPtbin = 0; ppPtbin < numppPtbins; ppPtbin++) {
        ppPtbins[ppPtbin] = 0;
        for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) ppJetSpectrum[ppPtbin + ppEtabin*numppPtbins] = 0;
    }

    string path = workPath + "HepData (#28325609)/";
    string dummyline;
    double pt_center, pt_left, pt_right, dsigma, stat_err;
    int ppPtbin, ppEtabin;
    vector<string> files = {"atlas_2012_jet_antiktr04_incljetpt_eta1_highmu_data.txt", "atlas_2012_jet_antiktr04_incljetpt_eta2_highmu_data.txt", "atlas_2012_jet_antiktr04_incljetpt_eta3_highmu_data.txt", "atlas_2012_jet_antiktr04_incljetpt_eta4_highmu_data.txt", "atlas_2012_jet_antiktr04_incljetpt_eta5_highmu_data.txt", "atlas_2012_jet_antiktr04_incljetpt_eta6_highmu_data.txt"};

    ppEtabin = 0;
    for (string file : files) {
        ifstream dataStream;
        dataStream.open(path + file);
        if (!dataStream.is_open()) {
            cout << "Error: In IdealRpPbAnalysis.C (breakpoint A): Cannot find input file!" << endl;
            throw runtime_error("File not found");
        }
        for (int i = 0; i < 641; i++) getline(dataStream, dummyline); // first 640 lines are systematic errors so skip them, 641st line is a table layout
        ppPtbin = 0;
        while (dataStream) {
            dataStream >> pt_center;
            dataStream >> pt_left;
            dataStream >> pt_right;
            dataStream >> dsigma;
            dataStream >> stat_err;
            //ppPtbins[ppPtbin] = pt_left;
            if (file == files[0]) {
                if (ppPtbin > 0) assert(ppPtbins[(ppPtbin-1)] == pt_left);
                else ppPtbins[ppPtbin] = pt_left;
                ppPtbins[(ppPtbin+1)] = pt_right;
            }
            numppPtEtabins[ppEtabin]++;
            
            
            ppJetSpectrum[ppPtbin + ppEtabin*numppPtbins] = dsigma;
            ppJetSpectrumStatError[ppPtbin + ppEtabin*numppPtbins] = stat_err*dsigma/100.;
            ppPtbin++;
        }
        ppEtabin++;

        dataStream.close();
    }

    TH1D** ppHistArr = new TH1D*[numppEtabins];
    for (ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) {
        TString histName = Form("pp_spectrum_ppEtabin%i", ppEtabin);
        ppHistArr[ppEtabin] = new TH1D (histName, ";#it{p}_{T}^{jet} #left[GeV#right];d#sigma/Ad#it{p}_{T} #left[nb GeV^{-1}#right]", numppPtEtabins[ppEtabin], ppPtbins);
        for (ppPtbin = 0; ppPtbin < numppPtbins; ppPtbin++) {
            ppHistArr[ppEtabin]->SetBinContent(ppPtbin+1, 1e-3*ppJetSpectrum[ppPtbin + ppEtabin*numppPtbins]); // 1e-3 to convert pb to nb
            ppHistArr[ppEtabin]->SetBinError(ppPtbin+1, 1e-3*ppJetSpectrumStatError[ppPtbin + ppEtabin*numppPtbins]);
        }
    }
    return ppHistArr; 
}


void IdealRpPbAnalysis(const int thisRunNumber, // Run number identifier.
                       double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;

    TH1D** ppHistArr = setupPPConfiguration();
    for (int ppEtabin = 0; ppEtabin < numppEtabins; ppEtabin++) delete ppHistArr[ppEtabin]; // don't really need these right now so just delete them

//    if (thisRunNumber == 313063) for (int ppPtbin = 0; ppPtbin < numppPtbins; ppPtbin++) cout << "bin " << ppPtbin << ": " << ppPtbins[ppPtbin] << ", " << ppPtbins[ppPtbin+1] << endl;

    initialize(thisRunNumber, true);
    vector<TF1*>* triggerEfficiencyFunctions = getTriggerEfficiencyFunctions();
    
    /**** Generate list of physics triggers ****/
    vector<Trigger*>* triggerSubvector = getTriggerSubvector(thisRunNumber);
    if (debugStatements) {
        cout << "Status: In IdealRpPbAnalysis.C (breakpoint B): Processing run " << thisRunNumber << " with triggers:" << endl;
        for (Trigger* trig : (*triggerSubvector)) {
            cout << "\t" << trig->name << endl;
        }
    }
    /**** End generate list of physics triggers ****/

    const int numhists = numtrigs * numppEtabins;

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
                    if (debugStatements) cout << "Status: In IdealRpPbAnalysis.C (breakpoint C): Found " << fname.Data() << endl; 
                    if (fname.Contains(to_string(thisRunNumber))) {
                        tree = (TTree*)(new TFile(dataPath+fname, "READ"))->Get("tree");
                        break;
                    }
                }
            }
        }
    }
    if (tree == NULL) {
        cout << "Error: In IdealRpPbAnalysis.C (breakpoint D): TTree not obtained for given run number. Quitting." << endl;
        return;
    }
    /**** End find TTree ****/


    /**** Disable loading of unimportant branch values - speeds up entry retrieval ****/
    {
        vector<string> interestingBranchNames = {"njet", "j_pt", "j_eta", "j_phi", "j_e", "vert_type", "nvert"};
        TObjArray* branches = (TObjArray*)(tree->GetListOfBranches());
        bool interestingBranch;
        for (TObject* obj : *branches) {
            TString branchName = (TString)obj->GetName();
            if (debugStatements) cout << "Status: In IdealRpPbAnalysis.C (breakpoint E): Tree contains branch \"" << branchName.Data() << "\"" << endl;
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
    // Create branching addresses:  
    // Create arrays to store trigger values for each event
    //bool m_trig_bool[numtrigs];   // stores whether trigger was triggered
    //float m_trig_prescale[numtrigs];      // stores the prescaling factor for the trigger
    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    float j_phi[60] = {};
    float j_e[60] = {};
    int njet = 0;
    int nvert = 0;
    int vert_type[60] = {};

    // Set branch addresses
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("j_phi", j_phi);
    tree->SetBranchAddress("j_e", j_e);
    tree->SetBranchAddress("njet", &njet);
    tree->SetBranchAddress("nvert", &nvert);
    tree->SetBranchAddress("vert_type", vert_type);
    for (Trigger* trig : (*triggerSubvector)) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &(trig->m_trig_bool));
    }
    /**** End set branch addresses ****/


    /**** Histogram initialization ****/
    TH1D* pPbCoMHistArr[numhists];
    int pPbCoMEtabin, index;
    for (pPbCoMEtabin = 0; pPbCoMEtabin < numpPbCoMEtabins; pPbCoMEtabin++) {
        TString histName = Form("pPb_spectrum_run%i_pPbCoMEtabin%i", thisRunNumber, pPbCoMEtabin);
        int ppEtabin_equiv = getppEtabin(TMath::Abs(0.5*(pPbCoMEtabins[pPbCoMEtabin]+pPbCoMEtabins[pPbCoMEtabin+1]) - etaCoM));
        pPbCoMHistArr[pPbCoMEtabin] = new TH1D(histName, ";#it{p}_{T}^{jet} #left[GeV#right];d#sigma/Ad#it{p}_{T} #left[nb GeV^{-1}#right]", numppPtEtabins[ppEtabin_equiv], ppPtbins);
        pPbCoMHistArr[pPbCoMEtabin]->Sumw2(); // instruct each histogram to propagate errors
    }


    /**** Iterate over each event ****/
    const int numentries = tree->GetEntries();
    const bool periodA = (thisRunNumber < 313500);

    double jpt, jeta, jphi, je, jy, eff, lumi, scale;
    int pbin, etabin, actetabin;
    TLorentzVector tlv;
    Trigger* bestTrigger = NULL;
    for (long long entry = 0; entry < numentries; entry++) {
        tree->GetEntry(entry); // stores trigger values and data in the designated branch addresses
        
        if ((nvert == 0) || (nvert > 0 && vert_type[0] != 1)) continue; // Basic event selection: require a primary vertex

        for (int j = 0; j < njet; j++) {
            jpt = (double)j_pt[j];
            jeta = (double)j_eta[j];
            jphi = (double)j_phi[j];
            je = (double)j_e[j];

            etabin = getEtabin(jeta);
            pbin = getPbin(jpt);
            if (pbin < 0 || etabin < 0 || pbin > numpbins || etabin > numetabins) continue; // this checks that the jets fall within the pt, eta bins

            bestTrigger = kinematicTriggerVec[pbin + etabin*numpbins];
            if (bestTrigger == NULL || !bestTrigger->m_trig_bool) continue; // make sure we're not trying to look at a null trigger and that it actually fired, and that the jet met the minimum pt cut

            eff = (*triggerEfficiencyFunctions)[bestTrigger->index]->Eval(jpt);
            if (periodA) actetabin = numetabins - etabin - 1;
            else actetabin = etabin; 
            lumi = kinematicLumiVec[pbin + actetabin*numpbins];
            if (eff == 0. || lumi == 0.) continue; // make sure we're not dividing by 0 for some reason

            if (periodA) jy *= -1.; // flips period A pseudorapidities into period B kinematics, so that when we shift into CoM frame we get the correct kinematics
            tlv.SetPtEtaPhiE(jpt, jeta, jphi, je);
            jy = tlv.Rapidity();

            pPbCoMEtabin = getpPbCoMEtabin(jy);
            if (pPbCoMEtabin < 0 || pPbCoMEtabin > numpPbCoMEtabins) continue; // this checks that the jets fall within the pp eta bins, which cover a smaller range (so this is an important check)

            scale = ((2*pi)/(2*pi- (upperPhiCut - lowerPhiCut)))*(1./(eff*lumi));
            
            if (jphi <= lowerPhiCut || jphi >= upperPhiCut) pPbCoMHistArr[pPbCoMEtabin]->Fill(jpt, scale);
            //pPbCoMHistArr[pPbCoMEtabin]->Fill(jpt, 1./(eff*lumi));
        }
    }
    /**** End event iteration ****/


    /**** Write output histograms to a root file ****/
    TFile* output = new TFile(Form("%srun_%i.root", RpPbPath.c_str(), thisRunNumber), "RECREATE");
    for (pPbCoMEtabin = 0; pPbCoMEtabin < numpPbCoMEtabins; pPbCoMEtabin++) {
        pPbCoMHistArr[pPbCoMEtabin]->Scale(1e3/A); // divide by A to normalize the histogram to pp results, and scale by 1e3 to convert ub to nb
        pPbCoMHistArr[pPbCoMEtabin]->Write();
    }
//    TVectorD lum_vec(1);
//    lum_vec[0] = luminosity;
//    lum_vec.Write(Form("lum_vec_%i", thisRunNumber));

    TVectorD run_vec(6);
    run_vec[0] = thisRunNumber;
    run_vec[1] = numppEtabins;
    run_vec[2] = numpPbCoMEtabins;
    run_vec[3] = numtrigs;
    run_vec[4] = numppPtbins;
    run_vec[5] = luminosity;
    run_vec.Write(Form("run_vec_%i", thisRunNumber));

    output->Close();
    /**** End write output ****/

    if (debugStatements) cout << "Status: In IdealRpPbAnalysis.C (breakpoint F): Finished calculating pt spectrum for run " << thisRunNumber << endl;
    return;
}
