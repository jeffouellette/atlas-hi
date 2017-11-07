#include "../triggerUtil.C"

void jets_pt(const int thisRunNumber, // Run number identifier.
             double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{

    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    initialize(thisRunNumber, true);
    const int numhists = numetabins;

    TTree* tree = (TTree*)(new TFile(Form("%srun_%i_raw.root", dataPath.c_str(), thisRunNumber)))->Get("tree");

    //Useful arrays for binning and cutting jets between bins
    const double harr_scales[numhists] = {0.005, 0.03, 0.1, 0.5, 1, 0.3, 0.05, 0.01};   // rescaling factors so the histograms don't overlap

    // Create branching addresses:  
    // Create arrays to store trigger values for each event
    bool m_trig_bool[numtrigs];   // stores whether trigger was triggered
    float m_trig_prescale[numtrigs];      // stores the prescaling factor for the trigger
    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    int njet = 0;

    // Set branch addresses
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("njet", &njet);
    for (Trigger* trig : trigger_vec) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &m_trig_bool[trig->index]);
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &m_trig_prescale[trig->index]);
    }

    TH1D* harr[numhists];
    for (int i = 0; i < numhists; i++) {
        TString histname = Form("%ieta%i", thisRunNumber, i);
        harr[i] = new TH1D(histname, Form("%g < #eta < %g (#times %g);#it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", etabins[i], etabins[i+1], harr_scales[i]), numpbins, pbins);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }


    // Iterate over each event
    const int numentries = tree->GetEntries();

    double jpt;
    bool takeEvent;
    Trigger* trig;
    int pbin, ebin, index;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses
        
        for (int j = 0; j < njet; j++) { // iterate over each jet momentum
            jpt = (double)j_pt[j];

            ebin = 0;
            while (etabins[ebin] <= j_eta[j]) ebin++;
            ebin--;
            pbin = 0;
            while (pbins[pbin] <= jpt) pbin++;
            pbin--;
            if (pbin == -1 || pbin >= numpbins || ebin == -1 || ebin >= numetabins) continue; // Verify that the jet momentum and pseudorapidity make sense.
 
            trig = (*((*trigger_pt_eta_bin_map)[ebin]))[pbin];
            index = trig->index;
//            if (thisRunNumber == 313063 && jpt > 300) cout << Form("jpt=%f, triggered=%i, prescale=%f", jpt, m_trig_bool[index], m_trig_prescale[index]) << endl;
            //takeEvent = (enabledTriggers->GetBinContent(pbin+1, ebin+1) == 1);
            takeEvent = trig->enabled;
            takeEvent = takeEvent && m_trig_bool[index] && m_trig_prescale[index] > 0; // If the trigger for this bin is enabled, we know the run must have been sensitive to this bin.
            if (takeEvent) {
                harr[ebin]->Fill(jpt, m_trig_prescale[index]);
     //           if (thisRunNumber == 313063 && jpt > 300) cout << Form("Plotting jpt = %f...", jpt) << endl;
            }
        }
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("./pt_data/run_%i.root", thisRunNumber), "RECREATE");
    for (int i = 0; i< numhists; i++) {
        harr[i]->Scale(1/A, "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        harr[i]->Write();
    }
    TH2D* luminosity_hist = new TH2D(Form("luminosity_hist_run%i", thisRunNumber), "", numpbins, -0.5, numpbins+0.5, numetabins, -0.5, numetabins+0.5);
    for (int pbin = 0; pbin < numpbins; pbin++) {
        for (int ebin = 0; ebin < numetabins; ebin++) {
            if ((*((*trigger_pt_eta_bin_map)[ebin]))[pbin]->enabled) luminosity_hist->SetBinContent(pbin+1, ebin+1, luminosity);
            //if (enabledTriggers->GetBinContent(pbin+1, ebin+1) == 1) luminosity_hist->SetBinContent(pbin+1, ebin+1, luminosity);
            else luminosity_hist->SetBinContent(pbin+1, ebin+1, 0.);
        }
    }
    luminosity_hist->Write();
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    output->Close();
}
