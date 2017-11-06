#include "../triggerUtil.C"

void jets_xa_xp(int thisRunNumber, // Run number identifier.
                double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{
    if (skipRun(thisRunNumber)) return;
    initialize(thisRunNumber);

    const int numbins = 100;
    const int numhists = 2*numetabins;
    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    TTree* tree = (TTree*)(new TFile(Form("%srun_%i_raw.root", dataPath.c_str(), thisRunNumber)))->Get("tree");

    const double harr_scales[8] = {1, 1, 1, 1, 1, 1, 1, 1};   // rescaling factors so the histograms don't overlap

    const double* xbins = logspace(2e-4, 1.6, numbins);
    // Create an array of 16 histograms, one for each rapidity region and one for x_p, x_a. 
    TH1D* harr[numhists];
    for (int i = 0; i < numhists/2; i++) {
        harr[i] = new TH1D(Form("%ieta%i", thisRunNumber, i), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d^{2}#sigma/d#it{x}_{p} dy #left[pb#right]", etabins[i], etabins[i+1]), numbins, xbins);
        harr[i]->Sumw2();
    }
    for (int i = numhists/2; i < numhists; i++) {
        harr[i] = new TH1D(Form("%ieta%i", thisRunNumber, i), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d^{2}#sigma/d#it{x}_{a} dy #left[pb#right]", etabins[i%(numhists/2)], etabins[(i%(numhists/2))+1]), numbins, xbins);
        harr[i]->Sumw2();
    }
    TH2D* xaxpcorr = new TH2D(Form("xaxpcorr_run%i", thisRunNumber), ";#it{x}_{a};#it{x}_{p};d^{2}#sigma/d#it{x}_{p}d#it{x}_{a}", numbins, xbins, numbins, xbins);
//    cout << numtrigs << endl;

    // Create arrays to store trigger values for each event
    bool m_trig_bool[numtrigs];
    float m_trig_prescale[numtrigs];

    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    float j_phi[60] = {};
    float j_e[60] = {};
    int eventNumber = 0;
    int njet = 0;
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("j_e", j_e);
    tree->SetBranchAddress("njet", &njet);
    tree->SetBranchAddress("j_phi", j_phi);
    tree->SetBranchAddress("eventNumber", &eventNumber);

    // Set branch addresses
    for (Trigger* trig : trigger_vec) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &m_trig_bool[trig->index]);
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &m_trig_prescale[trig->index]);
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();

    double leadingjpt, leadingjeta, subleadingjpt, subleadingjeta, xp, xa, extra_jpt_sum;
    bool takeEvent;
    int ebin, pbin, index;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses

        leadingjpt = (double)j_pt[0];
        subleadingjpt = (double)j_pt[1];
        leadingjeta = (double)j_eta[0];
        subleadingjeta = (double)j_eta[1];
        for (int j = 0; j < njet; j++) {
            if (j_pt[j] > leadingjpt) {
                subleadingjpt = leadingjpt;
                subleadingjeta = leadingjeta;
                leadingjpt = (double)j_pt[j];
                leadingjeta = (double)j_eta[j];
            }
        }

        if (leadingjpt > 2000) cout << Form("High pt (%.0f GeV) jet detected in run %i, event %i!", leadingjpt, thisRunNumber, eventNumber) << endl;
       /* 
        extra_jpt_sum = 0;
        for (int j = 2; j < njet; j++) {
            extra_jpt_sum += j_pt[j];
        }
        takeEvent = extra_jpt_sum / (j_pt[0] + j_pt[1] + extra_jpt_sum) <= dijet_pt_frac_cutoff;
*/
        takeEvent = subleadingjpt/leadingjpt >= dijet_pt_ratio_cutoff;
        if (takeEvent) {
            
            ebin = 0;
            while (etabins[ebin] < leadingjeta) ebin++;
            ebin--;
            pbin = 0;
            while (pbins[pbin] < leadingjpt) pbin++;
            pbin--;
            if (pbin == -1 || pbin >= numpbins || ebin == -1 || ebin >= numetabins) continue;            

            index = best_bins[pbin + ebin*numpbins];
            takeEvent = m_trig_bool[index] && m_trig_prescale[index] > 0 && kinematic_lumi_vec[pbin + ebin*numpbins] > 0;
            if (!takeEvent) continue;
            xp = get_xp(leadingjpt, subleadingjpt, leadingjeta, subleadingjeta, periodA);
            xa = get_xa(leadingjpt, subleadingjpt, leadingjeta, subleadingjeta, periodA);
            double lumi = kinematic_lumi_vec[pbin + ebin*numpbins];

            if (periodA) {
                leadingjeta *= -1;
                ebin = 0;
                while (etabins[ebin] < leadingjeta) ebin++;
                ebin--;
            }
 
            harr[ebin]->Fill(xp, m_trig_prescale[index]/lumi);
            harr[ebin+numetabins]->Fill(xa, m_trig_prescale[index]/lumi);
            xaxpcorr->Fill(xa, xp, m_trig_prescale[index]/lumi);
        }
    }
    // Save to root file
    TFile* output = new TFile(Form("%sxdata/run_%i.root", rootPath.c_str(), thisRunNumber), "RECREATE");
    for (int i = 0; i < numhists; i++) {
        harr[i]->Scale(harr_scales[i%(numetabins)]/(A * (etabins[(i%numetabins)+1] - etabins[(i%numetabins)])), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        harr[i]->Write();
    }
    xaxpcorr->Write();
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    output->Close();
}
