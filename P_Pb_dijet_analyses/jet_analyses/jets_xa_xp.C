#include "../triggerUtil.C"

void jets_xa_xp(int thisRunNumber, // Run number identifier.
                double luminosity, // Integrated luminosity for this run. Presumed constant over the run period.
                bool periodA)
{
    initialize(thisRunNumber);

    const int numbins = 40;
    const int numhists = 2*numetabins;
    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    TTree* tree = (TTree*)(new TFile(Form("../rundata/run_%i_raw.root", thisRunNumber)))->Get("tree");

    const double harr_scales[8] = {1, 1, 1, 1, 1, 1, 1, 1};   // rescaling factors so the histograms don't overlap

    const double* xbins = logspace(0, 1.6, numbins);
    // Create an array of 16 histograms, one for each rapidity region and one for x_p, x_a. 
    TH1D* harr[numhists];
    for (int i = 0; i < numhists/2; i++) {
        harr[i] = new TH1D(Form("%ieta%i", thisRunNumber, i), Form("%1.1f < #eta < %1.1f;#it{x}_{p};d#sigma^{2}/d#it{x}_{p} dy #left[pb#right]", etabins[i], etabins[i+1]), numbins, xbins);
        harr[i]->Sumw2();
    }
    for (int i = numhists/2; i < numhists; i++) {
        harr[i] = new TH1D(Form("%ieta%i", thisRunNumber, i), Form("%1.1f < #eta < %1.1f;#it{x}_{a};d#sigma^{2}/d#it{x}_{a} dy #left[pb#right]", etabins[i%(numhists/2)], etabins[(i%(numhists/2))+1]), numbins, xbins);
        harr[i]->Sumw2();
    }
//    cout << numtrigs << endl;

    // Create arrays to store trigger values for each event
    bool m_trig_bool[numtrigs];
    float m_trig_prescale[numtrigs];

    // Create arrays to store jet data for each event
    float j_pt[60] = {};
    float j_eta[60] = {};
    float j_phi[60] = {};
    float j_e[60] = {};
    int njet = 0;
    tree->SetBranchAddress("j_pt", j_pt);
    tree->SetBranchAddress("j_eta", j_eta);
    tree->SetBranchAddress("j_e", j_e);
    tree->SetBranchAddress("njet", &njet);
    tree->SetBranchAddress("j_phi", j_phi);

    // Set branch addresses
    for (Trigger* trig : trigger_vec) {
        tree->SetBranchAddress(Form("%s", trig->name.c_str()), &m_trig_bool[trig->index]);
        tree->SetBranchAddress(Form("%s_prescale", trig->name.c_str()), &m_trig_prescale[trig->index]);
    }

    // Iterate over each event
    const int numentries = tree->GetEntries();

    double leading_jpt, leading_jeta, jeta0, jpt0, jeta1, jpt1, xp, xa, extra_jpt_sum;
    bool takeEvent;
    int ebin, pbin, index;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses

        extra_jpt_sum = 0;
        for (int j = 2; j < njet; j++) {
            extra_jpt_sum += j_pt[j];
        }
        takeEvent = extra_jpt_sum / (j_pt[0] + j_pt[1] + extra_jpt_sum) <= dijet_pt_frac_cutoff;
        if (takeEvent) {
            
            jpt0 = (double)j_pt[0];
            jpt1 = (double)j_pt[1];
            
            jeta0 = (double)j_eta[0];
            jeta1 = (double)j_eta[1];

            if (jpt0 > jpt1) {
                leading_jpt = jpt0;
                leading_jeta = jeta0;
            } else {
                leading_jpt = jpt1;
                leading_jeta = jeta1;
            }

            ebin = 0;
            while (etabins[ebin] < leading_jeta) ebin++;
            ebin--;
            pbin = 0;
            while (pbins[pbin] < leading_jpt) pbin++;
            pbin--;
            if (pbin == -1 || pbin >= numpbins || ebin == -1 || ebin >= numetabins) continue;            

            index = (*best_trig_indices)[pbin + ebin*numpbins];
            takeEvent = m_trig_bool[index] && m_trig_prescale[index] > 0 && (*integrated_luminosity_vec)[index + ebin*numtrigs] > 0;
            if (!takeEvent) continue;

            xp = get_xp(jpt0, jpt1, jeta0, jeta1, periodA);
            xa = get_xa(jpt0, jpt1, jeta0, jeta1, periodA);

            if (!periodA) {
                harr[ebin]->Fill(xp, m_trig_prescale[index]/(*integrated_luminosity_vec)[index + ebin*numtrigs]);
                harr[ebin+numetabins]->Fill(xa, m_trig_prescale[index]/(*integrated_luminosity_vec)[index + ebin*numtrigs]);
            } else {
                harr[numetabins-1-ebin]->Fill(xp, m_trig_prescale[index]/(*integrated_luminosity_vec)[index + ebin*numtrigs]);
                harr[numhists-1-ebin]->Fill(xa, m_trig_prescale[index]/(*integrated_luminosity_vec)[index + ebin*numtrigs]);
            }
        }
    }
    // Save to root file
    TFile* output = new TFile(Form("../rootFiles/xdata/run_%i.root", thisRunNumber), "RECREATE");
    for (int i = 0; i < numhists; i++) {
        harr[i]->Scale(harr_scales[i%(numetabins)]/(A * (etabins[(i%numetabins)+1] - etabins[(i%numetabins)])), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        harr[i]->Write();
    }
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    output->Close();
}
