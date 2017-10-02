#include "triggerUtil.C"

void jets_pt(const int runNumber, // Run number identifier.
             double luminosity) // Integrated luminosity for this run. Presumed constant over the run period.
{

    const int numhists = 8;
    luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

    initialize();

    TTree* tree = (TTree*)(new TFile(Form("./rundata/run_%i_raw.root", runNumber)))->Get("tree");

    //Useful arrays for binning and cutting jets between bins
    const double xbins[42] = {25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 220., 240., 260., 280., 300., 350., 400., 500., 600., 800., 1100., 1500., 2000., 2500., 6000.};
//    const double xbins[41] = {25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 350, 400, 500, 600, 800, 1100, 1500, 2500, 6000};
    const int numbins = sizeof(xbins)/sizeof(xbins[0]) - 1;
    /*const double xbins_n200eta490[16] = {25, 40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};           // Eta range: -4.9 < eta <= -2
    const double xbins_0eta200[16] = {25, 40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};              // Eta range: -2 < eta < 2
    const double xbins_p200eta320[18] = {25, 40, 50, 55, 60, 70, 75, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};   // Eta range: 2 <= eta < 3.2
    const double xbins_p320eta490[17] = {25, 40, 50, 60, 65, 70, 85, 110, 150, 200, 280, 400, 600, 850, 1100, 2000, 6000};   // Eta range: 3.2 <= eta < 4.9
    const double* xbins[numhists] = {xbins_n200eta490, xbins_n200eta490, xbins_0eta200, xbins_0eta200, xbins_0eta200, xbins_0eta200, xbins_p200eta320, xbins_p320eta490};
    const int len_xbins[numhists] = {16, 16, 16, 16, 16, 16, 18, 17};
*/
    const double d_eta[numhists] = {1.7, 1.2, 1, 1, 1, 1, 1.2, 1.7};
    const double eta_cuts[numhists+1] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
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
        harr[i] = new TH1D(Form("%ieta%i", runNumber, i), Form("%g < #eta < %g (#times %g); #it{p}_{T}^{jet} #left[GeV/#it{c}#right];d^{2}#sigma/d#it{p}_{T}dy #left[pb (GeV/#it{c})^{-1}#right]", eta_cuts[i], eta_cuts[i+1], harr_scales[i]), numbins, xbins);
        harr[i]->Sumw2(); // instruct each histogram to propagate errors
    }


    // Iterate over each event
    const int numentries = tree->GetEntries();

    double jeta, jpt;
    bool takeEvent;
    for (int i = 0; i < numentries; i++) {
        tree->GetEntry(i); // stores trigger values and data in the designated branch addresses
        for (int j = 0; j < njet; j++) { // iterate over each jet momentum
            jeta = (double)j_eta[j];
            jpt = (double)j_pt[j];
            if (0 < jeta && jeta < 1) {
                for (Trigger* trig : triggers_p0eta100) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0; // Only take the event if triggered and trigger was not disabled
                    if (takeEvent && trig->min_pt <= jpt && jpt < trig->max_pt) {
                        harr[4]->Fill((double)jpt, m_trig_prescale[trig->index]);
                        break; // Break to ensure that only one trigger allows the event to be recorded
                    }
                }
            }
            else if (-1 < jeta && jeta < 0) {
                for (Trigger* trig : triggers_n100eta0) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt && jpt < trig->max_pt) {
                        harr[3]->Fill((double)jpt, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (1 < jeta && jeta < 2) {
                for (Trigger* trig : triggers_p100eta200) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt && jpt < trig->max_pt) {
                        harr[5]->Fill((double)jpt, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (-2 < jeta && jeta < -1) {
                for (Trigger* trig : triggers_n200eta100) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt && jpt < trig->max_pt) {
                        harr[2]->Fill((double)jpt, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (2 < jeta && jeta < 3.2) {
                for (Trigger* trig : triggers_p200eta320) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt && jpt < trig->max_pt) {
                        harr[6]->Fill((double)jpt, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (-3.2 < jeta && jeta < -2) {
                for (Trigger* trig : triggers_n320eta200) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt && jpt < trig->max_pt) {
                        harr[1]->Fill((double)jpt, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (3.2 < jeta && jeta < 4.9) {
                for (Trigger* trig : triggers_p320eta490) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt && jpt < trig->max_pt) {
                        harr[7]->Fill((double)jpt, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
            else if (-4.9 < jeta && jeta < -3.2) {
                for (Trigger* trig : triggers_n490eta320) {
                    takeEvent = m_trig_bool[trig->index] && m_trig_prescale[trig->index] > 0;
                    if (takeEvent && trig->min_pt <= jpt && jpt < trig->max_pt) {
                        harr[0]->Fill((double)jpt, m_trig_prescale[trig->index]);
                        break;
                    }
                }
            }
        }
    }

    // Write histograms to a root file
    TFile* output = new TFile(Form("./pt_data/run_%i.root", runNumber), "RECREATE");
    for (int i = 0; i< numhists; i++) {
        harr[i]->Scale((harr_scales[i]) / (A * luminosity * d_eta[i]), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
        harr[i]->Write();
    }
    TVectorD lum_vec(1);
    lum_vec[0] = luminosity;
    lum_vec.Write("lum_vec");
    output->Close();
}
