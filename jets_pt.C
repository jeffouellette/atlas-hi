
void jets_pt(int runNumber, float luminosity) {

        luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

	float A = 208; // value of A for Pb

        /* ETA RANGES:
                -4.9 < eta <= -3.2
                -3.2 < eta <= -2
                -2 < eta <= -1
                -1 < eta <= 0
                0 <= eta < 1
                1 <= eta < 2
                2 <= eta < 3.2
                3.2 <= eta < 4.9
        */

	// Store trigger names as an array of strings and loop over them to set the branch addresses

	static const int trigLength = 24;
	const char* m_trig_string[trigLength] = {
                "HLT_j15_p320eta490_L1MBTS_1_1",
                "HLT_j15_ion_p320eta490_L1MBTS_1_1",
                "HLT_j30_0eta490_L1TE10",
                "HLT_j30_ion_0eta490_L1TE10",
                "HLT_j40_L1J5",
                "HLT_j40_ion_L1J5",
                "HLT_j45_p200eta320",
                "HLT_j45_ion_p200eta320",
                "HLT_j50_L1J10",
                "HLT_j50_ion_L1J10",
                "HLT_j55_p320eta490",
                "HLT_j55_ion_p320eta490",
                "HLT_j60",
                "HLT_j60_ion_L1J20",
                "HLT_j65_p200eta320",
                "HLT_j65_ion_p200eta320",
                "HLT_j75_L1J20",
                "HLT_j75_ion_L1J20",
                "HLT_j100_L1J20",
                "HLT_j100_ion_L1J20",
                "HLT_2j10_p320eta490_L1TE10",
                "HLT_2j10_ion_p320eta490_L1TE10",
                "HLT_2j30_p320eta490",
                "HLT_2j30_ion_p320eta490"
        };


        // Initialize maps of trigger numbers (as ordered above) to the appropriate jet cutoffs for a specific eta range. This is done to obtain continuous coverage over the p_t spectrum.

        std::map<int, float> trig_lower_n320eta490; // Eta range : -4.9 < eta <= -3.2
        std::map<int, float> trig_upper_n320eta490;
        trig_lower_n320eta490[3] = 40;
        trig_upper_n320eta490[3] = 50;
        trig_lower_n320eta490[5] = 50;
        trig_upper_n320eta490[5] = 60;
        trig_lower_n320eta490[9] = 60;
        trig_upper_n320eta490[9] = 70;
        trig_lower_n320eta490[13] = 70;
        trig_upper_n320eta490[13] = 85;
        trig_lower_n320eta490[17] = 85;
        trig_upper_n320eta490[17] = 110;
        trig_lower_n320eta490[19] = 110;
        trig_upper_n320eta490[19] = 2000;
        std::map<int, float> trig_lower_n200eta320; // Eta range: -3.2 < eta <= -2
        std::map<int, float> trig_upper_n200eta320;
        trig_lower_n200eta320[3] = 40;
        trig_upper_n200eta320[3] = 50;
        trig_lower_n200eta320[5] = 50;
        trig_upper_n200eta320[5] = 60;
        trig_lower_n200eta320[9] = 60;
        trig_upper_n200eta320[9] = 70;
        trig_lower_n200eta320[13] = 70;
        trig_upper_n200eta320[13] = 85;
        trig_lower_n200eta320[17] = 85;
        trig_upper_n200eta320[17] = 110;
        trig_lower_n200eta320[19] = 110;
        trig_upper_n200eta320[19] = 2000;
        std::map<int, float> trig_lower_0eta200; // Eta range -2 < eta < 2
        std::map<int, float> trig_upper_0eta200;
        trig_lower_0eta200[3] = 40;
        trig_upper_0eta200[3] = 50;
        trig_lower_0eta200[5] = 50;
        trig_upper_0eta200[5] = 60;
        trig_lower_0eta200[9] = 60;
        trig_upper_0eta200[9] = 70;
        trig_lower_0eta200[13] = 70;
        trig_upper_0eta200[13] = 85;
        trig_lower_0eta200[17] = 85;
        trig_upper_0eta200[17] = 110;
        trig_lower_0eta200[19] = 110;
        trig_upper_0eta200[19] = 2000;
        std::map<int, float> trig_lower_p200eta320; // Eta range: 2 <= eta < 3.2
        std::map<int, float> trig_upper_p200eta320;
        trig_lower_p200eta320[3] = 40;
        trig_upper_p200eta320[3] = 50;
        trig_lower_p200eta320[5] = 50;
        trig_upper_p200eta320[5] = 55;
        trig_lower_p200eta320[7] = 55;
        trig_upper_p200eta320[7] = 60;
        trig_lower_p200eta320[9] = 60;
        trig_upper_p200eta320[9] = 70;
        trig_lower_p200eta320[13] = 70;
        trig_upper_p200eta320[13] = 75;
        trig_lower_p200eta320[15] = 75;
        trig_upper_p200eta320[15] = 85;
        trig_lower_p200eta320[17] = 85;
        trig_upper_p200eta320[17] = 110;
        trig_lower_p200eta320[19] = 110;
        trig_upper_p200eta320[19] = 2000;
        std::map<int, float> trig_lower_p320eta490; // Eta range: 3.2 <= eta < 4.9
        std::map<int, float> trig_upper_p320eta490;
        trig_lower_p320eta490[1] = 25;
        trig_upper_p320eta490[1] = 40;
        trig_lower_p320eta490[3] = 40;
        trig_upper_p320eta490[3] = 50;
        trig_lower_p320eta490[5] = 50;
        trig_upper_p320eta490[5] = 60;
        trig_lower_p320eta490[9] = 60;
        trig_upper_p320eta490[9] = 65;
        trig_lower_p320eta490[11] = 65;
        trig_upper_p320eta490[11] = 70;
        trig_lower_p320eta490[13] = 70;
        trig_upper_p320eta490[13] = 85;
        trig_lower_p320eta490[17] = 85;
        trig_upper_p320eta490[17] = 110;
        trig_lower_p320eta490[19] = 110;
        trig_upper_p320eta490[19] = 2000;

        // End trigger-jet cutoff map initialization



        TTree* tree = (TTree*)(new TFile(Form("./rundata/run_%i_raw.root", runNumber)))->Get("tree");

        //Useful arrays for binning and cutting jets between bins
	const float xbins_n200eta490[13] = {40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 1000, 2000};              // Eta range: -4.9 < eta <= -2
        const float jet_cuts_n200eta490[7] = {40, 50, 60, 70, 85, 110, 2000};
        const float xbins_0eta200[13] = {40, 50, 60, 70, 85, 110, 150, 200, 280, 400, 600, 1000, 2000};                 // Eta range: -2 < eta < 2
	const float jet_cuts_0eta200[7] = {40, 50, 60, 70, 85, 110, 2000};
        const float xbins_p200eta320[15] = {40, 50, 55, 60, 70, 75, 85, 110, 150, 200, 280, 400, 600, 1000, 2000};	// Eta range: 2 <= eta < 3.2
        const float jet_cuts_p200eta320[9] = {40, 50, 55, 60, 70, 75, 85, 110, 2000};
        const float xbins_p320eta490[15] = {25, 40, 50, 60, 65, 70, 85, 110, 150, 200, 280, 400, 600, 1000, 2000};      // Eta range: 3.2 <= eta < 4.9
        const float jet_cuts_p320eta490[9] = {25, 40, 50, 60, 65, 70, 85, 110, 2000};

        const float d_eta[8] = {1.7, 1.2, 1, 1, 1, 1, 1.2, 1.7};
	const float eta_cuts[9] = {-4.9, -3.2, -2, -1, 0, 1, 2, 3.2, 4.9};  // cuts for each eta range
	const float harr_scales[8] = {0.005, 0.03, 0.1, 0.5, 1, 0.3, 0.05, 0.01};   // rescaling factors so the histograms don't overlap

        // Create branching addresses:	
	// Create arrays to store trigger values for each event
	bool m_trig_bool[trigLength];   // stores whether trigger was triggered
	float m_trig_prescale[trigLength];      // stores the prescaling factor for the trigger
	// Create arrays to store jet data for each event
	float j_pt[60] = {};
	float j_eta[60] = {};

        // Set branch addresses
	tree->SetBranchAddress("j_pt", j_pt);
	tree->SetBranchAddress("j_eta", j_eta);
	for (int i = 0; i < trigLength; i++) {
		tree->SetBranchAddress(m_trig_string[i], &m_trig_bool[i]); 
		tree->SetBranchAddress(Form("%s_prescale", m_trig_string[i]), &m_trig_prescale[i]);
	}
        
        // Create an array of 9 histograms, one for each rapidity region.
	TH1F* harr[8];
	harr[0] = new TH1F(Form("%ieta0", runNumber), "-4.9 < #eta < -3.2 (#times 0.005);#it{p}_{T}^{jet} [GeV];d#sigma/d#it{p}_{T} [pb/GeV]", sizeof(xbins_n200eta490)/sizeof(xbins_n200eta490[0])-1, xbins_n200eta490);
	harr[1] = new TH1F(Form("%ieta1", runNumber), "-3.2 < #eta < -2 (#times 0.03)", sizeof(xbins_n200eta490)/sizeof(xbins_n200eta490[0])-1, xbins_n200eta490);
	harr[2] = new TH1F(Form("%ieta2", runNumber), "-2 < #eta < -1 (#times 0.1)", sizeof(xbins_0eta200)/sizeof(xbins_0eta200[0])-1, xbins_0eta200);
	harr[3] = new TH1F(Form("%ieta3", runNumber), "-1 < #eta < 0 (#times 0.5)", sizeof(xbins_0eta200)/sizeof(xbins_0eta200[0])-1, xbins_0eta200);
	harr[4] = new TH1F(Form("%ieta4", runNumber), "0 < #eta < 1 (#times 1)", sizeof(xbins_0eta200)/sizeof(xbins_0eta200[0])-1, xbins_0eta200);
	harr[5] = new TH1F(Form("%ieta5", runNumber), "1 < #eta < 2 (#times 0.3)", sizeof(xbins_0eta200)/sizeof(xbins_0eta200[0])-1, xbins_0eta200);
	harr[6] = new TH1F(Form("%ieta6", runNumber), "2 < #eta < 3.2 (#times 0.05)", sizeof(xbins_p200eta320)/sizeof(xbins_p200eta320[0])-1, xbins_p200eta320);
	harr[7] = new TH1F(Form("%ieta7", runNumber), "3.2 < #eta < 4.9 (#times 0.01)", sizeof(xbins_p320eta490)/sizeof(xbins_p320eta490[0])-1, xbins_p320eta490);
	for (int i = 0; i< sizeof(harr)/sizeof(harr[0]); i++) {
		harr[i]->Sumw2(); // instruct each histogram to propagate errors
	}


	// Iterate over each event
	int numentries = tree->GetEntries();

        std::vector<int> trig_n200eta490 = {3, 5, 9, 13, 17, 19};
        std::vector<int> trig_0eta200 = {3, 5, 9, 13, 17, 19};
        std::vector<int> trig_p200eta320 = {3, 5, 7, 9, 13, 15, 17, 19};
        std::vector<int> trig_p320eta490 = {1, 3, 5, 9, 11, 13, 17, 19};
        float jeta;
        for (int i = 0; i < numentries; i++) {
		tree->GetEntry(i); // stores trigger values and data in the designated branch addresses
                int j = 0;
		while (j_pt[j] != 0) { // iterate over each jet momentum
                        jeta = j_eta[j];

                        if (TMath::Abs(jeta) < 2 && TMath::Abs(jeta) >= 0) {
		                for (int trig_num : trig_0eta200) { // iterate over each trigger
		        		if (m_trig_bool[trig_num] && j_pt[j] >= trig_lower_0eta200[trig_num] && j_pt[j] < trig_upper_0eta200[trig_num]) { // if triggered, check whether the jet momentum falls in the correct range
		        			for (int k = 2; k < 6; k++) {
		        				if (jeta >= eta_cuts[k] && jeta < eta_cuts[k+1]) {
		        					harr[k]->Fill(j_pt[j], m_trig_prescale[trig_num]);
		        					break;
		        				}
		        			}
		        			break; // any jet should only be plotted once.
		        		}
		        	}
		        }
                        else if (jeta < 3.2 && jeta >= 2) {
		                for (int trig_num : trig_p200eta320) {
		        		if (m_trig_bool[trig_num] && j_pt[j] >= trig_lower_p200eta320[trig_num] && j_pt[j] < trig_upper_p200eta320[trig_num]) {
		        			harr[6]->Fill(j_pt[j], m_trig_prescale[trig_num]);
		        			break;
		        		}
		        	}
		        }
                        else if (jeta < 4.9 && jeta >= 3.2) {
		                for (int trig_num : trig_p320eta490) {
		        		if (m_trig_bool[trig_num] && j_pt[j] >= trig_lower_p320eta490[trig_num] && j_pt[j] < trig_upper_p320eta490[trig_num]) {
		        			harr[7]->Fill(j_pt[j], m_trig_prescale[trig_num]);
		        			break;
		        		}
		        	}
		        }
                        else if (jeta <= -2 && jeta > -3.2) {
		                for (int trig_num : trig_n200eta490) {
		        		if (m_trig_bool[trig_num] && j_pt[j] >= trig_lower_n200eta320[trig_num] && j_pt[j] < trig_upper_n200eta320[trig_num]) {
		        			harr[1]->Fill(j_pt[j], m_trig_prescale[trig_num]);
		        			break;
		        		}
		        	}
		        }
                        else if (jeta <= -3.2 && jeta > -4.9) {
		                for (int trig_num : trig_n200eta490) {
		        		if (m_trig_bool[trig_num] && j_pt[j] >= trig_lower_n320eta490[trig_num] && j_pt[j] < trig_upper_n320eta490[trig_num]) {
		        			harr[0]->Fill(j_pt[j], m_trig_prescale[trig_num]);
		        			break;
		        		}
		        	}
		        }
                	j++;
		}
	}

	// Write histograms to a root file
	TFile* output = new TFile(Form("./pt_data/run_%i.root", runNumber), "RECREATE");
	for (int i = 0; i< sizeof(harr)/sizeof(harr[0]); i++) {
		harr[i]->Scale((harr_scales[i]) / (A * luminosity * d_eta[i]), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
		harr[i]->Write();
	}
	output->Close();
}
