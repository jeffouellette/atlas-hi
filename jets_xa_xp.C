void jets_xa_xp(int runNumber, float luminosity) {

        luminosity = luminosity/1000; // convert from nb^(-1) to pb^(-1)

        float Z = 82; // value of Z for typical Pb
	float A = 208; // value of A for Pb
        float sqrt_s_nn = 8160;
	float d_eta = 0.5;
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
       

        // Map each trigger to a specific jet momentum cutoff range for each eta range 
        std::map<int, float> trig_lower_0eta200;
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
        std::map<int, float> trig_lower_200eta320;
        std::map<int, float> trig_upper_200eta320;
        trig_lower_200eta320[3] = 40;
        trig_upper_200eta320[3] = 50;
        trig_lower_200eta320[5] = 50;
        trig_upper_200eta320[5] = 55;
        trig_lower_200eta320[7] = 55;
        trig_upper_200eta320[7] = 60;
        trig_lower_200eta320[9] = 60;
        trig_upper_200eta320[9] = 70;
        trig_lower_200eta320[13] = 70;
        trig_upper_200eta320[13] = 75;
        trig_lower_200eta320[15] = 75;
        trig_upper_200eta320[15] = 85;
        trig_lower_200eta320[17] = 85;
        trig_upper_200eta320[17] = 110;
        trig_lower_200eta320[19] = 110;
        trig_upper_200eta320[19] = 2000;
        std::map<int, float> trig_lower_320eta490;
        std::map<int, float> trig_upper_320eta490;
        trig_lower_320eta490[1] = 25;
        trig_upper_320eta490[1] = 40;
        trig_lower_320eta490[3] = 40;
        trig_upper_320eta490[3] = 50;
        trig_lower_320eta490[5] = 50;
        trig_upper_320eta490[5] = 60;
        trig_lower_320eta490[9] = 60;
        trig_upper_320eta490[9] = 65;
        trig_lower_320eta490[11] = 65;
        trig_upper_320eta490[11] = 70;
        trig_lower_320eta490[13] = 70;
        trig_upper_320eta490[13] = 85;
        trig_lower_320eta490[17] = 85;
        trig_upper_320eta490[17] = 110;
        trig_lower_320eta490[19] = 110;
        trig_upper_320eta490[19] = 2000;

/*
	// Store file names to process as a TChain
	const char* files[10] = {
	//	"user.khill.jetTreeMaker.2.4.30hi.001.00313187.f774_m1736_myOutput_hadd.root",
		"user.khill.jetTreeMaker.2.4.30hi.001.00313572.f774_m1736_myOutput_hadd.root",
		"user.khill.jetTreeMaker.2.4.30hi.001.00313574.f774_m1736_myOutput_hadd.root",
		"user.khill.jetTreeMaker.2.4.30hi.001.00313575.f774_m1736_myOutput_hadd.root",
		"user.khill.jetTreeMaker.2.4.30hi.001.00313629.f781_m1741_myOutput_hadd.root",
		"user.khill.jetTreeMaker.2.4.30hi.001.00313630.f781_m1741_myOutput_hadd.root",
		"user.khill.jetTreeMaker.2.4.30hi.001.00313833.f781_m1741_myOutput_hadd.root",
		"user.khill.jetTreeMaker.2.4.30hi.001.00313929.f781_m1741_myOutput_hadd.root",
		"user.khill.jetTreeMaker.2.4.30hi.001.00314014.f781_m1741_myOutput_hadd.root",
		"user.khill.jetTreeMaker.2.4.30hi.001.00314157.f781_m1741_myOutput_hadd.root"
	};

	// Store luminosities of each run; must be of same length as files
	const float luminosities[sizeof(files)/sizeof(files[0])] = {3.24, 0.0051, 1.2, 7.059, 6.251, 6.632, 4.70, 0.0506, 6.850, 9.153}; //units in nb^-1
	float luminosity = 0;
	for (int i = 0; i<sizeof(files)/sizeof(files[0]); i++) {
		luminosity += luminosities[i]/1000; // divide by 1000 to convert nb^-1 to pb^-1
	}
*/

        TTree* tree = (TTree*)(new TFile(Form("./rundata/run_%i_raw.root", runNumber)))->Get("tree");

	// Create arrays to store trigger values for each event
	bool m_trig_bool[trigLength];
	float m_trig_prescale[trigLength];

	// Create a TChain which combines the tree of each file into one tree
	// Note that TChain inherits from TTree so it is a TTree as well
	TChain* tree = new TChain("tree");
	for (int i = 0; i < sizeof(files)/sizeof(files[0]); i++) {
		tree->Add(files[i]); // Add each file to the TChain
	}

	// Create arrays to store jet data for each event
	float j_pt[5] = {};
	float j_eta[5] = {};
        float j_phi[5] = {};
        float j_e[5] = {};
        int njet = 0;
	tree->SetBranchAddress("j_pt", j_pt);
	tree->SetBranchAddress("j_eta", j_eta);
        tree->SetBranchAddress("j_e", j_e);
        tree->SetBranchAddress("njet", &njet);
        tree->SetBranchAddress("j_phi", j_phi);

	const float jet_cuts[11] = {25, 40, 50, 55, 60, 65, 70, 75, 85, 110, 2000}; // cutoff momenta for each jet	
        
	const float eta_cuts[7] = {0, 0.5, 1, 1.5, 2, 2.5, 3};  // cutoff pseudorapidity for each bin
	
	TH1F* harr[2];

	harr[0] = new TH1F("xp", "#it{x}_{p};#it{x}_{i};d#sigma/d#it{x}_{i} [pb/GeV]", 80, 0, 1.6);
	harr[1] = new TH1F("xa", "#it{x}_{a}", 80, 0, 1.6);
	for (int i = 0; i< sizeof(harr)/sizeof(harr[0]); i++) {
		harr[i]->Sumw2();  // tell each histogram to propagate errors
	}
	
	// Set branch addresses
	for (int i = 0; i < trigLength; i++) {
		tree->SetBranchAddress(m_trig_string[i], &m_trig_bool[i]); 
		tree->SetBranchAddress(Form("%s_prescale", m_trig_string[i]), &m_trig_prescale[i]);
	}

	// Iterate over each event
	int numentries = tree->GetEntries();
	for (int i = 0; i < numentries; i++) {
		tree->GetEntry(i); // stores trigger values and data in the designated branch addresses
                if (njet == 2 || (njet == 3 && (j_pt[2]-0.5*(j_pt[0]+j_pt[1]))/(j_pt[2]+0.5*j_pt[0]+0.5*j_pt[1]) <= 0.1)) {	// select 2 jet events or 3 jet events with 2 dominating jets
		        for (int trig_num = 1; trig_num < 21; trig_num+=2) { // iterate over each trigger
			        if (m_trig_bool[trig_num] && j_pt[0] >= jet_cuts[(trig_num-1)/2] && j_pt[0] < jet_cuts[(trig_num+1)/2]) { // if triggered, check whether the jet momentum falls in the correct range
                                        harr[0]->Fill((TMath::Sqrt(Z/A) / sqrt_s_nn) * (j_pt[0]*TMath::Exp(j_eta[0]+0.465)+j_pt[1]*TMath::Exp(j_eta[1]+0.465)), m_trig_prescale[trig_num]);
                                        harr[1]->Fill((TMath::Sqrt(A/Z) / sqrt_s_nn) * (j_pt[0]*TMath::Exp(-j_eta[0]-0.465)+j_pt[1]*TMath::Exp(-j_eta[1]-0.465)), m_trig_prescale[trig_num]);
                                }
				break; // probably unnecessary, but any jet should only be plotted once.
			}
		}
	}

	// Save to root file
	TFile* output = new TFile(Form("./xdata/run%i.root", runNumber), "RECREATE");
	for (int i = 0; i< sizeof(harr)/sizeof(harr[0]); i++) {
		harr[i]->Scale(1/(A*luminosity), "width"); // each bin stores dN, so the cross section should be the histogram rescaled by the total luminosity, then divided by the pseudorapidity width
		harr[i]->Write();
	}
	output->Close();
}
