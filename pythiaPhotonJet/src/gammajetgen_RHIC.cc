#include <TFile.h>
#include <TTree.h>

#include <sstream>

#include "Pythia8/Pythia.h"
using namespace Pythia8;

float deltaR( float eta1, float eta2, float phi1, float phi2 ) {

  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  if (dphi > 3.14159) dphi -= 2*3.14159;
  if (dphi < -3.14159) dphi += 2*3.14159;

  return sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );

}

int main(int argc, char *argv[]) {

  /*
  if (argc != 4) {

    std::cout << " usage: ./gammajetgen MINPT NEVT FILENAMEOUT" << std::endl;

    return 0;
  }
  */

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;

  pythia.readString("Beams:eCM = 200.");
  //pythia.readString("PromptPhoton:gg2gammagamma = on");
  pythia.readString("PromptPhoton:all  = on");

  ostringstream ss; ss << "PhaseSpace:pTHatMin = 30.";
  //pythia.readString("PhaseSpace:pTHatMin = 10.");
  pythia.readString( ss.str().c_str() );

  pythia.init();

  int NEVT = 1000;

  //TFile *f = new TFile( argv[3] , "RECREATE");

  int b_id1;
  int b_id2;
  float b_x1pdf;
  float b_x2pdf;
  float b_Q;
  bool b_isValence1;
  bool b_isValence2;

  int b_photon_n;
  float b_photon_pt[20];
  float b_photon_eta[20];
  float b_photon_iso[20];

  TTree *t = new TTree("tree","a shambling vine tree");

  t->Branch("id1",&b_id1);
  t->Branch("id2",&b_id2);
  t->Branch("x1pdf",&b_x1pdf);
  t->Branch("x2pdf",&b_x2pdf);
  t->Branch("Q",&b_Q);
  t->Branch("isValence1",&b_isValence1);
  t->Branch("isValence2",&b_isValence2);

  t->Branch("photon_n",&b_photon_n);
  t->Branch("photon_pt",b_photon_pt,"photon_pt[photon_n]/F");
  t->Branch("photon_eta",b_photon_eta,"photon_eta[photon_n]/F");
  t->Branch("photon_iso",b_photon_iso,"photon_iso[photon_n]/F");
  
  for (int iEvent = 0; iEvent < NEVT; ++iEvent) {
    if (!pythia.next()) continue;

    if ( iEvent % 1000 == 0 ) std::cout << iEvent << std::endl;

    //std::cout << " --> " << pythia.info.id1pdf()  << ", " <<  pythia.info.id2pdf() << " with x1/x2 = " <<  pythia.info.x1pdf() << " / " <<  pythia.info.x2pdf() << ", sqrt(Q2) = " << pythia.info.QFac() << std::endl;

    b_id1 = pythia.info.id1pdf();
    b_id2 = pythia.info.id2pdf();
    b_x1pdf = pythia.info.x1pdf();
    b_x2pdf = pythia.info.x2pdf();
    b_Q =  pythia.info.QFac();
    
    b_isValence1 = pythia.info.isValence1();
    b_isValence2 = pythia.info.isValence2();

    b_photon_n = 0;

    for (int i = 0; i < pythia.event.size(); ++i) {
      if (!pythia.event[i].isFinal()) continue;
      
      if (pythia.event[i].pT() < 10 || pythia.event[i].id() != 22 ) continue;
      
      b_photon_pt[ b_photon_n ] = pythia.event[i].pT();
      b_photon_eta[ b_photon_n ] = pythia.event[i].eta();
      
      //std::cout << "i = " << i << ", id = " << pythia.event[i].id() << ", pT/eta/phi = " <<  pythia.event[i].pT() << " / " <<  pythia.event[i].eta() << " / " <<  pythia.event[i].phi() << std::endl;
      
      b_photon_iso[ b_photon_n ] = 0;
      
      for (int i2 = 0; i2 < pythia.event.size(); ++i2) {
	
	if (!pythia.event[i2].isFinal()) continue;
	if (i == i2) continue;
	
	if (deltaR( pythia.event[i2].eta(), pythia.event[i].eta(), pythia.event[i2].phi(), pythia.event[i].phi() ) < 0.4)  b_photon_iso[ b_photon_n ] += pythia.event[i2].pT();
	
      }

      b_photon_n++;
      
    }
    
    t->Fill();

  }
  pythia.stat();

  
  //t->Write();
  //f->Write();
  //f->Close();

  return 0;
}
