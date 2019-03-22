#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <sstream>

#include "Pythia8/Pythia.h"
#include "Pythia8/FJcore.h"

#include <Utilities.h>
#include <GlobalParams.h>

using namespace atlashi;
using namespace Pythia8;

int main (int argc, char *argv[]) {

  if (argc != 4) {

    std::cout << " usage: ./zgen MINPT NEVT FILENAMEOUT" << std::endl;

    return 0;
  }

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;

  pythia.readString ("Beams:eCM = 5020.");
  pythia.readString("PromptPhoton:all  = on");

  ostringstream ss; ss << "PhaseSpace:pTHatMin = " << argv[1] << ".";
  //pythia.readString("PhaseSpace:pTHatMin = 10.");
  pythia.readString (ss.str ().c_str ());

  pythia.init ();

  SlowJet *antikT4 = new SlowJet (-1, 1.0, 1, 5, 2, 1);

  int NEVT = atoi (argv[2]);

  TFile *f = new TFile (argv[3] , "RECREATE");

  int b_code;
  int b_id1;
  int b_id2;
  float b_x1pdf;
  float b_x2pdf;
  float b_Q;
  bool b_isValence1;
  bool b_isValence2;

  int b_photon_n, b_jet_n;
  vector<float> b_photon_pt, b_photon_eta, b_photon_phi, b_photon_m;
  vector<float> b_jet_pt, b_jet_eta, b_jet_phi, b_jet_e;

  TTree *t = new TTree("tree","a shambling vine tree");

  t->Branch ("code",  &b_code);
  t->Branch ("id1",   &b_id1);
  t->Branch ("id2",   &b_id2);
  t->Branch ("x1pdf", &b_x1pdf);
  t->Branch ("x2pdf", &b_x2pdf);
  t->Branch ("Q",     &b_Q);

  t->Branch ("photon_n",   &b_photon_n);
  t->Branch ("photon_pt",  &b_photon_pt);
  t->Branch ("photon_eta", &b_photon_eta);
  t->Branch ("photon_phi", &b_photon_phi);

  t->Branch ("jet_n",   &b_jet_n);
  t->Branch ("jet_pt",  &b_jet_pt);
  t->Branch ("jet_eta", &b_jet_eta);
  t->Branch ("jet_phi", &b_jet_phi);
  t->Branch ("jet_e",   &b_jet_e);

  TLorentzVector l1, l2;
  
  for (int iEvent = 0; iEvent < NEVT; iEvent++) {
    if (!pythia.next ())
      continue;

    b_photon_n = 0;
    b_photon_pt.clear ();
    b_photon_eta.clear ();
    b_photon_phi.clear ();

    for (int i = 0; i < pythia.event.size (); i++) {

      if (abs (pythia.event[i].id ()) != 22) continue; // check if photon

      b_photon_pt.push_back (pythia.event[i].pT ());
      b_photon_eta.push_back (pythia.event[i].eta ());
      b_photon_phi.push_back (pythia.event[i].phi ());
      b_photon_n++;
    }

    if (b_photon_n == 0) {
      iEvent--;
      continue;
    }

    antikT4->analyze (pythia.event);

    b_code = pythia.info.code ();
    b_id1 = pythia.info.id1pdf ();
    b_id2 = pythia.info.id2pdf ();
    b_x1pdf = pythia.info.x1pdf ();
    b_x2pdf = pythia.info.x2pdf ();
    b_Q =  pythia.info.QFac ();
    
    b_isValence1 = pythia.info.isValence1 ();
    b_isValence2 = pythia.info.isValence2 ();

    b_jet_n = 0;
    b_jet_pt.clear ();
    b_jet_eta.clear ();
    b_jet_phi.clear ();
    b_jet_e.clear ();

    for (int i = 0; i < antikT4->sizeJet (); i++) {

      bool matchesPhoton = false;
      for (int j = 0; !matchesPhoton && j < b_photon_n; j++) {
        matchesPhoton = DeltaR (antikT4->p (i).eta (), b_photon_eta.at (j), antikT4->phi (i), b_photon_phi.at (j)) < 0.2;
      }
      if (matchesPhoton) continue;

      b_jet_pt.push_back (antikT4->pT (i));
      b_jet_eta.push_back (antikT4->p (i).eta ());
      b_jet_phi.push_back (antikT4->phi (i));
      b_jet_e.push_back (antikT4->p (i).e ());
      b_jet_n++;
    }

    t->Fill();

    if (iEvent % (NEVT/100) == 0)
      std::cout << iEvent / (NEVT/100) << "\% done...\r" << std::flush;
  }

  pythia.stat();
  
  f->Write();
  f->Close();

  return 0;
}
