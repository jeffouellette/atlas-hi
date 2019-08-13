#include "Params.h"

#include <ArrayTemplates.h>

#include <AtlasUtils.h>

void QuickZMassPlot () {

  TFile* dataFile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/savedHists.root", "read");
  TH1D**** h_z_m = Get3DArray <TH1D*> (numCentBins, 2, 3);          // iCent, iSpc, iReg
  for (int iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 2; iSpc++) {
      const char* spc = (iSpc == 0 ? "ee" : "mumu");
      for (short iReg = 0; iReg < 3; iReg++) {
        h_z_m[iCent][iSpc][iReg] = (TH1D*) dataFile->Get (Form ("h_z_m_%s_iCent%i_iReg%i_data", spc, iCent, iReg));
        h_z_m[iCent][iSpc][iReg]->Scale (1./ h_z_m[iCent][iSpc][iReg]->Integral ());
      }
    }
  }

  for (short iSpc = 0; iSpc < 2; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : "mumu");
    for (short iReg = 0; iReg < 2; iReg++) {
      for (short iCent = numCentBins-1; iCent >= 1; iCent--) {
        const char* canvasName = Form ("c_z_m_%s_iCent%i_iReg%i", spc, iCent, iReg);
        TCanvas* c = nullptr;
        TPad* uPad = nullptr, *dPad = nullptr;
        c = new TCanvas (canvasName, "", 800, 800);

        c->cd ();
        uPad = new TPad (Form ("%s_uPad", canvasName), "", 0.0, 0.4, 1.0, 1.0);
        dPad = new TPad (Form ("%s_dPad", canvasName), "", 0.0, 0.0, 1.0, 0.4);
        uPad->SetBottomMargin (0);
        dPad->SetTopMargin (0);
        dPad->SetBottomMargin (0.25);
        uPad->Draw ();
        dPad->Draw ();
        gDirectory->Add (c);
        gDirectory->Add (uPad);
        gDirectory->Add (dPad);


        uPad->cd ();

        TH1D* h = h_z_m[0][iSpc][iReg];

        h->SetFillColorAlpha (kAzure+10, fillAlpha);
        h->SetLineColor (kBlack);
        h->SetMarkerSize (0);
        h->SetLineWidth (0);
        //h->GetYaxis ()->SetRangeUser (0, 1.3);
        h->GetYaxis ()->SetRangeUser (0, 0.12);

        h->GetXaxis ()->SetTitle (Form ("m_{%s} [GeV]", (iSpc == 0 ? "ee" : "#mu#mu")));
        //h->GetYaxis ()->SetTitle ("Arb. Units");
        h->GetYaxis ()->SetTitle ("Counts / Total");
        h->GetXaxis ()->SetTitleSize (0.04/0.6);
        h->GetYaxis ()->SetTitleSize (0.04/0.6);
        h->GetXaxis ()->SetLabelSize (0.04/0.6);
        h->GetYaxis ()->SetLabelSize (0.04/0.6);
        h->GetXaxis ()->SetTitleOffset (1.5*0.6);
        h->GetYaxis ()->SetTitleOffset (1.5*0.6);

        h->DrawCopy ("bar");
        h->SetLineWidth (1);
        h->Draw ("hist same");

        gPad->RedrawAxis ();

        TGraphAsymmErrors* g = make_graph (h_z_m[iCent][iSpc][iReg]);
        ResetXErrors (g);
        //deltaize (g, 0.1*(-1.5+iCent));

        const int markerStyle = kFullCircle;
        g->SetMarkerStyle (markerStyle);
        g->SetMarkerSize (1);
        g->SetLineWidth (1);
        g->SetLineColor (kBlack);
        g->SetMarkerColor (kBlack);
        //g->GetYaxis ()->SetRangeUser (0, 1.3);
        g->GetYaxis ()->SetRangeUser (0, 0.12);

        g->GetXaxis ()->SetTitle (Form ("m_{%s} [GeV]", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))));
        //g->GetYaxis ()->SetTitle ("Arb. Units");
        g->GetYaxis ()->SetTitle ("Counts / Total");
        g->GetXaxis ()->SetTitleSize (0.04/0.6);
        g->GetYaxis ()->SetTitleSize (0.04/0.6);
        g->GetXaxis ()->SetLabelSize (0.04/0.6);
        g->GetYaxis ()->SetLabelSize (0.04/0.6);
        g->GetXaxis ()->SetTitleOffset (1.5*0.6);
        g->GetYaxis ()->SetTitleOffset (1.5*0.6);
        g->Draw ("P");

        myText (0.22, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/0.6);
        const char* spc = iSpc == 0 ? "Z #rightarrow e^{+}e^{-}" : (iSpc == 1 ? "Z #rightarrow #mu^{+}#mu^{-}" : "Z #rightarrow l^{+}l^{-}");
        myText (0.66, 0.85, kBlack, spc, 0.04/0.6);
      
        myOnlyBoxText (0.71, 0.67, 1.2, kAzure+10, kBlack, 1, "#it{pp}", 0.04/0.6, 1001);
      
        //TVirtualPad* cPad = gPad; // store current pad
        //TBox* b = TBoxNDC (0.4+0.6*(0.598-0.025), 0.67-0.06*numPhiBins-0.018, 0.4+0.6*(0.598+0.025), 0.67-0.06*numPhiBins+0.018);
        //b->SetFillColorAlpha (fillColors[iCent], fillAlpha);
        //b->Draw ("l");
        //cPad->cd ();
        //myText (0.753, 0.67, kBlack, "MC", 0.04/0.6);
        myMarkerText (0.703, 0.76, kBlack, kFullCircle, Form ("Pb+Pb %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 1.25, 0.04/0.6);
      
        if (iReg == 0)
          myText (0.22, 0.67, kBlack, Form ("#left|y^{%s}#right| < 1", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))), 0.04/0.6);
        else if (iReg == 1)
          myText (0.22, 0.67, kBlack, Form ("#left|y^{%s}#right| > 1", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))), 0.04/0.6);


        dPad->cd ();

        h = (TH1D*) h_z_m[iCent][iSpc][iReg]->Clone (TString (h_z_m[iCent][iSpc][iReg]->GetName ()) + "_ratio");
        h->Divide (h_z_m[0][iSpc][iReg]);
        if (h) {
          TGraphAsymmErrors* g = make_graph (h);
          ResetXErrors (g);

          const int markerStyle = kFullCircle;
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerStyle (markerStyle);
          g->SetMarkerSize (1);
          g->SetLineWidth (1);
          g->SetLineColor (kBlack);
          g->SetMarkerColor (kBlack);
          g->GetYaxis ()->SetRangeUser (0.1, 2.4);

          g->GetXaxis ()->SetTitle (Form ("m_{%s} [GeV]", (iSpc == 0 ? "ee" : (iSpc == 1 ? "#mu#mu" : "ll"))));
          g->GetYaxis ()->SetTitle ("Pb+Pb / #it{pp}");
          g->GetXaxis ()->SetTitleSize (0.04/0.4);
          g->GetYaxis ()->SetTitleSize (0.04/0.4);
          g->GetXaxis ()->SetLabelSize (0.04/0.4);
          g->GetYaxis ()->SetLabelSize (0.04/0.4);
          g->GetXaxis ()->SetTitleOffset (2.5*0.4);
          g->GetYaxis ()->SetTitleOffset (1.5*0.4);
          g->GetYaxis ()->CenterTitle ();
          g->Draw ("AP");

          TLine* l = new TLine (76, 1, 106, 1);
          l->SetLineColor (46);
          l->SetLineWidth (2);
          l->SetLineStyle (5);
          l->Draw ("same");
          
        }
        else {
          cout << "Warning in FullAnalysis :: PlotZMassSpectra: Z mass spectra ratio not stored, needs to be calculated!" << endl;
        }

        if (iReg == 2)
          c->SaveAs (Form ("Plots/z%s_mass_spectrum_iCent%i.pdf", iSpc == 0 ? "ee" : "mumu", iCent));
        else
          c->SaveAs (Form ("Plots/z%s_mass_spectrum_iCent%i_iReg%i.pdf", iSpc == 0 ? "ee" : "mumu", iCent, iReg));
      }
    }
  }

}
