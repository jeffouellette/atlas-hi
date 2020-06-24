#ifndef __AtlasUtils_cxx__
#define __AtlasUtils_cxx__

#include "AtlasUtils.h"
#include "AtlasStyle.h"
#include "GlobalParams.h"

#include <Utilities.h>

#include <TROOT.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TPave.h>
#include <TH1.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TF1.h>

#include <iostream>
#include <cmath>


void SetStyle () {
  TStyle* myStyle = AtlasStyle();
  //myStyle->SetPalette(kRainBow);
  myStyle->SetPalette(kBird);
  myStyle->SetPadRightMargin(0.18);
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
  return;
}

TGraphErrors* myTGraphErrorsDivide (TGraphErrors* g1, TGraphErrors* g2) {
 
  const int debug=0; 

  if (!g1) printf("**myTGraphErrorsDivide: g1 does not exist !  \n"); 
  if (!g2) printf("**myTGraphErrorsDivide: g2 does not exist !  \n"); 


  int n1=g1->GetN();
  int n2=g2->GetN();

  if (n1!=n2) {
   printf("**myTGraphErrorsDivide: vector do not have same number of entries !  \n"); 
  }

  TGraphErrors* g3= new TGraphErrors();

  double  x1=0., y1=0., x2=0., y2=0.;
  double dx1=0.,dy1=0.,       dy2=0.;

  int iv=0;
  for (int i1=0; i1<n1; i1++) {
   for (int i2=0; i2<n2; i2++) {
     //if (debug) printf("**myTGraphErrorsDivide: %d  %d !  \n",i1,i2);

    g1->GetPoint(i1,x1,y1);
    g2->GetPoint(i2,x2,y2);
    if (x1!=x2) {
      //printf("**myTGraphErrorsDivide: %d x1!=x2  %f %f  !  \n",iv,x1,x2);
    }else{
      //if (debug) printf("**myTGraphErrorsDivide: %d x1=x2  %f %f  !  \n",iv,x1,x2);
     dx1  = g1->GetErrorX(i1);
     if (y1!=0) dy1  = g1->GetErrorY(i1)/y1;
     if (y2!=0) dy2  = g2->GetErrorY(i2)/y2;
   
     if (debug)
      printf("**myTGraphErrorsDivide: %d x1=%f x2=%f y1=%f y2=%f  \n",iv,x1,x2,y1,y2);

     if (y2!=0.) g3->SetPoint(iv, x1,y1/y2);
     else        g3->SetPoint(iv, x1,y2);
   
     double e=0.;
     if (y1!=0 && y2!=0) e=std::sqrt(dy1*dy1+dy2*dy2)*(y1/y2); 
     g3->SetPointError(iv,dx1,e);


     if (debug) {
       //double g3y, g3x,g3e;
       //g3->GetPoint(iv, g3y,g3x);
       //g3e=g3->GetErrorY(iv);
       //printf("%d g3y= %f g3e=%f  \n",iv,g3y,g3e);
     }
     iv++;
    }
    //    printf("**myTGraphErrorsDivide: ...next  \n");
   }
  }  
  return g3;

}




TGraphAsymmErrors* multiplyTGraphAsymmErrors (TGraphAsymmErrors* g1, TGraphAsymmErrors* g2) {
    TGraphAsymmErrors* g3 = new TGraphAsymmErrors();

    int n1 = g1->GetN();
    int n2 = g2->GetN();
    assert(n1 == n2);

    double x1, y1, x2, y2, x3, y3; // g1, g2, and g3 points
    double dx1h, dx1l, dy1h, dy1l; // g1 errors
    double dx2h, dx2l, dy2h, dy2l; // g2 errors
    double dx3h, dx3l, dy3h, dy3l; // g3 errors

    double* x1arr = g1->GetX();
    double* y1arr = g1->GetY();
    double* x2arr = g2->GetX();
    double* y2arr = g2->GetY();

    double* dx1harr = g1->GetEXhigh();
    double* dx1larr = g1->GetEXlow();
    double* dy1harr = g1->GetEYhigh();
    double* dy1larr = g1->GetEYlow();

    double* dx2harr = g2->GetEXhigh();
    double* dx2larr = g2->GetEXlow();
    double* dy2harr = g2->GetEYhigh();
    double* dy2larr = g2->GetEYlow();
    
    for (int bin = 0; bin < n1; bin++) {
        x1 = x1arr[bin];
        y1 = y1arr[bin];
        x2 = x2arr[bin];
        y2 = y2arr[bin];
        dx1h = dx1harr[bin];
        dx1l = dx1larr[bin];
        dy1h = dy1harr[bin];
        dy1l = dy1larr[bin];
        dx2h = dx2harr[bin];
        dx2l = dx2larr[bin];
        dy2h = dy2harr[bin];
        dy2l = dy2larr[bin];
        assert(x1 == x2);
        x3 = x1;
        y3 = y1*y2;
        dx3h = dx1h;
        dx3l = dx1l;
        dy3h = sqrt(dy1h*dy1l*y2*y2 + y1*y1*dy2h*dy2h);
        dy3l = sqrt(dy1l*dy1l*y2*y2 + y1*y1*dy2l*dy2l);

        g3->SetPoint(bin, x3, y3);
        g3->SetPointError(bin, dx3l, dx3h, dy3l, dy3h);
    }
    return g3;
}




TGraphAsymmErrors* myTGraphErrorsDivide (TGraphAsymmErrors* g1, TGraphAsymmErrors* g2) {

  const int debug=0; 

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();
  int n1=g1->GetN();
  int n2=g2->GetN();

  if (n1!=n2) {
    printf(" vectors do not have same number of entries !  \n");
   return g3;
  }

  double   x1=0.,   y1=0., x2=0., y2=0.;
  double dx1h=0., dx1l=0.;
  double dy1h=0., dy1l=0.;
  double dy2h=0., dy2l=0.;

  double* X1 = g1->GetX();
  double* Y1 = g1->GetY();
  double* EXhigh1 = g1->GetEXhigh();
  double* EXlow1 =  g1->GetEXlow();
  double* EYhigh1 = g1->GetEYhigh();
  double* EYlow1 =  g1->GetEYlow();

  double* X2 = g2->GetX();
  double* Y2 = g2->GetY();
  double* EXhigh2 = g2->GetEXhigh();
  double* EXlow2 =  g2->GetEXlow();
  double* EYhigh2 = g2->GetEYhigh();
  double* EYlow2 =  g2->GetEYlow();

  for (int i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i,x1,y1);
    g2->GetPoint(i,x2,y2);
    dx1h  = EXhigh1[i];
    dx1l  = EXlow1[i];
    if (y1!=0.) dy1h  = EYhigh1[i]/y1;
    else        dy1h  = 0.;
    if (y2!=0.) dy2h  = EYhigh2[i]/y2;
    else        dy2h  = 0.;
    if (y1!=0.) dy1l  = EYlow1 [i]/y1;
    else        dy1l  = 0.;
    if (y2!=0.) dy2l  = EYlow2 [i]/y2;
    else        dy2l  = 0.;
   
    //if (debug)
    //printf("%d x1=%f x2=%f y1=%f y2=%f  \n",i,x1,x2,y1,y2);
    if (debug)
      printf("%d dy1=%f %f dy2=%f %f sqrt= %f %f \n",i,dy1l,dy1h,dy2l,dy2h,
	     std::sqrt(dy1l*dy1l+dy2l*dy2l), std::sqrt(dy1h*dy1h+dy2h*dy2h));

    if (y2!=0.) g3->SetPoint(i, x1,y1/y2);
    else       g3->SetPoint(i, x1,y2);
    double el=0.; double eh=0.;

    if (y1!=0. && y2!=0.) el=std::sqrt(dy1l*dy1l+dy2l*dy2l)*(y1/y2);
    if (y1!=0. && y2!=0.) eh=std::sqrt(dy1h*dy1h+dy2h*dy2h)*(y1/y2);

    if (debug) printf("dx1h=%f  dx1l=%f  el=%f  eh=%f \n",dx1h,dx1l,el,eh);
    g3->SetPointError(i,dx1h,dx1l,el,eh);

  }  
  return g3;

}




TGraphAsymmErrors* myMakeBand (TGraphErrors* g0, TGraphErrors* g1, TGraphErrors* g2) {
  // default is g0
    //const int debug=0;

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();

  double  x1=0., y1=0., x2=0., y2=0., y0=0, x3=0.;
  //double dx1=0.;
  double dum;
  for (int i=0; i<g1->GetN(); i++) {
    g0->GetPoint(i, x1,y0);
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    // if (y1==0) y1=1;
    //if (y2==0) y2=1;

    if (i==g1->GetN()-1) x2=x1;
    else                 g2->GetPoint(i+1,x2,dum);

    if (i==0)            x3=x1;
    else                 g2->GetPoint(i-1,x3,dum);

    double tmp=y2;
    if (y1<y2) {y2=y1; y1=tmp;}
    //double y3=1.;
    double y3=y0;
    g3->SetPoint(i,x1,y3);

    double binwl=(x1-x3)/2.;
    double binwh=(x2-x1)/2.;
    if (binwl==0.)  binwl= binwh;
    if (binwh==0.)  binwh= binwl;
    g3->SetPointError(i,binwl,binwh,(y3-y2),(y1-y3));

  }
  return g3;

}




void myAddtoBand (TGraphErrors* g1, TGraphAsymmErrors* g2) {

  double  x1=0., y1=0.,  y2=0., y0=0;
  //double dx1=0.;
  //double dum;

  if (g1->GetN()!=g2->GetN())
    std::cout << " graphs have not the same # of elements " << std::endl;

  double* EYhigh = g2-> GetEYhigh();
  double* EYlow  = g2-> GetEYlow();

  for (int i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);
    
    if ( y1==0 || y2==0 ) { 
      std::cerr << "check these points very carefully : myAddtoBand() : point " << i << std::endl;  
    }
    //    if (y1==0) y1=1;
    //    if (y2==0) y2=1;

    //    if (i==g1->GetN()-1) x2=x1;
    //    else                 g2->GetPoint(i+1,x2,dum);
    //    if (i==0)            x3=x1;
    //    else                 g2->GetPoint(i-1,x3,dum);

    double eyh=0., eyl=0.;
    //if (y1<y2) {y2=y1; y1=tmp;}
    //double y3=1.;

    //printf("%d: y1=%f y2=%f Eyhigh= %f Eylow= %f \n",i,y1,y2,EYhigh[i],EYlow[i]);

    y0=y1-y2;
    if (y0!=0) {
     if (y0>0){
      eyh=EYhigh[i];
      eyh=std::sqrt(eyh*eyh+y0*y0);
      //printf("high: %d: y0=%f eyh=%f  \n",i,y0,eyh);
      g2->SetPointEYhigh(i,eyh);
     } else {
      eyl=EYlow[i];
      eyl=std::sqrt(eyl*eyl+y0*y0);
      // printf("low: %d: y0=%f eyl=%f  \n",i,y0,eyl);
      g2->SetPointEYlow (i,eyl);
     }
    }
  }
  return;

}



/**
 * Returns true iff this eta, phi coordinate lies in the disabled HEC region.
 */
bool InDisabledHEC (const double eta, double phi, const double dr) {
  phi = InTwoPi (phi);
  return 1.5-dr < eta && eta < 3.2+dr &&
         M_PI-dr < phi && phi < 3*M_PI/2+dr;
}


/**
 * Returns true iff this eta lies within the EMCal.
 */
bool InEMCal (const float eta) {
  return TMath::Abs(eta) < 1.37 || (1.56 < TMath::Abs(eta) && TMath::Abs(eta) < 2.47);
}


/**
 * Returns true iff this object is within a given radius in the HCal.
 */
bool InHadCal (const float eta, const float R) {
  return TMath::Abs(eta) <= 4.9 - R;
}


///**
// * Modifies the directory strings to point to the correct locations.
// */
//void SetupDirectories (const TString dataSubDir, const TString thisWorkPath) {
//  ResetDirectories ();
//
//  workPath = homePath + "/" + thisWorkPath;
//  externalWorkPath = drivePath + "/" + thisWorkPath;
//  intWorkPath = intPath + "/" + thisWorkPath;
//
//  rootPath = intWorkPath + "/rootFiles/" + dataSubDir;
//  dataPath = externalWorkPath + "/data/" + dataSubDir;
//  plotPath = workPath + "/Plots/" + dataSubDir;
//  ptPath = rootPath + "/ptData/";
//  trigPath = rootPath + "/trigData/";
//  effPath = rootPath + "/effData/";
//  xPath = rootPath + "/xData/";
//  RpPbPath = rootPath + "/RpPbData/";
//
//  return;
//}
//
//
///**
// * Clears sub-directory information from the directory strings
// */
//void ResetDirectories () {
//  workPath = ""; // Home analysis directory, should be modified in code outside this path structure
//  externalWorkPath = ""; // External drive storage directory, should be modified in code below
//  intWorkPath = ""; // Base directory for intermediate root files, should be modified in code outside this path structure
//  rootPath = ""; // Where analyzed *.root files are stored. Different analysis modules have different subdirectories here.
//  dataPath = ""; // Where the *.root raw data files (from the CERN grid) are stored.
//  plotPath = ""; // Where plots are stored.
//  ptPath = ""; // Where the pt analysis module output is stored.
//  trigPath = "";  // Where the trigger fire count module output is stored.
//  effPath = ""; // Where the trigger efficiency module output is stored.
//  xPath = ""; // Where the xa/xp module output is stored.
//  RpPbPath = ""; // Where the R_pPb module output is stored.
//}
//




///**
// * Returns the appropriate file in the given directory.
// * For MC, inFileName MUST be specified.
// */
//TFile* GetFile (const char* directory, const int dataSet, const bool isMC, const char* inFileName) {
//  TFile* file = NULL;
//
//  // First figure out the file we are looking for
//  TString fileIdentifier;
//  if (TString (inFileName) == "") {
//   if (!isMC) fileIdentifier = to_string (dataSet);
//   else {
//    cout << "Error: In Utilities.cxx: Cannot identify this MC file! Will return null!" << endl;
//    return NULL;
//   }
//  }
//  else fileIdentifier = inFileName;
//
//  // Now get the list of files
//  const TString dataPathTemp = dataPath + "/" + directory + "/";
//  TSystemDirectory dir (dataPathTemp.Data (), dataPathTemp.Data ());
//  TList* sysfiles = dir.GetListOfFiles ();
//  if (!sysfiles) {
//   cout << "Error: In Utilities.cxx: Cannot get list of files! Will return null!" << endl;
//   return NULL;
//  }
//  TSystemFile* sysfile;
//  TString fname;
//  TIter next (sysfiles);
//
//  while ( (sysfile = (TSystemFile*)next ())) {
//   fname = sysfile->GetName ();
//   if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
//    if (debugStatements) cout << "Status: In Utilities.cxx: Found " << fname.Data () << endl;
//    
//    if (fname.Contains (fileIdentifier)) {
//     file = new TFile (dataPathTemp+fname, "READ");
//     break;
//    }
//   }
//  }
//
//  if (!file) {
//   cout << "Error: In Utilities.cxx: TFile not obtained for given data set. Will return null!" << endl;
//   return NULL;
//  }
//  else return file;
//}
  

/**
 * Returns an abbreviated, unique identifier for a given dataset.
 */
TString GetIdentifier (const int dataSet, const char* inFileName, const bool isMC, const bool isSignalOnlySample, const bool periodA) {
  if (!isMC) return to_string (dataSet);

  TString id = (periodA ? "pPb_" : "Pbp_");

  id = id + (isSignalOnlySample ? "Signal_" : "Overlay_");

  if (TString (inFileName).Contains ("42310") && TString (inFileName).Contains ("Slice")) { // gamma+jet
   if (dataSet <= 0) return "";
   id = id + "GammaJet_Slice" + to_string (dataSet);
  }
  else if (TString (inFileName).Contains ("ZeeJet")) { // Zee+jet
   if (dataSet < 0) return "";
   id = id + "ZeeJet" + (dataSet == 0 ? "" : "_Slice" + to_string (dataSet));
  }
  else if (TString (inFileName).Contains ("ZmumuJet")) { // Zmumu+jet
   if (dataSet != 0) return "";
   id = id + "ZmumuJet";
  }
  else if (TString (inFileName).Contains ("Dstar")) { // D* samples
   id = id + "Dstar";
   if (TString (inFileName).Contains ("Minus"))
    id = id + "Minus";
   else if (TString (inFileName).Contains ("Plus"))
    id = id + "Plus";

   if (TString (inFileName).Contains ("JZ1WA"))
    id = id + "_JZ1WA";
   else if (TString (inFileName).Contains ("JZRW1B"))
    id = id + "_JZRW1B";
  }
  else if (TString (inFileName).Contains ("jetjet")) { // dijet
   if (dataSet <= 0) return "";
   id = id + "Dijet_Slice" + to_string (dataSet);
  }

  return id;
}




/**
 * Returns the TProfile of an input histogram along the x axis. Can use either statistical mean or gaussian mean.
 */
TH1D* GetProfileX (const TString name, TH2D* hist, const int nbinsx, const double* xbins, const bool useFit, const double* xlos, const double* xhis) {

  TH1D* prof = new TH1D (name, "", nbinsx, xbins);

  const bool useCustomBounds = (xlos && xhis);

  for (int xbin = 1; xbin <= nbinsx; xbin++) {
    TH1D* proj = hist->ProjectionY ("proj", xbin, xbin);
    if (proj->Integral () != 0)
      proj->Scale (1.0 / proj->Integral ());

    double mean, mean_err;

    // Calculate gaussian mean
    if (useFit) {
      TF1* gaus;
      if (useCustomBounds)
        gaus = new TF1 ("gaus", "gaus(0)", xlos[xbin-1], xhis[xbin-1]);
      else
        gaus = new TF1 ("gaus", "gaus(0)", proj->GetMean ()-2.8*proj->GetStdDev (), proj->GetMean ()+2.8*proj->GetStdDev ());
      gaus->SetParameter (1, proj->GetMean ());
      gaus->SetParameter (2, proj->GetStdDev ());
      proj->Fit (gaus, "Q0R");
      mean = gaus->GetParameter (1);
      mean_err = gaus->GetParError (1);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // Calculate statistical mean
    else {
      mean = proj->GetMean ();
      mean_err = proj->GetMeanError ();
    }

    prof->SetBinContent (xbin, mean);
    prof->SetBinError (xbin, mean_err);
    if (proj) { delete proj; proj = NULL; }
  }

  return prof;
}


/**
 * Returns the TProfile of an input histogram along the y axis. Can use either statistical mean or gaussian mean.
 */
TH1D* GetProfileY (const TString name, TH2D* hist, const int nbinsy, const double* ybins, const bool useFit, const double* ylos, const double* yhis) {

  TH1D* prof = new TH1D (name, "", nbinsy, ybins);

  const bool useCustomBounds = (ylos && yhis);

  for (int ybin = 1; ybin <= nbinsy; ybin++) {
    TH1D* proj = hist->ProjectionX ("proj", ybin, ybin);
    if (proj->Integral () != 0)
      proj->Scale (1.0 / proj->Integral ());

    double mean, mean_err;

    // Calculate gaussian mean
    if (useFit) {
      TF1* gaus;
      if (useCustomBounds)
        gaus = new TF1 ("gaus", "gaus(0)", ylos[ybin-1], yhis[ybin-1]);
      else
        gaus = new TF1 ("gaus", "gaus(0)", proj->GetMean ()-2.8*proj->GetStdDev (), proj->GetMean ()+2.8*proj->GetStdDev ());
      gaus->SetParameter (1, proj->GetMean ());
      gaus->SetParameter (2, proj->GetStdDev ());
      proj->Fit (gaus, "Q0R");
      mean = gaus->GetParameter (1);
      mean_err = gaus->GetParError (1);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // Calculate statistical mean
    else {
      mean = proj->GetMean ();
      mean_err = proj->GetMeanError ();
    }

    prof->SetBinContent (ybin, mean);
    prof->SetBinError (ybin, mean_err);
    if (proj) { delete proj; proj = NULL; }
  }

  return prof;
}


/**
 * Returns a histogram with the TProfile of data over the TProfile of MC along either the x or y axes. Can use either the statistical or gaussian mean.
 */
TH1D* GetDataOverMC (const TString name, TH2D* data, TH2D* mc, const int numbins, const double* bins, const bool useFit, const TString axis, const double* los, const double* his) {

  // figure out which axis to use
  if (axis != "x" && axis != "y" && axis != "X" && axis != "Y") {
    cout << "JetCalibration::GetDataOverMC: Invalid axis specified!" << endl;
    return NULL;
  }
  const bool useXaxis = (axis == "x" || axis == "X");
  const bool useCustomBounds = (los && his);

  TH1D* dataOverMC = new TH1D (name, "", numbins, bins);

  TH1D* proj = NULL;

  double dataAvg, dataErr, mcAvg, mcErr;

  for (int bin = 1; bin <= numbins; bin++) {

    // first calculate the data value (numerator)
    if (useXaxis) proj = data->ProjectionY (name + Form ("data_bin%i", bin), bin, bin);
    else proj = data->ProjectionX (name + Form ("data_bin%i", bin), bin, bin);
    if (proj->Integral () != 0)
      proj->Scale (1.0 / proj->Integral ());

    // Calculate gaussian mean
    if (useFit) {
      TF1* gaus;
      if (useCustomBounds)
        gaus = new TF1 ("gaus", "gaus(0)", los[bin-1], his[bin-1]);
      else
        gaus = new TF1 ("gaus", "gaus(0)", proj->GetMean ()-2.8*proj->GetStdDev (), proj->GetMean ()+2.8*proj->GetStdDev ());
      gaus->SetParameter (1, proj->GetMean ());
      gaus->SetParameter (2, proj->GetStdDev ());
      proj->Fit (gaus, "Q0R");
      dataAvg = gaus->GetParameter (1);
      dataErr = gaus->GetParError (1);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // Calculate statistical mean
    if (!useFit) {
      dataAvg = proj->GetMean ();
      dataErr = proj->GetMeanError ();
    }

    if (proj) { delete proj; proj = NULL; }

    // next calculate the MC value (denominator)
    if (useXaxis) proj = mc->ProjectionY (name + Form ("mc_bin%i", bin), bin, bin);
    else proj = mc->ProjectionX (name + Form ("mc_bin%i", bin), bin, bin);
    if (proj->Integral () != 0)
      proj->Scale (1.0 / proj->Integral ());

    // Calculate gaussian mean
    if (useFit) {
      TF1* gaus;
      if (useCustomBounds)
        gaus = new TF1 ("gaus", "gaus(0)", los[bin-1], his[bin-1]);
      else
        gaus = new TF1 ("gaus", "gaus(0)", proj->GetMean ()-2.8*proj->GetStdDev (), proj->GetMean ()+2.8*proj->GetStdDev ());
      gaus->SetParameter (1, proj->GetMean ());
      gaus->SetParameter (2, proj->GetStdDev ());
      proj->Fit (gaus, "Q0R");
      mcAvg = gaus->GetParameter (1);
      mcErr = gaus->GetParError (1);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // Calculate statistical mean
    else {
      mcAvg = proj->GetMean ();
      mcErr = proj->GetMeanError ();
    }

    if (proj) { delete proj; proj = NULL; }

    // now set the nominal value to the numerator over denominator
    if (!(dataAvg == 0 || isnan (dataAvg) || mcAvg == 0 || isnan (mcAvg))) {
      const double dataOverMCavg = dataAvg/mcAvg;
      const double dataOverMCerr = sqrt (pow (dataErr/mcAvg, 2) + pow (dataOverMCavg * mcErr/mcAvg, 2));
      dataOverMC->SetBinContent (bin, dataOverMCavg);
      dataOverMC->SetBinError (bin, dataOverMCerr);
    }
  }

  if (useXaxis) dataOverMC->GetXaxis ()->SetTitle (data->GetXaxis ()->GetTitle ()); // copy x axis title
  else dataOverMC->GetXaxis ()->SetTitle (data->GetYaxis ()->GetTitle ()); // copy y axis title
  return dataOverMC;
}


/**
 * Converts a TProfile to a TH1D.
 */
TH1D* TProfile2TH1D (const char* name, TProfile* p, const int nx, const double* x) {
  TH1D* h = new TH1D (name, "", nx, x);
  h->Sumw2 ();

  for (int ix = 1; ix <= nx; ix++) {
    h->SetBinContent (ix, p->GetBinContent (ix));
    h->SetBinError (ix, p->GetBinError (ix));
  }
  return h;
}


/**
 * Reflects the contents of h around the n-th bin in x.
 */
void GetReflectionX (TH1D* h, const int n) {
  if (n < 1 || h->GetNbinsX () < n) {
    cout << "Invalid choice of n! Exiting gracefully." << endl;
    return;
  }

  const int d = std::max (n, h->GetNbinsX () - n);
  int start = n - d + (h->GetNbinsX () % 2) + 1;
  int end = n + d;
  while (start < end) {
    if (1 <= start && end <= h->GetNbinsX ()) {
      const double temp = h->GetBinContent (start); 
      const double temperr = h->GetBinError (start);

      h->SetBinContent (start, h->GetBinContent (end)); 
      h->SetBinError (start, h->GetBinError (end));

      h->SetBinContent (end, temp);
      h->SetBinError (end, temperr);
    }
    else {
      if (start < 1) {
        h->SetBinContent (end, 0);
        h->SetBinError (end, 0);
      }
      if (end > h->GetNbinsX ()) {
        h->SetBinContent (start, 0);
        h->SetBinError (start, 0);
      }
    }
    start++; 
    end--; 
  } 

  return;
}


#endif // __AtlasUtils_cxx__
