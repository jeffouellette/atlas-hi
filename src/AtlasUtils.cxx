
#include <iostream>
#include <cmath>

#include "../include/AtlasUtils.h"
#include "../include/AtlasStyle.h"

#include "TROOT.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TPave.h"
#include "TH1.h"
#include "TH1D.h"
#include "TStyle.h"

bool IsOpenMarker (const Style_t ms) {
  return ms == kOpenCircle ||
         ms == kOpenSquare ||
         ms == kOpenTriangleUp ||
         ms == kOpenDiamond ||
         ms == kOpenCross ||
         ms == kOpenStar ||
         ms == kOpenTriangleDown ||
         ms == kOpenDiamondCross ||
         ms == kOpenSquareDiagonal ||
         ms == kOpenThreeTriangles ||
         ms == kOpenFourTrianglesX ||
         ms == kOpenDoubleDiamond ||
         ms == kOpenFourTrianglesPlus ||
         ms == kOpenCrossX;
}

bool IsFullMarker (const Style_t ms) {
  return ms == kFullCircle ||
         ms == kFullSquare ||
         ms == kFullTriangleUp ||
         ms == kFullDiamond ||
         ms == kFullCross ||
         ms == kFullStar ||
         ms == kFullTriangleDown ||
         ms == kFullThreeTriangles ||
         ms == kFullFourTrianglesX ||
         ms == kFullDoubleDiamond ||
         ms == kFullFourTrianglesPlus ||
         ms == kFullCrossX;
}

Style_t FullToOpenMarker (const Style_t ms) {
  switch (ms) {
    case kFullCircle: return kOpenCircle;
    case kFullSquare: return kOpenSquare;
    case kFullDiamond: return kOpenDiamond;
    case kFullCross: return kOpenCross;
    case kFullTriangleUp: return kOpenTriangleUp;
    case kFullTriangleDown: return kOpenTriangleDown;
    case kFullStar: return kOpenStar;
    case kFullCrossX: return kOpenCrossX;
    case kFullFourTrianglesPlus: return kOpenFourTrianglesPlus;
    case kFullFourTrianglesX: return kOpenFourTrianglesX;
    case kFullThreeTriangles: return kOpenThreeTriangles;
    case kFullDoubleDiamond: return kOpenDoubleDiamond;
    default: return kDot;
  }
}


Style_t OpenToFullMarker (const Style_t ms) {
  switch (ms) {
    case kOpenCircle: return kFullCircle;
    case kOpenSquare: return kFullSquare;
    case kOpenDiamond: return kFullDiamond;
    case kOpenCross: return kFullCross;
    case kOpenTriangleUp: return kFullTriangleUp;
    case kOpenTriangleDown: return kFullTriangleDown;
    case kOpenStar: return kFullStar;
    case kOpenCrossX: return kFullCrossX;
    case kOpenFourTrianglesPlus: return kFullFourTrianglesPlus;
    case kOpenFourTrianglesX: return kFullFourTrianglesX;
    case kOpenThreeTriangles: return kFullThreeTriangles;
    case kOpenDoubleDiamond: return kFullDoubleDiamond;
    default: return kDot;
  }
}

void SetStyle () {
  TStyle* myStyle = AtlasStyle();
  //myStyle->SetPalette(kRainBow);
  myStyle->SetPalette(kBird);
  myStyle->SetPadRightMargin(0.18);
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
  return;
}

void FormatTH2Canvas (TCanvas* c, const bool zAxisSpace) {
  SetStyle();
  if (zAxisSpace) c->SetRightMargin (0.18);
  else c->SetRightMargin (0.06);

  c->SetLeftMargin (0.15); 
}


void ATLAS_LABEL (double x, double y,Color_t color) {
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x,y,"ATLAS");
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




TGraphErrors* TH1TOTGraph (TH1 *h1) {


  if (!h1) std::cout << "TH1TOTGraph: histogram not found !" << std::endl;

 TGraphErrors* g1= new TGraphErrors();

 double x, y, ex, ey;
 for (int i=1 ; i<=h1->GetNbinsX(); i++) {
   y=h1->GetBinContent(i);
   ey=h1->GetBinError(i);
   x=h1->GetBinCenter(i);
   ex=h1->GetBinWidth(i);
   
  //   cout << " x,y = " << x << " " << y << " ex,ey = " << ex << " " << ey << endl;

   g1->SetPoint(i-1,x,y);
   g1->SetPointError(i-1,ex,ey);

 }

 //g1->Print();

 return g1;
}




void BinomialDivide (TH1D* n, TH1D* d) {
  for (int ix = 1; ix <= n->GetNbinsX (); ix++) {
    float eff = n->GetBinContent (ix);
    if (d->GetBinContent (ix) != 0)
      eff = eff / d->GetBinContent (ix);

    float var = (eff - eff*eff);
    if (d->GetBinContent (ix) != 0)
      var = var / d->GetBinContent (ix);

    n->SetBinContent (ix, eff);
    n->SetBinError (ix, sqrt (fabs (var)));
  }
  return;
}




void myText (double x, double y, Color_t color, const char *text, double tsize) {

//  double tsize=0.04;
  TLatex l; /*l.SetTextAlign(12);*/ l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
 



void myBoxText (double x, double y, double boxsize, int mcolor, const char *text) {

  double tsize=0.06;

  TLatex l; l.SetTextAlign(12); //l.SetTextSize(tsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text);

  double y1=y-0.25*tsize;
  double y2=y+0.25*tsize;
  double x2=x-0.3*tsize;
  double x1=x2-boxsize;

  printf("x1= %f x2= %f y1= %f y2= %f \n",x1,x2,y1,y2);

  TPave *mbox= new TPave(x1,y1,x2,y2,0,"NDC");

  mbox->SetFillColor(mcolor);
  mbox->SetFillStyle(1001);
  mbox->Draw();

  TLine mline;
  mline.SetLineWidth(4);
  mline.SetLineColor(1);
  mline.SetLineStyle(1);
  double y_new=(y1+y2)/2.;
  mline.DrawLineNDC(x1,y_new,x2,y_new);

}




/*void myMarkerText(double x,double y,int color,int mstyle, const char *text) 
{
  double tsize=0.032;
  float msize = 0.75;
  TMarker *marker = new TMarker(x-(0.4*tsize),y,8);
  marker->SetMarkerColor(color);  marker->SetNDC();
  marker->SetMarkerStyle(mstyle);
  marker->SetMarkerSize(msize);
  marker->Draw();

  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text);
}*/




void myLineText (double x, double y, int color, int lstyle, const char *text, float lsize, double tsize) {
  TLine *markerLine = new TLine(x-(0.8*tsize)-0.02*lsize, y, x-(0.8*tsize)+0.02*lsize, y);
  markerLine->SetNDC();
  markerLine->SetLineColor(color);
  markerLine->SetLineStyle(lstyle);
  markerLine->SetLineWidth(2);
  markerLine->Draw();

  if (text[0] != '\0') {
    TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); /*l.SetTextColor (color);*/
    l.SetNDC();
    l.DrawLatex(x,y,text);
  }
}




void myLineColorText (double x, double y, int color, int lstyle, const char *text, float lsize, double tsize) {
  TLine *markerLine = new TLine(x-(0.8*tsize)-0.02*lsize, y, x-(0.8*tsize)+0.02*lsize, y);
  markerLine->SetNDC();
  markerLine->SetLineColor(color);
  markerLine->SetLineStyle(lstyle);
  markerLine->SetLineWidth(2);
  markerLine->Draw();

  if (text[0] != '\0') {
    TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); l.SetTextColor (color);
    l.SetNDC();
    l.DrawLatex(x,y,text);
  }
}




void myMarkerText (double x, double y, int color, int mstyle, const char *text, float msize, double tsize) {
//  double tsize=0.032;
  //TMarker *marker = new TMarker(x-(0.44*tsize),y,8);
  TMarker *marker = new TMarker(x-(0.8*tsize),y,8);
  marker->SetMarkerColor(color);  marker->SetNDC();
  marker->SetMarkerStyle(mstyle);
  marker->SetMarkerSize(msize);
  marker->Draw();

  //TLine *markerLine = new TLine(marker->GetX()-0.018*msize, marker->GetY(), marker->GetX()+0.018*msize, marker->GetY());
  TLine *markerLine = new TLine(x-(0.8*tsize)-0.02, y, x-(0.8*tsize)+0.02, y);
  markerLine->SetNDC();
  markerLine->SetLineColor(color);
  markerLine->SetLineStyle(1);
  markerLine->SetLineWidth(2);
  markerLine->Draw();

  if (text[0] != '\0') {
    TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); 
    l.SetNDC();
    l.DrawLatex(x,y,text);
  }
}




void myMarkerTextNoLine (double x, double y, int color, int mstyle, const char *text, float msize, double tsize) {
//  double tsize=0.032;
  //TMarker *marker = new TMarker(x-(0.44*tsize),y,8);
  TMarker *marker = new TMarker(x-(0.8*tsize),y,8);
  marker->SetMarkerColor(color);  marker->SetNDC();
  marker->SetMarkerStyle(mstyle);
  marker->SetMarkerSize(msize);
  marker->Draw();

  if (text[0] != '\0') {
    TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); 
    l.SetNDC();
    l.DrawLatex(x,y,text);
  }
}




void myOnlyBoxText (double x, double y, double boxsize, int mcolor, int lcolor, int lstyle, const char *text, double tsize, int bstyle, double balpha) {
  double y1=y-0.25*tsize;
  double y2=y+0.25*tsize;
  double x2=x-0.8*tsize+0.02*boxsize;
  double x1=x-0.8*tsize-0.02*boxsize;
  //double x2=x-0.15*tsize;
  //double x1=x-0.95*tsize;
  //printf("x1= %f x2= %f y1= %f y2= %f \n",x1,x2,y1,y2);
  TPave *mbox= new TPave(x1,y1,x2,y2,0,"NDC");
  mbox->SetFillColorAlpha(mcolor,balpha);
  mbox->SetFillStyle(bstyle);
  mbox->Draw();
  TLine mline;
  mline.SetLineWidth(1);
  mline.SetLineColor(lcolor);
  //mline.SetLineStyle(lstyle);
  mline.SetLineStyle(lstyle);
  //double y_new=(y1+y2)/2.;
  //mline.DrawLineNDC(x1,y_new,x2,y_new);
  mline.DrawLineNDC(x1,y1,x2,y1);
  mline.DrawLineNDC(x1,y2,x2,y2);
  mline.DrawLineNDC(x1,y1,x1,y2);
  mline.DrawLineNDC(x2,y1,x2,y2);

  if (text[0] != '\0') {
    TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize);
    l.SetNDC();
    l.DrawLatex(x,y,text);
  }
}




TBox* TBoxNDC (const double x1, const double y1, const double x2, const double y2) {
  TPad* p;
  if (gDirectory->Get ("box_pad"))
    p = (TPad*)gDirectory->Get ("box_pad");
  else {
    p = new TPad ("box_pad", "box_pad", 0., 0., 1., 1.);
    p->SetFillStyle (0);
  } 
  p->Draw ();
  p->cd ();
  TBox* b = new TBox (x1, y1, x2, y2);
  return b;
}




void myMarkerAndBoxAndLineText (double x, double y, const double bsize, const int bstyle, const int bcolor, const double balpha, const int mcolor, const int mstyle, const double msize, const char* text, const double tsize) {
  const double y1 = y - (0.25*tsize) - (0.004*bsize) + 0.25*tsize;
  const double y2 = y + (0.25*tsize) + (0.004*bsize) + 0.25*tsize;
  const double x2 = x - (0.8*tsize) + 0.02;
  const double x1 = x - (0.8*tsize) + 0.02 - (0.04*bsize);

  TPave *mbox= new TPave (x1, y1, x2, y2, 0, "NDC");
  mbox->SetFillColorAlpha (bcolor, balpha);
  mbox->SetFillStyle (bstyle);
  mbox->Draw ();

  TLine mline;
  mline.SetLineWidth (1);
  mline.SetLineColor (mcolor);
  mline.SetLineStyle (1);
  mline.DrawLineNDC (x1, y1, x2, y1);
  mline.DrawLineNDC (x1, y2, x2, y2);
  mline.DrawLineNDC (x1, y1, x1, y2);
  mline.DrawLineNDC (x2, y1, x2, y2);

  if (mstyle != -1) {

    TMarker *marker = new TMarker(x - (0.8*tsize)+0.02-0.02*bsize, y+0.25*tsize, 8);
    marker->SetNDC();
    marker->SetMarkerColor (IsOpenMarker (mstyle) ? kBlack : mcolor);
    marker->SetMarkerStyle (mstyle);
    marker->SetMarkerSize (msize);

    //TLine *markerLine = new TLine(marker->GetX()-0.018*msize, marker->GetY(), marker->GetX()+0.018*msize, marker->GetY());
    TLine *markerLine = new TLine ();
    markerLine->SetNDC();
    markerLine->SetLineColor(mcolor);
    markerLine->SetLineStyle(1);
    markerLine->SetLineWidth(2);

    markerLine->DrawLineNDC (0.9*x1+0.1*x2, 0.5*(y1+y2), 0.1*x1+0.9*x2, 0.5*(y1+y2));
    markerLine->DrawLineNDC (0.5*(x1+x2), 0.9*y1+0.1*y2, 0.5*(x1+x2), 0.1*y1+0.9*y2);

    //if (IsOpenMarker (mstyle) {
    //  Rectangle_t bb = marker->GetBBox ();

    //  const double xbbl_user = gPad->PixelToUser (bb.fX);
    //  const short xbbl = bb.fX- 0.5*bb.fWidth ();
    //  const short ybbl = bb.fY + 0.5*bb.fWidth ();

    //  
    //}

    marker->Draw();

    if (FullToOpenMarker (mstyle) != kDot) {
      TMarker *marker2 = new TMarker(x - (0.8*tsize)+0.02-0.02*bsize, y+0.25*tsize, 8);
      marker2->SetNDC();
      marker2->SetMarkerColor (kBlack);
      marker2->SetMarkerStyle (FullToOpenMarker (mstyle));
      marker2->SetMarkerSize (msize);
      marker2->Draw();
    }
  }
  
  TLatex l;
  l.SetTextAlign (11);
  l.SetTextSize (tsize);
  l.SetNDC ();
  l.DrawLatex (x, y, text);
}
