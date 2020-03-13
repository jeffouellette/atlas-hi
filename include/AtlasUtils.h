//
//   @file    AtlasUtils.h         
//   
//
//   @author M.Sutton
// 
//   Copyright (C) 2010 Atlas Collaboration
//
//   $Id: AtlasUtils.h, v0.0   Thu 25 Mar 2010 10:34:20 CET $


#ifndef __ATLASUTILS_H
#define __ATLASUTILS_H

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TCanvas.h"

void SetStyle ();

void FormatTH2Canvas (TCanvas* c, const bool zAxisSpace=true);

void ATLAS_LABEL(Double_t x,Double_t y,Color_t color=1); 

TGraphErrors* myTGraphErrorsDivide(TGraphErrors* g1,TGraphErrors* g2);

TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);

TGraphAsymmErrors* myMakeBand(TGraphErrors* g0, TGraphErrors* g1,TGraphErrors* g2);

void myAddtoBand(TGraphErrors* g1, TGraphAsymmErrors* g2);

TGraphErrors* TH1TOTGraph(TH1 *h1);

void BinomialDivide (TH1D* n, TH1D* d);

void myText(Double_t x,Double_t y,Color_t color,const char *text,Double_t tsize=0.04);

void myBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,const char *text);

void myLineText(Double_t x,Double_t y,Int_t color, Int_t lstyle,const char *text,Float_t lsize,Double_t tsize);

void myLineColorText(Double_t x,Double_t y,Int_t color, Int_t lstyle,const char *text,Float_t lsize,Double_t tsize);

void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle,const char *text,Float_t msize=1.25,Double_t tsize=0.032);

void myMarkerTextNoLine(Double_t x,Double_t y,Int_t color,Int_t mstyle,const char *text,Float_t msize=1.25,Double_t tsize=0.032);

void myOnlyBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,Int_t lcolor,Int_t lstyle, const char *text, Double_t tsize, Int_t bstyle, Double_t balpha);

TBox* TBoxNDC (const double x1, const double y1, const double x2, const double y2);

#endif // __ATLASUTILS_H
