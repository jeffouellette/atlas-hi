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

bool IsOpenMarker (const Style_t ms);

bool IsFullMarker (const Style_t ms);

Style_t FullToOpenMarker (const Style_t ms);

Style_t OpenToFullMarker (const Style_t ms);

void SetStyle ();

void FormatTH2Canvas (TCanvas* c, const bool zAxisSpace=true);

void ATLAS_LABEL (double x, double y, Color_t color=1); 

TGraphErrors* myTGraphErrorsDivide (TGraphErrors* g1, TGraphErrors* g2);

TGraphAsymmErrors* myTGraphErrorsDivide (TGraphAsymmErrors* g1, TGraphAsymmErrors* g2);

TGraphAsymmErrors* myMakeBand (TGraphErrors* g0, TGraphErrors* g1, TGraphErrors* g2);

void myAddtoBand (TGraphErrors* g1, TGraphAsymmErrors* g2);

TGraphErrors* TH1TOTGraph (TH1 *h1);

void BinomialDivide (TH1D* n, TH1D* d);

void myText (double x, double y,  Color_t color, const char *text, double tsize=0.04);

void myBoxText (double x, double y, double boxsize, int mcolor, const char *text);

void myLineText (double x, double y, int color, int lstyle, const char *text, float lsize, double tsize);

void myLineColorText (double x, double y, int color, int lstyle, const char *text, float lsize, double tsize);

void myMarkerText (double x, double y, int color, int mstyle, const char *text, float msize=1.25, double tsize=0.032);

void myMarkerTextNoLine (double x, double y, int color, int mstyle, const char *text, float msize=1.25, double tsize=0.032);

void myOnlyBoxText (double x, double y, double boxsize, int mcolor, int lcolor, int lstyle, const char *text, double tsize, int bstyle, double balpha);

TBox* TBoxNDC (const double x1, const double y1, const double x2, const double y2);

void myMarkerAndBoxAndLineText (double x, double y, const double bsize, const int bstyle, const int bcolor, const double balpha, const int mcolor, const int mstyle, const double msize, const char* text, const double tsize=0.032);

#endif // __ATLASUTILS_H
