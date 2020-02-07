#ifndef __header_h__
#define __header_h__

#include <string>

using namespace std;

namespace pullstudy {

float GetPtHatMin (short iPtZ) {
  switch (iPtZ) {
    case 2: return 5;
    case 3: return 15;
    case 4: return 30;
    default: return 0;
  }
}


float GetMinZPt (short iPtZ) {
  switch (iPtZ) {
    case 2: return 15;
    case 3: return 30;
    case 4: return 60;
    default: return 0;
  }
}


float GetMaxZPt (short iPtZ) {
  switch (iPtZ) {
    case 2: return 30;
    case 3: return 60;
    case 4: return 10000;
    default: return 0;
  }
}


int GetMinCent (short iCent) {
  switch (iCent) {
    case 1: return 80;
    case 2: return 30;
    case 3: return 10;
    default: return -1;
  }
}


int GetMaxCent (short iCent) {
  switch (iCent) {
    case 1: return 30;
    case 2: return 10;
    case 3: return 0;
    default: return -1;
  }
}


string GetMinTrkStr (short iX, bool useTrkPt) {
  switch (iX) {
    case 0: return useTrkPt ? "1" : "1/60";
    case 1: return useTrkPt ? "2" : "1/30";
    case 2: return useTrkPt ? "4" : "1/15";
    case 3: return useTrkPt ? "8" : "1/8";
    case 4: return useTrkPt ? "15" : "1/4";
    case 5: return useTrkPt ? "30" : "1/2";
    default: return "";
  }
}


string GetMaxTrkStr (short iX, bool useTrkPt) {
  switch (iX) {
    case 0: return useTrkPt ? "2" : "1/30";
    case 1: return useTrkPt ? "4" : "1/15";
    case 2: return useTrkPt ? "8" : "1/8";
    case 3: return useTrkPt ? "15" : "1/4";
    case 4: return useTrkPt ? "30" : "1/2";
    case 5: return useTrkPt ? "60" : "1";
    default: return "";
  }
}


int GetNEvents (short iPtZ, short iCent) {
  if (iPtZ == 2 && iCent == 1) return 657;
  if (iPtZ == 2 && iCent == 2) return 1501;
  if (iPtZ == 2 && iCent == 3) return 1587;
  if (iPtZ == 3 && iCent == 1) return 354;
  if (iPtZ == 3 && iCent == 2) return 848;
  if (iPtZ == 3 && iCent == 3) return 849;
  if (iPtZ == 4 && iCent == 1) return 141;
  if (iPtZ == 4 && iCent == 2) return 318;
  if (iPtZ == 4 && iCent == 3) return 318;
  return 0;
}


string FormatCounts (int counts) {
  if (counts < 1000) return "";
  else if (1000 <= counts && counts < 10000) {
    string countsStr = atlashi::FormatMeasurement (counts, 0, 1);
    countsStr = countsStr.substr(0, 1) + "k";
    return countsStr;
  }
  else if (10000 <= counts && counts < 100000) {
    string countsStr = atlashi::FormatMeasurement (counts, 0, 2);
    countsStr = countsStr.substr(0, 2) + "k";
    return countsStr;
  }
  else if (100000 <= counts && counts < 1000000) {
    string countsStr = atlashi::FormatMeasurement (counts, 0, 3);
    countsStr = countsStr.substr(0, 3) + "k";
    return countsStr;
  }
  else return "";
}


int GetNumInGroup1 (short iPtZ, short iCent) {
  if (iPtZ == 2) {
    if (iCent == 1)       return 260;
    else if (iCent == 2)  return 673;
    else if (iCent == 3)  return 782;
  }
  else if (iPtZ == 3) {
    if (iCent == 1)       return 148;
    else if (iCent == 2)  return 364;
    else if (iCent == 3)  return 443;
  }
  else if (iPtZ == 4) {
    if (iCent == 1)       return 52;
    else if (iCent == 2)  return 147;
    else if (iCent == 3)  return 145;
  }
  return 0;
}


int GetNumInGroup2 (short iPtZ, short iCent) {
  if (iPtZ == 2) {
    if (iCent == 1)       return 397;
    else if (iCent == 2)  return 828;
    else if (iCent == 3)  return 805;
  }
  else if (iPtZ == 3) {
    if (iCent == 1)       return 206;
    else if (iCent == 2)  return 484;
    else if (iCent == 3)  return 406;
  }
  else if (iPtZ == 4) {
    if (iCent == 1)       return 89;
    else if (iCent == 2)  return 171;
    else if (iCent == 3)  return 173;
  }
  return 0;
}


} // end namespace pullstudy


#endif
