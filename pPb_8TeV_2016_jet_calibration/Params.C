#include "../GlobalParams.C"

/** User defined parameters **/

const bool runPeriodA = true; // Analyze period A data
const bool runPeriodB = true; // Analyze period B data

/** End user defined parameters **/


/** General (non-user defined) paramters **/

const double pbins[15] = {20., 30., 40., 50., 60., 75., 90., 105., 120., 140., 160., 180., 200., 250., 300.};
const int numpbins = sizeof(pbins)/sizeof(pbins[0]) - 1;
const double etabins[7] = {-1.3, -0.8, -0.3, 0, 0.3, 0.8, 1.3};
const int numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;

/** End general parameters **/
