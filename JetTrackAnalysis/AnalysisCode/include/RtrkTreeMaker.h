#ifndef __RtrkTreeMaker_h__
#define __RtrkTreeMaker_h__

namespace JetTrackAnalysis {

void RtrkTreeMaker (const char* directory,
                    const int dataSet,
                    const bool isMC,
                    const bool isPeriodA,
                    const char* inFileName = "",
                    const double crossSection_microbarns = 0,
                    const double filterEfficiency = 0,
                    const int numberEvents = 0);

} // end namespace

#endif
