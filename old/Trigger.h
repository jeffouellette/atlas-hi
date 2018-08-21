/**
 * Implements a trigger class inspired by a linked list design. A trigger points at its reference trigger, which is used to define its efficiency factor.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

#include <string>

/**
 * Trigger stores information about a trigger, including a jet momenta range, pseudorapidity interval,
 * name, and a (unique) branching index for event analysis.
 */

#ifndef __Trigger_h__
#define __Trigger_h__

struct Trigger {
    
  public:
   string name;

   int minPt;
   int maxPt;
   int thresPt;
   double lowerEta;    
   double upperEta;
   int lowerRunNumber;
   int upperRunNumber;
   int index;
   bool disabled;
   bool isBootstrapped;
   Trigger* referenceTrigger;

   bool trigBool;
   float trigPrescale;

   Trigger(const string _name, const int _thresPt, const double _lowerEta, const double _upperEta, const int _lRN=0, const int _uRN=10000000, const int _minPt=0);
   Trigger(const string _name, const int _thresPt, const double _lowerEta, const double _upperEta, const bool _disabled, const int _lRN=0, const int _uRN=10000000, const int _minPt=0);
   Trigger(const Trigger* t);

};

#endif
