#include "Trigger.h"

namespace atlashi {

/**
 * Creates a Trigger object. By default, the maximum momentum and branching index are both 0. It is
 * expected that these values will be nonzero by the time the object is used purposefully.
 */
Trigger::Trigger (const string _name, const int _thresPt, const double _lowerEta, const double _upperEta, const int _lRN, const int _uRN, const int _minPt) {
  name = _name;
  thresPt = _thresPt;
  minPt = _minPt;
  maxPt = 0;
  lowerEta = _lowerEta;
  upperEta = _upperEta;
  lowerRunNumber = _lRN;
  upperRunNumber = _uRN;
  index = 0;
  disabled = false;
  isBootstrapped = false;
  referenceTrigger = NULL;
}

/**
 * Creates a Trigger object. By default, the maximum momentum and branching index are both 0. It is
 * expected that these values will be nonzero by the time the object is used purposefully.
 */
Trigger::Trigger (const string _name, const int _thresPt, const double _lowerEta, const double _upperEta, const bool _disabled, const int _lRN, const int _uRN, const int _minPt) {
  name = _name;
  thresPt = _thresPt;
  minPt = _minPt;
  maxPt = 0;
  lowerEta = _lowerEta;
  upperEta = _upperEta;
  lowerRunNumber = _lRN;
  upperRunNumber = _uRN;
  index = 0;
  disabled = _disabled;
  isBootstrapped = false;
  referenceTrigger = NULL;
}

/**
 * Creates a copy of trigger t.
 */
Trigger::Trigger(const Trigger* t) {
  name = t->name;
  thresPt = t->thresPt;
  minPt = t->minPt;
  maxPt = t->maxPt;
  lowerEta = t->lowerEta;
  upperEta = t->upperEta;
  lowerRunNumber = t->lowerRunNumber;
  upperRunNumber = t->upperRunNumber;
  index = t->index;
  disabled = t->disabled;
  isBootstrapped = t->isBootstrapped;
  referenceTrigger = t->referenceTrigger;
}

} // end namespace
