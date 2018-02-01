#include "params.C"
#include <string>

/**
 * Trigger stores information about a trigger, including a jet momenta range, pseudorapidity interval,
 * name, and a (unique) branching index for event analysis.
 */
class Trigger {
    
    public:
    string name;

    int min_pt;
    double lower_eta;    
    double upper_eta;
    int lowerRunNumber;
    int upperRunNumber;
    int index;
    bool disabled;
    bool isBootstrapped;
    Trigger* referenceTrigger;

    bool m_trig_bool;
    float m_trig_prescale;

    Trigger(string, int, double, double, int, int);
    Trigger(string, int, double, double, bool, int, int);
    Trigger(const Trigger* t);

};

/**
 * Creates a Trigger object. By default, the maximum momentum and branching index are both 0. It is
 * expected that these values will be nonzero by the time the object is used purposefully.
 */
Trigger::Trigger(string thisname, int thismin_pt, double etal, double etau, int lRN=0, int uRN=INT_MAX) {
    name = thisname;
    min_pt = thismin_pt;
    lower_eta = etal;
    upper_eta = etau;
    lowerRunNumber = lRN;
    upperRunNumber = uRN;
    index = 0;
    disabled = false;
    isBootstrapped = false;
    referenceTrigger = NULL;
}

/**
 * Creates a Trigger object. By default, the maximum momentum and branching index are both 0. It is
 * expected that these values will be nonzero by the time the object is used purposefully.
 */
Trigger::Trigger(string thisname, int thismin_pt, double etal, double etau, bool thisdisabled, int lRN=0, int uRN=INT_MAX) {
    name = thisname;
    min_pt = thismin_pt;
    lower_eta = etal;
    upper_eta = etau;
    lowerRunNumber = lRN;
    upperRunNumber = uRN;
    index = 0;
    disabled = thisdisabled;
    isBootstrapped = false;
    referenceTrigger = NULL;
}

/**
 * Creates a copy of trigger t.
 */
Trigger::Trigger(const Trigger* t) {
    name = t->name;
    min_pt = t->min_pt;
    lower_eta = t->lower_eta;
    upper_eta = t->upper_eta;
    lowerRunNumber = t->lowerRunNumber;
    upperRunNumber = t->upperRunNumber;
    index = t->index;
    disabled = t->disabled;
    isBootstrapped = t->isBootstrapped;
    referenceTrigger = t->referenceTrigger;
}
