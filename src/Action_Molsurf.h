#ifndef INC_ACTION_MOLSURF_H
#define INC_ACTION_MOLSURF_H
#include "Action.h"
#include "molsurf.h"
// Class: Molsurf
/// Wrapper for the molsurf routine in molsurf.c
/** This is the cpptraj wrapper for the molsurf routine originally written
  * by Paul Beroza. This calculates the Connolly surface area of the
  * molecule. See "M.L. Connolly, Analytical molecular surface calculation,
  * J. Appl. Cryst., 16, p. 548-558, 1983."
  */
class Molsurf: public Action {
    DataSet *sasa;
    AtomMask Mask1;
    ATOM *atom;
    double probe_rad;
  public:
    Molsurf();
    ~Molsurf();

    int init();
    int setup();
    int action();
};
#endif
