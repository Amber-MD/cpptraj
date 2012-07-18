#ifndef INC_ACTION_MRT_H
#define INC_ACTION_MRT_H
#include "Action.h"
/** Perform mean residence time calculations.
  * \author Original code by: Hannes H. Loeffler
  *         2005-2007: Lab. of Molecular Design, ICMB, the University of Tokyo, Japan
  *         2008-2010: STFC Daresbury Laboratory, UK
  * \author Adapted to C++ by Daniel R. Roe.
  */
class Action_MRT : public Action {
  public:
    Action_MRT();
  private:
    int init();

    CpptrajFile outfile_;
    double time_; // darg1
    int nStar_;
    double lowerCutoff2_;
    double upperCutoff2_;
    std::string autoCorr_;
};
#endif
