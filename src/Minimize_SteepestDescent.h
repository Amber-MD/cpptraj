#ifndef INC_MINIMIZE_STEEPESTDESCENT_H
#define INC_MINIMIZE_STEEPESTDESCENT_H
#include <string>
class PotentialFunction;
class Frame;
class Minimize_SteepestDescent {
  public:
    Minimize_SteepestDescent();
    /** Set up minimization with optional output trajectory, RMS tolerance, # steps. */
    int SetupMin(std::string const&, double, int);
    /** Run minimization with given potential function and coordinates. */
    int RunMin(PotentialFunction*, Frame&) const;
  private:
    std::string trajoutName_;
    double min_tol_;
    int nMinSteps_;
};
#endif
