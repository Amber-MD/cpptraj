#ifndef INC_MINIMIZE_STEEPESTDESCENT_H
#define INC_MINIMIZE_STEEPESTDESCENT_H
#include <string>
class PotentialFunction;
class Frame;
class CpptrajFile;
class Minimize_SteepestDescent {
  public:
    Minimize_SteepestDescent();
    /** Set up minimization with optional output trajectory, RMS tolerance, step size, # steps. */
    int SetupMin(std::string const&, double, double,int);
    /** Run minimization with given potential function and coordinates. */
    int RunMin(PotentialFunction&, Frame&, CpptrajFile&) const;
  private:
    std::string trajoutName_;
    double min_tol_; ///< Min RMS tolerance
    double dx0_;     ///< Initial step size
    int nMinSteps_;  ///< Number of minimization steps.
};
#endif
