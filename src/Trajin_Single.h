#ifndef TRAJIN_SINGLE_H
#define TRAJIN_SINGLE_H
#include "Trajin.h"
class Trajin_Single : public Trajin {
  public:
    Trajin_Single();
    ~Trajin_Single();

    int SetupTrajRead(std::string const&, ArgList *, Topology *);
  private:
    TrajectoryIO* trajio_;   ///< Hold class that will interface with traj format.  
};
#endif
