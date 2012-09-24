#ifndef INC_TRAJIN_MULTI_H
#define INC_TRAJIN_MULTI_H
#include "Trajin.h"
class Trajin_Multi : public Trajin {
  public:
    Trajin_Multi();
    ~Trajin_Multi();

    int SetupTrajRead(std::string const&, ArgList *, Topology *);
  private:
    /// Define type that will hold REMD indices
    typedef std::vector<int> RemdIdxType;
    typedef std::vector<TrajectoryIO*> IOarrayType;
    enum TargetType { TEMP = 0, INDICES };

    double remdtrajtemp_;                    ///< Get frames with this temperature on read
    RemdIdxType remdtrajidx_;                ///< Get frames with these indices on read
    int* remd_indices_;                      ///< Space for reading in REMD indices.
    IOarrayType REMDtraj_;                   ///< Input replica trajectories
    int lowestRepnum_;                       ///< Hold the lowest replica number
    TargetType targetType_;

    std::vector<std::string> SearchForReplicas(bool);
};
#endif
