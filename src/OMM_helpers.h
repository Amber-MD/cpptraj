#ifndef INC_OMM_HELPERS_H
#define INC_OMM_HELPERS_H
#ifdef HAS_OPENMM
#include "ParameterTypes.h"
#include <vector>
// Forward declares
namespace OpenMM {
  class System;
  class Context;
  class HarmonicBondForce;
  class HarmonicAngleForce;
  class PeriodicTorsionForce;
}

namespace Cpptraj {
/** Contains functions used to interface between cpptraj and openmm. */
namespace OMM {

/// Add bonds to a HarmonicBondForce
void AddBonds(OpenMM::HarmonicBondForce*, OpenMM::System*,
              std::vector< std::pair<int,int> >&,
              BondArray const&, BondParmArray const&, std::vector<int> const&,
              bool);
/// Add angles to a HarmonicAngleForce
void AddAngles(OpenMM::HarmonicAngleForce*, AngleArray const&, AngleParmArray const&, 
               std::vector<int> const&);
/// Add dihedrals to a PeriodicTorsionForce
void AddDihedrals(OpenMM::PeriodicTorsionForce*, DihedralArray const&,
                  DihedralParmArray const&, std::vector<int> const&);

} // END namespace OMM
} // END namespace Cpptraj
#endif /* HAS_OPENMM */
#endif
