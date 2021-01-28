#ifdef HAS_OPENMM
#include "OMM_helpers.h"
#include "OpenMM.h"
/** Add BondArray to openmm */
void Cpptraj::OMM::AddBonds(OpenMM::HarmonicBondForce* bondStretch,
                                    OpenMM::System* system,
                                    std::vector< std::pair<int,int> >& bondPairs,
                                    BondArray const& bonds, BondParmArray const& BP,
                                    std::vector<int> const& oldToNew,
                                    bool useConstraints)
{
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    int a1 = oldToNew[bnd->A1()];
    int a2 = oldToNew[bnd->A2()];
    if (a1 != -1 && a2 != -1)
    {
      bondPairs.push_back(std::make_pair(a1, a2));
      if (useConstraints) {
        system->addConstraint( a1, a2, BP[bnd->Idx()].Req() * OpenMM::NmPerAngstrom );
      } else {
         // Note factor of 2 for stiffness below because Amber specifies the constant
         // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
         // it as used in the force term kx, with energy kx^2/2.
        bondStretch->addBond( a1, a2,
                              BP[bnd->Idx()].Req() * OpenMM::NmPerAngstrom,
                              BP[bnd->Idx()].Rk() * 2 * OpenMM::KJPerKcal
                                * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm );
      }
    }
  }
}

/** Add AngleArray to openmm. */
void Cpptraj::OMM::AddAngles(OpenMM::HarmonicAngleForce* angleStretch,
                                     AngleArray const& angles, AngleParmArray const& AP,
                                     std::vector<int> const& oldToNew)
{
  for (AngleArray::const_iterator ang = angles.begin(); ang != angles.end(); ++ang)
  {
    int a1 = oldToNew[ang->A1()];
    int a2 = oldToNew[ang->A2()];
    int a3 = oldToNew[ang->A3()];
    if (a1 != -1 && a2 != -1 && a3 != -1)
    {
      // CPPTRAJ angles are already in radians, no need for OpenMM::RadiansPerDegree
      // Note factor of 2 for stiffness below because Amber specifies the constant
      // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants
      // it as used in the force term kx, with energy kx^2/2.
      angleStretch->addAngle( a1, a2, a3,
                            AP[ang->Idx()].Teq(),
                            AP[ang->Idx()].Tk() * 2 * OpenMM::KJPerKcal);
    }
  }
}

/** Add DihedralArray to openmm */
void Cpptraj::OMM::AddDihedrals(OpenMM::PeriodicTorsionForce* ptorsion,
                                     DihedralArray const& dihedrals, DihedralParmArray const& DP,
                                     std::vector<int> const& oldToNew)
{
  for (DihedralArray::const_iterator dih = dihedrals.begin(); dih != dihedrals.end(); ++dih)
  {
    int a1 = oldToNew[dih->A1()];
    int a2 = oldToNew[dih->A2()];
    int a3 = oldToNew[dih->A3()];
    int a4 = oldToNew[dih->A4()];
    if (a1 != -1 && a2 != -1 && a3 != -1 && a4 != -1)
    {
      // CPPTRAJ dihedrals are already in radians, no need for OpenMM::RadiansPerDegree
      ptorsion->addTorsion( a1, a2, a3, a4,
                            DP[dih->Idx()].Pn(),
                            DP[dih->Idx()].Phase(),
                            DP[dih->Idx()].Pk() * OpenMM::KJPerKcal);
    }
  }
}

#endif /* HAS_OPENMM */
