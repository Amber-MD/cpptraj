#ifndef INC_ENERGY_ENE_ELEC_IMAGE_H
#define INC_ENERGY_ENE_ELEC_IMAGE_H
namespace Cpptraj {
namespace Energy {
template <typename T>
T Ene_Elec_Image(Frame const& fIn, Topology const& tIn, AtomMask const& mask,
                                 ExclusionArray const& Excluded,
                                 int n_points)
{
  // Sum over images.
  T Eimage = 0.0;
  // Cache npoints values, excluding this cell (0,0,0)
  std::vector<Vec3> Cells;
  int Ncells = (2*n_points)+1;
  Cells.reserve( (Ncells*Ncells*Ncells) - 1 );
  for (int ix = -n_points; ix <= n_points; ix++)
    for (int iy = -n_points; iy <= n_points; iy++)
      for (int iz = -n_points; iz <= n_points; iz++)
        if (ix != 0 || iy != 0 || iz != 0)
          Cells.push_back( Vec3(ix, iy, iz) );
  // Outer loop over atoms (i)
  for (AtomMask::const_iterator atom1 = mask.begin(); atom1 != mask.end(); ++atom1)
  {
//    mprintf("\nDEBUG: Atom %i\n", *atom1+1);
    Vec3 T1( fIn.XYZ(*atom1) );
    // Inner loop over atoms (j)
    for (AtomMask::const_iterator atom2 = mask.begin(); atom2 != mask.end(); ++atom2)
    {
      Vec3 frac2 = fIn.BoxCrd().FracCell() * Vec3(fIn.XYZ( *atom2 )); // atom j in fractional coords
      T qiqj = Constants::COULOMBFACTOR * tIn[*atom1].Charge() * tIn[*atom2].Charge();
      // Loop over images of atom j
      for (std::vector<Vec3>::const_iterator ixyz = Cells.begin(); ixyz != Cells.end(); ++ixyz)
      {
//        mprintf("DEBUG: Atom %4i to %4i Image %3i %3i %3i", *atom1+1, *atom2+1, ix, iy, iz);
        // atom j image back in Cartesian space minus atom i in Cartesian space.
        Vec3 dxyz = fIn.BoxCrd().UnitCell().TransposeMult(frac2 + *ixyz) - T1;
        T rij2 = dxyz.Magnitude2();
        T rij = sqrt(rij2);
//        mprintf(" Distance= %g\n", rij);
        T e_elec = qiqj / rij;
        Eimage += e_elec;
      }
    } // atom j
  } // atom i
  return (Eimage/2.0);
}
}
}
#endif
