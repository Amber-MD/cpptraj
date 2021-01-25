#include "PotentialTerm_Dihedral.h"
#include "Topology.h"
#include "CharMask.h"
#include "EnergyArray.h"
#include "CpptrajStdio.h" 

void PotentialTerm_Dihedral::addDihedrals(DihedralArray const& dihedrals, CharMask const& maskIn)
{
  for (DihedralArray::const_iterator dih = dihedrals.begin(); dih != dihedrals.end(); ++dih)
  {
    if (maskIn.AtomInCharMask( dih->A1() ) ||
        maskIn.AtomInCharMask( dih->A2() ) ||
        maskIn.AtomInCharMask( dih->A3() ) ||
        maskIn.AtomInCharMask( dih->A4() ))
    {
      mprintf("DEBUG: Dihedral %i to %i to %i to %i\n", dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1);
      activeDihs_.push_back( *dih );
    }
  }
}

/** Set up truncated Fourier series dihedral term. */
int PotentialTerm_Dihedral::SetupTerm(Topology const& topIn, Box const& boxIn,
                                      CharMask const& maskIn, EnergyArray& Earray)
{
  activeDihs_.clear();
  addDihedrals( topIn.Dihedrals(),  maskIn );
  addDihedrals( topIn.DihedralsH(), maskIn );

  dihParm_ = &(topIn.DihedralParm());
  Edih_ = Earray.AddType( EnergyArray::E_DIHEDRAL );

  return 0;
}

/** Calculate truncated Fourier series dihedral term. 
  * NOTE: Code adapted from AmberTools 20 SFF, sff2.c ephi2()
  */
void PotentialTerm_Dihedral::CalcForce(Frame& frameIn, CharMask const& maskIn) const {
  *Edih_ = 0.0;

}
