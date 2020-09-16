#include "PotentialTerm_Angle.h"
#include "Topology.h"
#include "CharMask.h"
#include "EnergyArray.h"

/** Set up Hooke's law angle term. */
int PotentialTerm_Angle::SetupTerm(Topology const& topIn, CharMask const& maskIn,
                                  EnergyArray& Earray)
{
  activeAngs_.clear();
  for (AngleArray::const_iterator ang = topIn.Angles().begin(); ang != topIn.Angles().end(); ++ang)
  {
    if (maskIn.AtomInCharMask( ang->A1() ) ||
        maskIn.AtomInCharMask( ang->A2() ) ||
        maskIn.AtomInCharMask( ang->A3() ))
    {
      //mprintf("DEBUG: Angle %i to %i to %i\n", ang->A1()+1, ang->A2()+1, ang->A3()+1);
      activeAngs_.push_back( *ang );
    }
  }
  angParm_ = &(topIn.AngleParm());
  Eang_ = Earray.AddType( EnergyArray::E_ANGLE );

  return 0;
}

