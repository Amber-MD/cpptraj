#include "Image_List_Unit.h"
#include "AtomMask.h"
#include "Topology.h"
#include "CpptrajStdio.h"

/** Set up list to contain units selected by the mask expression. */
int Image::List_Unit::SetupList(Topology const& Parm, std::string const& maskExpression)
{
  AtomMask mask;
  if (mask.SetMaskString(maskExpression)) return 1;
  if ( Parm.SetupIntegerMask( mask ) ) return 1;
  mask.MaskInfo();
  if (mask.None()) return -1;
  std::vector<int> molnums = Parm.MolnumsSelectedBy( mask );
  for (std::vector<int>::const_iterator it = molnums.begin(); it != molnums.end(); ++it)
  {
    units_.push_back( Parm.Mol(*it).MolUnit() );
  }
  return 0;
}

/** Add given unit to the list. */
void Image::List_Unit::AddUnit(Unit const& unitIn) {
  units_.push_back( unitIn );
}

/** \return A Unit containing all units in the list. */
Unit Image::List_Unit::AllEntities() const {
  Unit unitOut;
  for (Uarray::const_iterator unit = units_.begin(); unit != units_.end(); ++unit)
    for (Unit::const_iterator seg = unit->segBegin(); seg != unit->segEnd(); ++seg)
      for (int idx = seg->Begin(); idx != seg->End(); ++idx)
        unitOut.AddIndex( idx );
  return unitOut;
}

/** Print molecules to be imaged to STDOUT. */
void Image::List_Unit::PrintEntities() const {
  mprintf("\tImaging molecules defined by atom ranges:\n");
  for (Uarray::const_iterator unit = units_.begin(); unit != units_.end(); ++unit)
  {
    mprintf("\t ");
    for (Unit::const_iterator seg = unit->segBegin(); seg != unit->segEnd(); ++seg)
      mprintf(" %i-%i", seg->Begin()+1, seg->End());
    mprintf("\n");
  }
}
