#include "DataSet_Parameters.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataSet_Parameters::DataSet_Parameters() {}

void DataSet_Parameters::Debug() const {
  mprintf("Atom Types:\n");
  for (AtomTypeArray::const_iterator at = atomTypes_.begin(); at != atomTypes_.end(); ++at)
    mprintf("\t%s %i\n", *(at->first), at->second);
  mprintf("Bond parameters:\n");
  for (ParmHolder<BondParmType>::const_iterator bp = bondParm_.begin(); bp != bondParm_.end(); ++bp)
    mprintf("\t%s - %s : %f %f\n", *(bp->first[0]), *(bp->first[1]), bp->second.Rk(), bp->second.Req());
  mprintf("Angle parameters:\n");
  for (ParmHolder<AngleParmType>::const_iterator bp = angleParm_.begin(); bp != angleParm_.end(); ++bp)
    mprintf("\t%s - %s - %s : %f %f\n", *(bp->first[0]), *(bp->first[1]), *(bp->first[2]), bp->second.Tk(), bp->second.Teq());
  mprintf("UB parameters:\n");
  for (ParmHolder<BondParmType>::const_iterator bp = ubParm_.begin(); bp != ubParm_.end(); ++bp)
    mprintf("\t%s - %s : %f %f\n", *(bp->first[0]), *(bp->first[1]), bp->second.Rk(), bp->second.Req());

}
  
