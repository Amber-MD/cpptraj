#include "DataSet_Topology.h"
#include "ParmFile.h"
#include "CpptrajStdio.h"

/** This routine assumes that MetaData has been set.*/
int DataSet_Topology::LoadTopFromFile(ArgList const& argIn, int debugIn)
{
  if (Meta().Fname().empty()) {
    mprinterr("Internal Error: Topology DataSet file name has not been set.\n");
    return 1;
  }
  top_.SetDebug( debugIn );
  ParmFile pfile;
  if (pfile.ReadTopology(top_, Meta().Fname(), argIn, debugIn)) {
    mprinterr("Error: Could not open topology '%s'\n", Meta().Fname().full());
    return 1;
  }
  return 0;
}

int DataSet_Topology::StripTop( std::string const& maskexpr ) {
  AtomMask tempMask( maskexpr );
  // Since want to keep atoms outside mask, invert selection
  tempMask.InvertMaskExpression();
  if (top_.SetupIntegerMask( tempMask )) return 1;
  mprintf("\tStripping atoms in mask [%s] (%i) from %s\n",tempMask.MaskString(),
           top_.Natom() - tempMask.Nselected(), legend());
  Topology* tempParm = top_.modifyStateByMask(tempMask);
  if (tempParm==0) { 
    mprinterr("Error: Could not strip parm '%s'.\n", legend());
    return 1;
  } else { 
    // Replace parm with stripped version
    top_ = *tempParm;
    top_.Brief("Stripped parm:");
    delete tempParm;
  }
  return 0;
}
