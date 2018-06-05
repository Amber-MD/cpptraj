#include "Parm_CharmmParam.h"
#include "ParameterSet.h"
#include "CharmmParamFile.h"

int Parm_CharmmParam::WriteParm(FileName const& fname, Topology const& topOut) {
  // Convert topology to parameter set
  ParameterSet pout;

  for (Topology::atom_iterator it = topOut.begin(); it != topOut.end(); ++it)
  {
    pout.AT().AddAtomType( it->Type(), AtomType( it->Mass() ) );
  }

  // Write out
  CharmmParamFile outfile;
  return outfile.WriteParams(pout, fname, debug_);
}
