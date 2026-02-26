#include "DihedralParmSet.h"
#include "ParameterSet.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Parm;

/// Print all dihedral multiplicities to stdout
static inline void printdih(DihedralParmArray const& dpa)
{
  for (DihedralParmArray::const_iterator it = dpa.begin(); it != dpa.end(); ++it)
    mprintf("Warning:\t\t%3g PK= %g Phase= %g\n", it->Pn(), it->Pk(), it->Phase()*Constants::RADDEG);
}

/** Add all dihedral parameters into the given parameter set. */
int DihedralParmSet::ToDihParm(ParameterSet& prm) const
{
  for (DihedralParmSet::const_iterator it = begin(); it != end(); ++it)
  {
    // Check multiplicites; warn if one seems missing
    if (it->second.size() > 1) {
      DihedralParmArray::const_iterator dp = it->second.begin();
      unsigned int next_expected_mult = (unsigned int)dp->Pn() + 1;
      dp++;
      for (; dp != it->second.end(); dp++) {
        if ((unsigned int)dp->Pn() != next_expected_mult)
          mprintf("Warning: %s does not have dihedral term for multiplicity %u\n",
                  it->first.TypeNameStr("dihedral").c_str(), next_expected_mult);
        next_expected_mult = (unsigned int)dp->Pn() + 1;
      }
    }
    Cpptraj::Parm::RetType ret = prm.DP().AddParm( it->first, it->second, true );
    if (ret == Cpptraj::Parm::SAME)
      mprintf("Warning: Duplicated %s\n", it->first.TypeNameStr("dihedral").c_str());
    else if (ret == Cpptraj::Parm::UPDATED) {// TODO SCEE/SCNB?
      mprintf("Warning: Redefining %s from\n", it->first.TypeNameStr("dihedral").c_str());
      printdih(prm.DP().PreviousArray());
      mprintf("Warning: to\n");
      printdih(it->second);
    } else if (ret == Cpptraj::Parm::ERR) {
      mprinterr("Error: Reading %s\n", it->first.TypeNameStr("dihedral").c_str());
      return 1;
    }
  }
  return 0;
}
