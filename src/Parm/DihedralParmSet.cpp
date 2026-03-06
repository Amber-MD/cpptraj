#include "DihedralParmSet.h"
#include "ParameterSet.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Parm;

/// Print all multiplicities to stdout
static inline void printdih(DihedralParmArray const& dpa)
{
  for (DihedralParmArray::const_iterator it = dpa.begin(); it != dpa.end(); ++it)
    mprintf("Warning:\t\t%3g PK= %g Phase= %g\n", it->Pn(), it->Pk(), it->Phase()*Constants::RADDEG);
}

/** Print multiplicity warning */
void DihedralParmSet::check_mult(std::string const& desc, TypeNameHolder const& types, DihedralParmArray const& DPA) {
  if (DPA.size() > 1) {
    DihedralParmArray::const_iterator dp = DPA.begin();
    unsigned int next_expected_mult = (unsigned int)dp->Pn() + 1;
    dp++;
    for (; dp != DPA.end(); dp++) {
      if ((unsigned int)dp->Pn() != next_expected_mult)
        mprintf("Warning: %s does not have term for multiplicity %u\n",
                types.TypeNameStr(desc).c_str(), next_expected_mult);
      next_expected_mult = (unsigned int)dp->Pn() + 1;
    }
  }
}

/** Add all dihedral parameters into the given dihedral parameter holder. */
int DihedralParmSet::ToDihParm(DihedralParmHolder& DP) const
{
  for (DihedralParmSet::const_iterator it = begin(); it != end(); ++it)
  {
    // Check multiplicites; warn if one seems missing
    if (debug_ > 0)
      check_mult("dihedral", it->first, it->second);

    Cpptraj::Parm::RetType ret = DP.AddParm( it->first, it->second, true );
    if (ret == Cpptraj::Parm::SAME)
      mprintf("Warning: Duplicated %s\n", it->first.TypeNameStr("dihedral").c_str());
    else if (ret == Cpptraj::Parm::UPDATED) {// TODO SCEE/SCNB?
      mprintf("Warning: Redefining %s from\n", it->first.TypeNameStr("dihedral").c_str());
      printdih(DP.PreviousParm());
      mprintf("Warning: to\n");
      printdih(it->second);
    } else if (ret == Cpptraj::Parm::ERR) {
      mprinterr("Error: Reading %s\n", it->first.TypeNameStr("dihedral").c_str());
      return 1;
    }
  }
  return 0;
}

/** Add all improper parameters into the given improper parameter holder. */
int DihedralParmSet::ToImpParm(ImproperParmHolder& IP) const
{
  for (DihedralParmSet::const_iterator it = begin(); it != end(); ++it)
  {
    // Check multiplicites; warn if one seems missing
    if (debug_ > 0)
      check_mult("improper", it->first, it->second);

    Cpptraj::Parm::RetType ret = IP.AddParm( it->first, it->second, true );
    if (ret == Cpptraj::Parm::SAME)
      mprintf("Warning: Duplicated %s\n", it->first.TypeNameStr("improper").c_str());
    else if (ret == Cpptraj::Parm::UPDATED) {
      mprintf("Warning: Redefining %s from\n", it->first.TypeNameStr("improper").c_str());
      printdih(IP.PreviousParm());
      mprintf("Warning: to\n");
      printdih(it->second);
    } else if (ret == Cpptraj::Parm::ERR) {
      mprinterr("Error: Reading %s\n", it->first.TypeNameStr("improper").c_str());
      return 1;
    }
  }
  return 0;
}
