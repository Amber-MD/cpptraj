#include "Exec_Random.h"
#include "CpptrajStdio.h"
#include "Random.h"

// Exec_Random::Help()
void Exec_Random::Help() const
{
  mprintf(
# ifdef C11_SUPPORT
          "\t[setdefault {marsaglia|stdlib|mt|pcg32|xo128}]\n"
# else
          "\t[setdefault {marsaglia|stdlib|pcg32|xo128}]\n"
# endif
          "\t[createset <name> count <#> settype {int} [seed <#>]]\n");
}

// Exec_Random::Execute()
Exec::RetType Exec_Random::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string setarg = argIn.GetStringKey("setdefault");
  if (!setarg.empty()) {
    if (setarg == "marsaglia") Random_Number::SetDefaultRng( Random_Number::MARSAGLIAS );
    else if (setarg == "stdlib") Random_Number::SetDefaultRng( Random_Number::STDLIB );
    else if (setarg == "mt") {
#     ifdef C11_SUPPORT
      Random_Number::SetDefaultRng( Random_Number::MERSENNE_TWISTER );
#     else
      mprinterr("Error: Mersenne twister RNG requires C++11 support.\n");
      return CpptrajState::ERR;
#     endif
    } else if (setarg == "pcg32") Random_Number::SetDefaultRng( Random_Number::PCG32 );
    else if (setarg == "xo128") Random_Number::SetDefaultRng( Random_Number::XOSHIRO128PP );
    else {
      mprinterr("Error: Unrecognized RNG type for 'setdefault': %s\n", setarg.c_str());
      return CpptrajState::ERR;
    }
    mprintf("\tDefault RNG set to '%s'\n", Random_Number::CurrentDefaultRngStr());
  }

  std::string dsname = argIn.GetStringKey("createset");
  if (!dsname.empty()) {
    int iseed = argIn.getKeyInt("seed", -1);
    int count = argIn.getKeyInt("count", -1);
    if (count < 1) {
      mprinterr("Error: Must specify 'count' > 0 for 'createset'\n");
      return CpptrajState::ERR;
    }
    DataFile* outfile = State.DFL().AddDataFile( argIn.GetStringKey("out"), argIn );

    Random_Number rng;
    if (rng.rn_set( iseed )) return CpptrajState::ERR;
    std::string typestr = argIn.GetStringKey("settype");
    if (typestr == "int") {
      // Create integer set
      DataSet* ds = State.DSL().AddSet( DataSet::UNSIGNED_INTEGER, MetaData(dsname) );
      if (ds == 0) return CpptrajState::ERR;
      if (outfile != 0) outfile->AddDataSet( ds );
      for (int idx = 0; idx != count; idx++) {
        unsigned int rn = rng.rn_num();
        ds->Add(idx, &rn);
      }
    } else {
      mprinterr("Error: Unrecognized 'settype' for 'createset': %s\n", typestr.c_str());
      return CpptrajState::ERR;
    }
  }

  return CpptrajState::OK;
}
