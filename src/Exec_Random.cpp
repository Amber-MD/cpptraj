#include "Exec_Random.h"
#include "CpptrajStdio.h"
#include "Random.h"

// Exec_Random::Help()
void Exec_Random::Help() const
{

}

// Exec_Random::Execute()
Exec::RetType Exec_Random::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string setarg = argIn.GetStringKey("setdefault");
  if (!setarg.empty()) {
    if (setarg == "marsaglia") Random_Number::SetDefaultRng( Random_Number::MARSAGLIAS );
    else if (setarg == "stdlib") Random_Number::SetDefaultRng( Random_Number::STDLIB );
    else if (setarg == "mt") Random_Number::SetDefaultRng( Random_Number::MERSENNE_TWISTER );
    else if (setarg == "pcg32") Random_Number::SetDefaultRng( Random_Number::PCG32 );
    else if (setarg == "xo128") Random_Number::SetDefaultRng( Random_Number::XOSHIRO128PP );
    else {
      mprinterr("Error: Unrecognized RNG type for 'setdefault': %s\n", setarg.c_str());
      return CpptrajState::ERR;
    }
    mprintf("\tDefault RNG set to '%s'\n", Random_Number::CurrentDefaultRngStr());
  }

  return CpptrajState::OK;
}
