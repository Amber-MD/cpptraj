#include "MdOpts.h"
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "Constants.h"

MdOpts::MdOpts() :
  shakeType_(Constraints::OFF),
  scaleEE_(1.0/1.2), // Amber default
  scaleNB_(1.0/2.0), // Amber default
  cutEE_(8.0),       // in Ang., Amber default
  cutNB_(8.0),       // in Ang., Amber default
  qfac_(Constants::ELECTOAMBER * Constants::ELECTOAMBER), // Assume charges in elec. units, Amber default
  nExclude_(4)       // Exclude dihedral, angle, bond
{}

/** Display options recognized by MdOpts to stdout. */
void MdOpts::PrintHelp() {
  mprintf("\tcut <cutoff>     : Set nonbonded interaction cutoff in Ang. (electrostatics and vdW).\n"
          "\tcutee <cutoff>   : Set electrostatics interaction cutoff in Ang.\n"
          "\tcutnb <cutoff>   : Set vdW interaction cutoff in Ang.\n"
          "\tscaleee <factor> : Scaling factor to multiply 1-4 electrostatics interactions by.\n"
          "\tscalenb <factor> : Scaling factor to multiply 1-4 vdW interactions by.\n"
          "\tshake <type>     : Use SHAKE constraints of <type>;\n"
          "\t                   <type> = 'hydrogen' : Constrain bonds to hydrogen.\n"
          "\t                   <type> = 'all'      : Constrain all bonds.\n"
          "\tnexclude <#>     : Number of bonded atoms within which nonbonded interactions are excluded.\n"
          "\tqfac <factor>    : Factor to use in electrostatic calculation.\n"
  );
}

/** Get options from incoming argument list. */
int MdOpts::GetOptsFromArgs(ArgList& argIn)
{
  cutEE_ = argIn.getKeyDouble("cutee", 8.0);
  cutNB_ = argIn.getKeyDouble("cutnb", 8.0);
  if (argIn.Contains("cut")) {
    double cut = argIn.getKeyDouble("cut", 8.0);
    cutEE_ = cut;
    cutNB_ = cut;
  }
  scaleEE_ = argIn.getKeyDouble("scaleee", 1.0/1.2);
  scaleNB_ = argIn.getKeyDouble("scalenb", 1.0/2.0);
  shakeType_ = Constraints::OFF;
  std::string shakearg = argIn.GetStringKey("shake");
  if (!shakearg.empty()) {
    if (shakearg == "h" || shakearg == "hydrogen")
      shakeType_ = Constraints::BONDS_TO_H;
    else if (shakearg == "all")
      shakeType_ = Constraints::ALL_BONDS;
    else {
      mprinterr("Error: Unrecognized shake arg: %s\n", shakearg.c_str());
      return 1;
    }
  }
  nExclude_ = argIn.getKeyInt("nexclude", 4);
  qfac_ = argIn.getKeyDouble("qfac", Constants::ELECTOAMBER * Constants::ELECTOAMBER);
  
  return 0;
}

void MdOpts::PrintOpts() const {
  mprintf("\tElectrostatics cutoff : %g Ang.\n", cutEE_);
  mprintf("\tvdW cutoff            : %g Ang.\n", cutNB_);
  mprintf("\t1-4 Elec. scaling     : %g\n", scaleEE_);
  mprintf("\t1-4 vdW scaling       : %g\n", scaleNB_);
  mprintf("\tSHAKE constraints     : %s\n", Constraints::shakeString(shakeType_));
  mprintf("\tCoulomb factor        : %g\n", qfac_);
  mprintf("\tExclude beyond        : %i atoms\n", nExclude_);
}
