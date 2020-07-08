#include "PotentialTerm_InitOpts.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

PotentialTerm::InitOpts::InitOpts() :
  shakeType_(Constraints::OFF),
  scaleEE_(1.0/1.2), // Amber default
  scaleNB_(1.0/2.0), // Amber default
  cutEE_(8.0),       // in Ang., Amber default
  cutNB_(8.0)        // in Ang., Amber default
{}

PotentialTerm::InitOpts::InitOpts(ArgList& argIn)
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
    else
      mprinterr("Error: Unrecognized shake arg: %s\n", shakearg.c_str());
  }
}
