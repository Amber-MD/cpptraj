#include "EwaldOptions.h"
#include "CpptrajStdio.h"
#include "ArgList.h"

/** CONSTRUCTOR */
EwaldOptions::EwaldOptions() :
  allowLjPme_(true),
  type_(NOT_SET),
  cutoff_(8),
  dsumtol_(1E-5),
  ewcoeff_(0),
  erfcDx_(0),
  skinnb_(2),
  lwcoeff_(-1),
  ljswidth_(0),
  rsumtol_(5E-5),
  maxexp_(0),
  mlimits1_(0),
  mlimits2_(0),
  mlimits3_(0),
  npoints_(6),
  nfft1_(-1),
  nfft2_(-1),
  nfft3_(-1)
{}

/** Set whether we want to read LJ PME options. */
void EwaldOptions::AllowLjPme(bool allow) {
  allowLjPme_ = allow;
}

// Keywords lines for help text
const char* EwaldOptions::KeywordsCommon1() {
  return "[cut <cutoff>] [dsumtol <dtol>] [ewcoeff <coeff>]";
}

const char* EwaldOptions::KeywordsCommon2() {
  return "[erfcdx <dx>] [skinnb <skinnb>] [ljswidth <width>]";
}

const char* EwaldOptions::KeywordsRegEwald() {
  return "[rsumtol <rtol>] [maxexp <max>] [mlimits <X>,<Y>,<Z>]";
}

const char* EwaldOptions::KeywordsPME() {
  return "[order <order>] [nfft <nfft1>,<nfft2>,<nfft3>]";
}

const char* EwaldOptions::KeywordsLjpme() {
  return "[ljpme [ewcoefflj <ljcoeff>]]";
}

/** Get a comma-separated list of 3 integer values after a keyword. */
int EwaldOptions::GetCommaSeparatedArgs(ArgList& argIn, const char* keyword, int& i1, int& i2, int& i3, int defaultVal) {
  std::string marg = argIn.GetStringKey(keyword);
  if (!marg.empty()) {
    ArgList mlim(marg, ",");
    if (mlim.Nargs() != 3) {
      mprinterr("Error: Need 3 integers in comma-separated list for '%s'\n", keyword);
      return 1;
    }
    i1 = mlim.getNextInteger(0);
    i2 = mlim.getNextInteger(0);
    i3 = mlim.getNextInteger(0);
  } else {
    i1 = defaultVal;
    i2 = defaultVal;
    i3 = defaultVal;
  }
  return 0;
}

/** Parse Ewald options from ArgList. */
int EwaldOptions::GetOptions(OptType typeIn, ArgList& actionArgs, const char* desc) {
  type_ = typeIn;
  // Common options
  cutoff_ = actionArgs.getKeyDouble("cut", 8.0);
  dsumtol_ = actionArgs.getKeyDouble("dsumtol", 1E-5);
  ewcoeff_ = actionArgs.getKeyDouble("ewcoeff", 0.0);
  erfcDx_ = actionArgs.getKeyDouble("erfcdx", 0.0);
  skinnb_ = actionArgs.getKeyDouble("skinnb", 2.0);
  // LJ Options
  // NOTE: lwcoeff_ > 0 is LJPME on. An lwcoeff_ of -1 is off, and 0 is set from ewcoeff_.
  lwcoeff_ = -1;
  if (actionArgs.hasKey("ljpme")) {
    if (!allowLjPme_) {
      mprinterr("Error: LJ PME option 'ljpme' not allowed for '%s'\n", desc);
      return 1;
    }
    lwcoeff_ = 0.4;
  }
  lwcoeff_ = actionArgs.getKeyDouble("ewcoefflj", lwcoeff_);
  if (!allowLjPme_ && lwcoeff_ >= 0) {
    mprinterr("Error: LJ PME option 'ewcoefflj' not allowed for '%s'\n", desc);
    return 1;
  }
  ljswidth_ = actionArgs.getKeyDouble("ljswidth", 0.0);
  // Regular Ewald options
  if (type_ != PME) {
    rsumtol_ = actionArgs.getKeyDouble("rsumtol", 5E-5);
    maxexp_ = actionArgs.getKeyDouble("maxexp", 0.0);
    if (GetCommaSeparatedArgs(actionArgs, "mlimits", mlimits1_, mlimits2_, mlimits3_, 0)) return 1;
  }
  // PME options
  if (type_ != REG_EWALD) {
    npoints_ = actionArgs.getKeyInt("order", 6);
    if (GetCommaSeparatedArgs(actionArgs, "nfft", nfft1_, nfft2_, nfft3_, -1)) return 1;
  }
  return 0;
}

/** Print current options to STDOUT. */
void EwaldOptions::PrintOptions() const {
  // Common options
  mprintf("\tDirect space cutoff= %.4f\n", cutoff_);
  if (dsumtol_ != 0.0)
    mprintf("\tDirect sum tolerance= %g\n", dsumtol_);
  if (ewcoeff_ == 0.0)
    mprintf("\tWill determine Ewald coefficient from cutoff and direct sum tolerance.\n");
  else
    mprintf("\tEwald coefficient= %.4f\n", ewcoeff_);
  if (erfcDx_ > 0.0)
    mprintf("\tERFC table dx= %g\n", erfcDx_);
  if (skinnb_ > 0)
    mprintf("\tSize of non-bonded \"skin\"= %.4f\n", skinnb_);
  // Regular Ewald options
  if (type_ != PME) {
    if (rsumtol_ != 0.0)
      mprintf("\tReciprocal sum tolerance= %g\n", rsumtol_);
    if (maxexp_ == 0.0)
      mprintf("\tWill determine MaxExp from Ewald coefficient and direct sum tolerance.\n");
    else
      mprintf("\tMaxExp= %g\n", maxexp_);
    if (mlimits1_ < 1 && mlimits2_ < 1 && mlimits3_ < 1)
      mprintf("\tWill determine number of reciprocal vectors from MaxExp.\n");
    else
      mprintf("\tNumber of reciprocal vectors in each direction= {%i,%i,%i}\n",
              mlimits1_, mlimits2_, mlimits3_);
  }
  // PME options
  if (type_ != REG_EWALD) {
    mprintf("\tSpline order= %i\n", npoints_);
    if (nfft1_ < 1 && nfft2_ < 1 && nfft3_ < 1)
      mprintf("\tWill determine number of FFT grid points from box size.\n");
    else
      mprintf("\tNumber of FFT grid points in each direction= {%i,%i,%i}\n",
              nfft1_, nfft2_, nfft3_);
  }
  // LJ options
  if (type_ != REG_EWALD) {
    if (lwcoeff_ < 0)
      mprintf("\tUsing long range correction for nonbond VDW calc.\n");
    else if (lwcoeff_ > 0.0)
      mprintf("\tUsing Lennard-Jones PME with Ewald coefficient %.4f\n", lwcoeff_);
    else
      mprintf("\tLennard-Jones PME Ewald coefficient will be set to elec. Ewald coefficient.\n");
  }
  if (ljswidth_ > 0.0)
    mprintf("\tWidth of LJ switch region: %.4f Ang.\n", ljswidth_);
}
