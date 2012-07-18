#include "Action_MRT.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_MRT::Action_MRT() :
  time_(1),
  nStar_(0),
  lowerCutoff2_(0),
  upperCutoff2_(0)
{}

// Action_MRT::init()
/** Usage: 
  * mrt out <filename> ([autocorr <filename> [tcorr <time>] [toffset <time>]])
  *    [lower <dist>] [upper <dist>] [time <t>] [tstar <t>] [noimage]
  *    ([solvent <mask> | solute <mask>])
  *    (siteatoms <mask> | onemol <mask> | <sitemask1> ... <sitemaskN>)
  */
int Action_MRT::init() {
  useImage_ = !(actionArgs.hasKey("noimage"));
  std::string filename = actionArgs.GetStringKey("out");
  if (filename.empty()) {
    mprinterr("Error: MRT: No output filename specified 'out <filename>'\n");
    return 1;
  }
  time_ = actionArgs.getKeyDouble("time", 1.0);

  // Time that water can be inside/outside without counting it as
  // having left/entered.
  double tstar = actionArgs.getKeyDouble("tstar", 0.0);
  tstar /= (time_ + 0.5);
  nStar_ = (int)tstar;

  // Lower and upper distance limits
  double lowerCutoff = actionArgs.getKeyDouble("lower", 0.01);
  double upperCutoff = actionArgs.getKeyDouble("upper", 3.50);
  lowerCutoff2_ = lowerCutoff * lowerCutoff;
  upperCutoff2_ = upperCutoff * upperCutoff;

  // If specified, filename for autocorrelation fn
  autoCorr_ = actionArgs.GetStringKey("autocorr");

  return 0;
}
