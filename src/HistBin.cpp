#include <cmath> // ceil
#include "HistBin.h"
#include "CpptrajStdio.h"
#include "Constants.h" // SMALL

// CONSTRUCTOR
HistBin::HistBin() : max_(0.0), bins_(-1) {}

HistBin::HistBin(int b, double m, double s, std::string const& l) :
  Dimension(m, s, l),
  bins_(b)
{
  max_ = (Step() * (double)bins_) + Min();
}

// COPY CONSTRUCTOR
HistBin::HistBin(const HistBin& rhs) : Dimension(rhs), max_(rhs.max_), bins_(rhs.bins_) {}

// Assignment
HistBin& HistBin::operator=(const HistBin& rhs) {
  if (this != &rhs) {
    Dimension::operator=(rhs);
    max_ = rhs.max_;
    bins_ = rhs.bins_;
  }
  return *this;
}

// HistBin::CalcBinsOrStep()
/** If both have been defined, use the specified bins and calculate a new 
  * step. If neither have been defined, use default bins and calculate 
  * step. When calculating bins from a stepsize, round up.
  */
int HistBin::CalcBinsOrStep(double minIn, double maxIn, double stepIn, int binsIn,
                            std::string const& labelIn)
{
  if (maxIn - minIn < Constants::SMALL) {
    mprinterr("Error: HistBin: Max (%g) must be greater than min (%g)\n", maxIn, minIn);
    return 1;
  }
  double tempStep = stepIn;
  if (binsIn > 0 && tempStep != 0.0) {
    mprintf("Warning: Both bins (%i) and step (%g) have been specified. Recalculating step.\n",
            binsIn, tempStep);
    tempStep = 0.0;
  }
  if (binsIn < 1 && tempStep == 0.0) {
    mprinterr("Error: [%s] Bins and step undefined.\n", labelIn.c_str());
    return 1;
  }
  bins_ = binsIn;
  max_ = maxIn;
  if (bins_ > 0) {
    //if (debug>0)
    mprintf("\t\tCalculating step from min=%g max=%g bins=%i.\n", minIn, max_, bins_);
    tempStep = max_ - minIn;
    tempStep = tempStep / (double)(bins_);
  } else {
    //if (debug>0)
    mprintf("\t\tCalculating bins from min=%g max=%g step=%g.\n", minIn, max_, tempStep);
    double temp = ((max_ - minIn) / tempStep);
    temp = ceil(temp);
    bins_ = (int)temp;
  }
  SetDimension(minIn, tempStep, labelIn);
  return 0;
}

void HistBin::PrintHistBin() const {
  mprintf("\tDim %s: %f->%f, step %f, %i bins.\n", label(), Min(), max_, Step(), bins_);
}
