#include <cmath> // ceil
#include "Dimension.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Dimension::Dimension() :
  min_(0),
  max_(0),
  step_(-1),
  bins_(-1),
  offset_(0),
  minIsSet_(false),
  maxIsSet_(false)
{}

// COPY CONSTRUCTOR
Dimension::Dimension(const Dimension& rhs) :
  label_(rhs.label_),
  min_(rhs.min_),
  max_(rhs.max_),
  step_(rhs.step_),
  bins_(rhs.bins_),
  offset_(rhs.offset_),
  minIsSet_(rhs.minIsSet_),
  maxIsSet_(rhs.maxIsSet_)
{}

// Assignment
Dimension& Dimension::operator=(const Dimension& rhs) {
  if (this == &rhs) return *this;
  label_ = rhs.label_;
  min_ = rhs.min_;
  max_ = rhs.max_;
  step_ = rhs.step_;
  bins_ = rhs.bins_;
  offset_ = rhs.offset_;
  minIsSet_ = rhs.minIsSet_;
  maxIsSet_ = rhs.maxIsSet_;
  return *this;
}

// Dimension::CalcBinsOrStep()
/** If both have been defined, use the specified bins and calculate a new 
  * step. If neither have been defined, use default bins and calculate 
  * step. When calculating bins from a stepsize, round up.
  */
int Dimension::CalcBinsOrStep() {
  if (!minIsSet_ || !maxIsSet_) {
    mprinterr("Error: Dimension::CalcBinsOrStep: [%s] min or max not set.\n", label_.c_str());
    return 1;
  }
  if (bins_!=-1 && step_!=-1) {
    mprintf("\tHist: Bins and step have been specified. Recalculating step.\n");
    step_ = -1;
  }

  if (bins_==-1 && step_==-1) {
    mprinterr("Error: Hist: [%s] Bins and step undefined.\n", label_.c_str());
    return 1;
  }

  if (step_==-1) {
    //if (debug>0) mprintf("\t\tCalculating step.\n");
    if (bins_ <= 0) {
      mprinterr("Error: Hist: Dimension %s: bins <=0!\n", label_.c_str());
      return 1;
    }
    step_ = max_ - min_;
    step_ = step_ / bins_;
  } else if (bins_ == -1) {
    //if (debug>0) mprintf("\t\tCalculating bins.\n");
    if (step_ <= 0) {
      mprinterr("Error: Hist: Dimension %s: step <=0!\n",label_.c_str());
      return 1;
    } 
    double temp = ((max_ - min_) / step_);
    temp = ceil(temp);
    bins_ = (int)temp;
  }
  return 0;
}

void Dimension::PrintDim() const {
  mprintf("\tDim %s: %f->%f, step %f, %i bins.\n", label_.c_str(),
          min_, max_, step_, bins_);
}
