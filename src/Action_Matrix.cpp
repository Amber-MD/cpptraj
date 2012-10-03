#include "Action_Matrix.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Matrix::Action_Matrix() :
  Mat_(0),
  outtype_(BYATOM),
  start_(0),
  stop_(-1),
  offset_(1),
  order_(2)
{}

