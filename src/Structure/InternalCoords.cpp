#include <algorithm>
#include "InternalCoords.h"

using namespace Cpptraj::Structure;

const int InternalCoords::NO_ATOM = -1;

/** CONSTRUCTOR */
InternalCoords::InternalCoords() {
  std::fill(idx_, idx_+3, NO_ATOM);
  std::fill(val_, val_+3, 0);
}
