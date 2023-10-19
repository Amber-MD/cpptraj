#include <algorithm>
#include "InternalCoords.h"
#include "../Topology.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Structure;

const int InternalCoords::NO_ATOM = -1;

/** CONSTRUCTOR */
InternalCoords::InternalCoords() : 
  ati_(NO_ATOM)
{
  std::fill(idx_, idx_+3, NO_ATOM);
  std::fill(val_, val_+3, 0);
}

/** CONSTRUCTOR - Take indices and values */
InternalCoords::InternalCoords(int atI, int atJ, int atK, int atL,
                               double dist, double theta, double phi) :
  ati_(atI)
{
  if (ati_ < 0)
    ati_ = NO_ATOM;
  if (atJ < 0)
    idx_[0] = NO_ATOM;
  else
    idx_[0] = atJ;
  if (atK < 0)
    idx_[1] = NO_ATOM;
  else
    idx_[1] = atK;
  if (atL < 0)
    idx_[2] = NO_ATOM;
  else
    idx_[2] = atL;

  val_[0] = dist;
  val_[1] = theta;
  val_[2] = phi;
}

/** COPY CONSTRUCTOR */
InternalCoords::InternalCoords(InternalCoords const& rhs) :
  ati_(rhs.ati_)
{
  std::copy( rhs.idx_, rhs.idx_+3, idx_ );
  std::copy( rhs.val_, rhs.val_+3, val_ );
}

/** ASSIGNMENT */
InternalCoords& InternalCoords::operator=(InternalCoords const& rhs) {
  if (&rhs == this) return *this;
  ati_ = rhs.ati_;
  std::copy( rhs.idx_, rhs.idx_+3, idx_ );
  std::copy( rhs.val_, rhs.val_+3, val_ );
  return *this;
}

/** Print to stdout */
void InternalCoords::printIC(Topology const& top) const {
  mprintf(" %6i %6i %6i %6i [%s - %s - %s - %s] r=%g t=%g p=%g\n",
          AtI()+1, AtJ()+1, AtK()+1, AtL()+1,
          top.AtomMaskName(AtI()).c_str(),
          top.AtomMaskName(AtJ()).c_str(),
          top.AtomMaskName(AtK()).c_str(),
          top.AtomMaskName(AtL()).c_str(),
          val_[0], val_[1], val_[2]);
}
