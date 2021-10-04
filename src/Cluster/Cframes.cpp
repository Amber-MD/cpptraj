#include <algorithm>
#include "Cframes.h"

bool Cpptraj::Cluster::Cframes::HasFrame(int frame) const {
  Iarray::const_iterator it = std::find(frames_.begin(), frames_.end(), frame);
  return !(it == frames_.end());
}

int Cpptraj::Cluster::Cframes::FrameIdx(int frame) const {
  Iarray::const_iterator it = std::find(frames_.begin(), frames_.end(), frame);
  if (it == frames_.end()) return -1;
  return (int)(it - frames_.begin());
}

void Cpptraj::Cluster::Cframes::Remove(int frame) {
  Iarray::iterator pend = std::remove( frames_.begin(), frames_.end(), frame);
  std::size_t newsize = pend - frames_.begin();
  frames_.resize( newsize );
}

void Cpptraj::Cluster::Cframes::Sort() {
  std::sort(frames_.begin(), frames_.end());
}
