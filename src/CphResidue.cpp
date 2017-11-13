#include "CphResidue.h"
//#inc lude "CpptrajStdio.h"

CphResidue::CphResidue(NameType const& rname, int resid,
                             Iarray const& protcnts, int max_prots) :
  resname_(rname), resid_(resid), protcnts_(protcnts)
{
  protonated_.reserve(protcnts_.size());
  for (Iarray::const_iterator it = protcnts_.begin(); it != protcnts_.end(); ++it)
    protonated_.push_back(*it == max_prots);
}

CphResidue::CphResidue(CphResidue const& rhs) :
  resname_(rhs.resname_),
  resid_(rhs.resid_),
  protcnts_(rhs.protcnts_),
  protonated_(rhs.protonated_)
{}

CphResidue& CphResidue::operator=(CphResidue const& rhs) {
  if (this != &rhs) {
    resname_ = rhs.resname_;
    resid_ = rhs.resid_;
    protcnts_ = rhs.protcnts_;
    protonated_ = rhs.protonated_;
  }
  return *this;
}

/*void CphResidue::Print() const {
  mprintf("\t%s %8i (%zu frames)", *resname_, resid_, states_.size());
  for (unsigned int i = 0; i != protcnts_.size(); i++)
    mprintf(" %2i (%1i),", protcnts_[i], (int)protonated_[i]);
  mprintf("\n");
}*/
