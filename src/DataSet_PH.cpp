#include "DataSet_PH.h"
#include "CpptrajStdio.h" // DEBUG

DataSet_PH::DataSet_PH() :
  // 0 dim indicates DataSet-specific write
  DataSet(PH, GENERIC, TextFormat(TextFormat::DOUBLE, 10, 4), 0),
  nframes_(0)
{}

// =============================================================================
DataSet_PH::Residue::Residue(NameType const& rname, int resid,
                             Iarray const& protcnts, int max_prots) :
  resname_(rname), resid_(resid), protcnts_(protcnts)
{
  protonated_.reserve(protcnts_.size());
  for (Iarray::const_iterator it = protcnts_.begin(); it != protcnts_.end(); ++it)
    protonated_.push_back(*it == max_prots);
}

DataSet_PH::Residue::Residue(Residue const& rhs) :
  resname_(rhs.resname_),
  resid_(rhs.resid_),
  protcnts_(rhs.protcnts_),
  protonated_(rhs.protonated_),
  states_(rhs.states_)
{}

DataSet_PH::Residue& DataSet_PH::Residue::operator=(Residue const& rhs) {
  if (this != &rhs) {
    resname_ = rhs.resname_;
    resid_ = rhs.resid_;
    protcnts_ = rhs.protcnts_;
    protonated_ = rhs.protonated_;
    states_ = rhs.states_;
  }
  return *this;
}

void DataSet_PH::Residue::Print() const {
  mprintf("\t%s %8i (%zu frames)", *resname_, resid_, states_.size());
  for (unsigned int i = 0; i != protcnts_.size(); i++)
    mprintf(" %2i (%1i),", protcnts_[i], (int)protonated_[i]);
  mprintf("\n");
}
