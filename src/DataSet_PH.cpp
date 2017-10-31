#include "DataSet_PH.h"
#include "CpptrajStdio.h"

DataSet_PH::DataSet_PH() :
  // 0 dim indicates DataSet-specific write
  DataSet(PH, GENERIC, TextFormat(TextFormat::DOUBLE, 10, 4), 0)
{}

int DataSet_PH::Allocate(SizeArray const& sizeIn) {
  if (!sizeIn.empty()) {
    for (Rarray::iterator res = residues_.begin(); res != residues_.end(); ++res)
      res->Allocate( sizeIn[0] );
    solvent_pH_.reserve( sizeIn[0] );
  }
  return 0;
}

void DataSet_PH::Resize(size_t n) {
  for (Rarray::iterator res = residues_.begin(); res != residues_.end(); ++res)
    res->Resize( n );
  solvent_pH_.resize(n, 0.0);
}

void DataSet_PH::Info() const {
  mprintf(" (%zu residues", residues_.size());
  if (!residues_.empty())
    mprintf(", %zu frames", residues_[0].Nframes());
  mprintf(")");
}

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
