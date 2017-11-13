#include "DataSet_pH.h"
#include "CpptrajStdio.h"

DataSet_pH::DataSet_pH() :
  DataSet_1D(PH, TextFormat(TextFormat::DOUBLE, 10, 4))
{}

int DataSet_pH::Allocate(SizeArray const& sizeIn) {
  if (!sizeIn.empty()) {
    states_.reserve( sizeIn[0] );
    //solvent_pH_.reserve( sizeIn[0] );
  }
  return 0;
}

void DataSet_pH::Resize(size_t n) {
  states_.resize( n );
  //solvent_pH_.resize(n, 0.0);
}

void DataSet_pH::Info() const {
  mprintf(" (%s %i pH= %.2f", *resname_, resid_, solvent_pH_);
  if (!pH_Values_.empty())
    mprintf(", unsorted REMD");
  mprintf(")");
}
# ifdef MPI
void DataSet_pH::Consolidate(Parallel::Comm const& commIn, int rank)
{
  if (commIn.Rank() == rank) {
    //commIn.Barrier(); // DEBUG
    //rprintf("MASTER:\n");
    //for (std::vector<int>::const_iterator it = residues_.front().begin(); it != residues_.front().end(); ++it)
    //  rprintf("%4u %1i\n", it-residues_.front().begin()+1, *it);
    //mflush();
    //commIn.Barrier(); // DEBUG
    std::vector<float> ftmp = pH_Values_;
    commIn.Reduce(rank, &pH_Values_[0], &ftmp[0], pH_Values_.size(), MPI_FLOAT, MPI_SUM);
    ftmp.clear();
    std::vector<int> itmp = states_;
    commIn.Reduce(rank, &states_[0],    &itmp[0], states_.size(),    MPI_INT,   MPI_SUM);
  } else {
    //commIn.Barrier(); // DEBUG
    //commIn.Barrier();
    //rprintf("CHILD_:\n");
    //for (std::vector<int>::const_iterator it = residues_.front().begin(); it != residues_.front().end(); ++it)
    //  rprintf("%4u %1i\n", it-residues_.front().begin()+1, *it);
    commIn.Reduce(rank, 0, &pH_Values_[0], pH_Values_.size(), MPI_FLOAT, MPI_SUM);
    commIn.Reduce(rank, 0, &states_[0],    states_.size(),    MPI_INT,   MPI_SUM);
  }
}
#endif

// DataSet_pH::SetResidueInfo()
void DataSet_pH::SetResidueInfo(NameType const& rnameIn, int residIn,
                                Iarray const& protcntsIn, int max_prots)
{
  resname_ = rnameIn;
  resid_ = residIn;
  protcnts_ = protcntsIn;
  protonated_.reserve(protcnts_.size());
  for (Iarray::const_iterator it = protcnts_.begin(); it != protcnts_.end(); ++it)
    protonated_.push_back(*it == max_prots);
}

void DataSet_pH::SetResidueInfo(DataSet_pH const& rhs) {
  resname_ = rhs.resname_;
  resid_ = rhs.resid_;
  protcnts_ = rhs.protcnts_;
  protonated_ = rhs.protonated_;
}

/* DataSet_pH::Residue& DataSet_pH::Residue::operator=(Residue const& rhs) {
  if (this != &rhs) {
    resname_ = rhs.resname_;
    resid_ = rhs.resid_;
    protcnts_ = rhs.protcnts_;
    protonated_ = rhs.protonated_;
    states_ = rhs.states_;
  }
  return *this;
}*/
