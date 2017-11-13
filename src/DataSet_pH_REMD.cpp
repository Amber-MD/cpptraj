#include "DataSet_pH_REMD.h"
#include "CpptrajStdio.h"

DataSet_pH_REMD::DataSet_pH_REMD() :
  // 0 dim indicates DataSet-specific write
  DataSet(PH, GENERIC, TextFormat(TextFormat::DOUBLE, 10, 4), 0)
{}

int DataSet_pH_REMD::Allocate(SizeArray const& sizeIn) {
  if (!sizeIn.empty()) {
    for (Rarray::iterator res = residues_.begin(); res != residues_.end(); ++res)
      res->Allocate( sizeIn[0] );
    solvent_pH_.reserve( sizeIn[0] );
  }
  return 0;
}

void DataSet_pH_REMD::Resize(size_t n) {
  for (Rarray::iterator res = residues_.begin(); res != residues_.end(); ++res)
    res->Resize( n );
  solvent_pH_.resize(n, 0.0);
}

void DataSet_pH_REMD::Info() const {
  mprintf(" (%zu residues", residues_.size());
  if (!residues_.empty())
    mprintf(", %zu frames", residues_[0].Nframes());
  mprintf(")");
}
# ifdef MPI
void DataSet_pH_REMD::Consolidate(Parallel::Comm const& commIn, int rank)
{
  if (commIn.Rank() == rank) {
    //commIn.Barrier(); // DEBUG
    //rprintf("MASTER:\n");
    //for (std::vector<int>::const_iterator it = residues_.front().begin(); it != residues_.front().end(); ++it)
    //  rprintf("%4u %1i\n", it-residues_.front().begin()+1, *it);
    //mflush();
    //commIn.Barrier(); // DEBUG
    std::vector<float> ftmp = solvent_pH_;
    commIn.Reduce(rank, &solvent_pH_[0], &ftmp[0], solvent_pH_.size(), MPI_FLOAT, MPI_SUM);
    ftmp.clear();
    if (!residues_.empty()) {
      std::vector<int> itmp;
      itmp.reserve( Nframes() );
      for (Rarray::iterator res = residues_.begin(); res != residues_.end(); ++res) {
        itmp = res->States();
        commIn.Reduce(rank, res->StatesPtr(), &itmp[0], res->States().size(), MPI_INT, MPI_SUM);
      }
    }
  } else {
    //commIn.Barrier(); // DEBUG
    //commIn.Barrier();
    //rprintf("CHILD_:\n");
    //for (std::vector<int>::const_iterator it = residues_.front().begin(); it != residues_.front().end(); ++it)
    //  rprintf("%4u %1i\n", it-residues_.front().begin()+1, *it);
    commIn.Reduce(rank, 0, &solvent_pH_[0], solvent_pH_.size(), MPI_FLOAT, MPI_SUM);
    if (!residues_.empty()) {
      for (Rarray::iterator res = residues_.begin(); res != residues_.end(); ++res)
        commIn.Reduce(rank, 0, res->StatesPtr(), res->States().size(), MPI_INT, MPI_SUM);
    }
  }
}
#endif

// =============================================================================
DataSet_pH_REMD::Residue::Residue(NameType const& rname, int resid,
                             Iarray const& protcnts, int max_prots) :
  resname_(rname), resid_(resid), protcnts_(protcnts)
{
  protonated_.reserve(protcnts_.size());
  for (Iarray::const_iterator it = protcnts_.begin(); it != protcnts_.end(); ++it)
    protonated_.push_back(*it == max_prots);
}

DataSet_pH_REMD::Residue::Residue(Residue const& rhs) :
  resname_(rhs.resname_),
  resid_(rhs.resid_),
  protcnts_(rhs.protcnts_),
  protonated_(rhs.protonated_),
  states_(rhs.states_)
{}

DataSet_pH_REMD::Residue& DataSet_pH_REMD::Residue::operator=(Residue const& rhs) {
  if (this != &rhs) {
    resname_ = rhs.resname_;
    resid_ = rhs.resid_;
    protcnts_ = rhs.protcnts_;
    protonated_ = rhs.protonated_;
    states_ = rhs.states_;
  }
  return *this;
}

void DataSet_pH_REMD::Residue::Print() const {
  mprintf("\t%s %8i (%zu frames)", *resname_, resid_, states_.size());
  for (unsigned int i = 0; i != protcnts_.size(); i++)
    mprintf(" %2i (%1i),", protcnts_[i], (int)protonated_[i]);
  mprintf("\n");
}
