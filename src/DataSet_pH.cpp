#include "DataSet_pH.h"
#include "CpptrajStdio.h"
//#incl ude "StringRoutines.h" // DEBUG

DataSet_pH::DataSet_pH() :
  DataSet_1D(PH, TextFormat(TextFormat::INTEGER, 10)),
  t0_(-1.0),
  dt_(-10.0),
  mc_stepsize_(-1)
{}

int DataSet_pH::Allocate(SizeArray const& sizeIn) {
  if (!sizeIn.empty()) {
    states_.reserve( sizeIn[0] );
    recType_.reserve( sizeIn[0] );
  }
  return 0;
}

void DataSet_pH::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= states_.size())
    cbuffer.Printf(format_.fmt(), 0);
  else
    cbuffer.Printf(format_.fmt(), states_[pIn[0]]);
}

void DataSet_pH::Resize(size_t n) {
  states_.resize( n );
  recType_.resize( n );
}

void DataSet_pH::Info() const {
  mprintf(" (%s %i pH= %.2f)", *(res_.Name()), res_.Num(), solvent_pH_);
}
# ifdef MPI
void DataSet_pH::Consolidate(Parallel::Comm const& commIn, int rank)
{
  if (commIn.Rank() == rank) {
    /*commIn.Barrier(); // DEBUG
    std::string msg("MASTER:");
    for (std::vector<int>::const_iterator it = states_.begin(); it != states_.end(); ++it)
      msg.append(" " + integerToString(*it));
    rprintf("%s\n", msg.c_str());
    mflush();
    commIn.Barrier();
    commIn.Barrier(); // END DEBUG*/
    std::vector<int> itmp = states_;
    commIn.Reduce(rank, &states_[0], &itmp[0], states_.size(), MPI_INT, MPI_SUM);
    itmp = recType_;
    commIn.Reduce(rank, &recType_[0], &itmp[0], recType_.size(), MPI_INT, MPI_SUM);
  } else {
    /*commIn.Barrier(); // DEBUG
    commIn.Barrier();
    std::string msg("CHILD_:");
    for (std::vector<int>::const_iterator it = states_.begin(); it != states_.end(); ++it)
      msg.append(" " + integerToString(*it));
    rprintf("%s\n", msg.c_str());
    commIn.Barrier();*/
    commIn.Reduce(rank, 0, &states_[0], states_.size(), MPI_INT, MPI_SUM); 
    commIn.Reduce(rank, 0, &recType_[0], full.size(), MPI_INT, MPI_SUM);
  }
}
#endif
