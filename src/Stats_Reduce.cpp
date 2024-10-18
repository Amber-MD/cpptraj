#include "Stats_Reduce.h"
#ifdef MPI
/** Calculate overall average/stdev for array elements using the parallel
  * form of the Online algorithm ( Chan, T.F.; Golub, G.H.; LeVeque, R.J.
  * (1979), "Updating Formulae and a Pairwise Algorithm for Computing Sample
  * Variances.", Technical Report STAN-CS-79-773, Department of Computer
  * Science, Stanford University ).
  */
int Cpptraj::Stats_Reduce(Parallel::Comm const& trajComm_,
                          std::vector< Stats<double> >& histogram_,
                          unsigned long& maxbin_)
{
  if (trajComm_.Size() < 2) return 0;
  std::vector<double> buffer;
  unsigned long rank_size;
  if (trajComm_.Master()) {
    for (int rank = 1; rank < trajComm_.Size(); rank++) {
      // 1. Get size of histogram on rank.
      trajComm_.SendMaster(&rank_size, 1, rank, MPI_UNSIGNED_LONG);
      // 2. Receive histogram from rank.
      buffer.resize(3*rank_size, 0.0); // mean, m2, N
      unsigned long master_size = (unsigned long)histogram_.size();
      trajComm_.SendMaster(&buffer[0], buffer.size(), rank, MPI_DOUBLE);
      unsigned long idx = 0; // Index into buffer
      // Only sum for bins where master and rank both have data
      unsigned long Nbins = std::min( master_size, rank_size );
      for (unsigned long i = 0; i < Nbins; i++, idx += 3) {
        double mB = buffer[idx  ];
        double sB = buffer[idx+1];
        double nB = buffer[idx+2];
        histogram_[i].Combine( Stats<double>(nB, mB, sB) );
      }
      // If rank had more data than master, fill in data
      if (rank_size > master_size) {
        histogram_.resize( rank_size );
        maxbin_ = (unsigned long)histogram_.size() - 1;
        idx = master_size * 3;
        for (unsigned long i = master_size; i < rank_size; i++, idx += 3)
          histogram_[i] = Stats<double>( buffer[idx+2], buffer[idx], buffer[idx+1] );
      }
    }
  } else {
    // 1. Send size of histogram on this rank to master.
    rank_size = (unsigned long)histogram_.size();
    trajComm_.SendMaster(&rank_size, 1, trajComm_.Rank(), MPI_UNSIGNED_LONG);
    // 2. Place histogram data into a buffer and send to master.
    buffer.reserve(3*histogram_.size()); // mean, m2, N
    for (unsigned long i = 0; i < histogram_.size(); i++) {
      buffer.push_back( (double)histogram_[i].mean()  );
      buffer.push_back( (double)histogram_[i].M2()    );
      buffer.push_back( (double)histogram_[i].nData() );
    }
    trajComm_.SendMaster(&buffer[0], buffer.size(), trajComm_.Rank(), MPI_DOUBLE);
  }
  
  return 0;
}
#endif /* MPI */
