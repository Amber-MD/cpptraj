#include <cmath>

#include "Action_PairDist.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "DistRoutines.h"
#include "DataSet_Mesh.h"


/** Calculate pair distribution function P(r) between two masks.
  * \author Hannes H. Loeffler.
  */

// CONSTRUCTOR
Action_PairDist::Action_PairDist() :
  delta_(0.0),
  maxbin_(0),
  single_mask_(false),
  nframes_(0)
{}

void Action_PairDist::Help() const {
  mprintf("\t[out <filename>] mask <mask> [mask2 <mask>] [delta <resolution>]\n"
          "\t[maxdist <distance>] [noimage]\n"
          "  Calculate pair distribution function P(r) within a single mask\n"
          "  or between two masks.\n"
          "  If 'maxdist' is specified the initial histogram max size will be set to\n"
          "  <distance>; if larger distances are encountered the histogram will be\n"
          "  resized appropriately.\n");
}

// Action_PairDist::init()
Action::RetType Action_PairDist::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  nframes_ = 0;
  imageOpt_.InitImaging( !actionArgs.hasKey("noimage") );
  delta_ = actionArgs.getKeyDouble("delta", 0.01);
  double maxDist = actionArgs.getKeyDouble("maxdist", -1.0);
  if (maxDist > 0.0) {
    maxbin_ = (unsigned long)(maxDist / delta_);
    histogram_.resize( maxbin_ + 1 );
  }

  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Get first mask
  std::string mask1 = actionArgs.GetStringKey("mask");
  if (mask1.empty()) {
    mprinterr("Error: pairdist: No mask1 specified.\n");
    return Action::ERR;
  }
  if (mask1_.SetMaskString(mask1)) return Action::ERR;
  // Get second mask if specified
  std::string mask2 = actionArgs.GetStringKey("mask2");
  if (mask2.empty()) {
    single_mask_ = true;
  } else {
    if (mask2_.SetMaskString(mask2)) return Action::ERR;
    single_mask_ = false;
  }
  // Set up data sets
  std::string dsname = actionArgs.GetStringNext();
  if (dsname.empty())
    dsname = init.DSL().GenerateDefaultName("PDIST");
  MetaData md(dsname, "Pr", MetaData::NOT_TS );
  Pr_ = init.DSL().AddSet(DataSet::XYMESH, md);
  md.SetAspect("std");
  std_ = init.DSL().AddSet(DataSet::XYMESH, md);
  if (Pr_ == 0 || std_ == 0) return Action::ERR;
  if (outfile != 0) {
    outfile->AddDataSet( Pr_ );
    outfile->AddDataSet( std_ );
  }
# ifdef MPI
  // Do not need to be synced since not time series
  Pr_->SetNeedsSync( false );
  std_->SetNeedsSync( false );
# endif

  mprintf("    PAIRDIST: Calculate P(r)");
  if (!single_mask_)
    mprintf(" between atoms selected by '%s' and '%s'\n", mask1_.MaskString(), mask2_.MaskString());
  else
    mprintf(" for atoms selected by '%s'\n", mask1_.MaskString());
  if (outfile != 0)
    mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());
  mprintf("\tResolution is %f Ang.\n", delta_);
  if (maxDist > 0.0)
    mprintf("\tInitial histogram max distance= %f Ang\n", maxDist);
  if (imageOpt_.UseImage())
    mprintf("\tImaging enabled if box info present.\n");
  else
    mprintf("\tNo imaging.\n");
# ifdef MPI
  if (trajComm_.Size() > 1 && maxDist < 0.0)
    mprintf("Warning: Due to the way 'pairdist' currently accumulates data, the\n"
            "Warning:   final histogram may have numerical inaccuracies in the\n"
            "Warning:   tail in parallel.\n"
            "Warning: To avoid this try setting the histogram initial max\n"
            "Warning:   distance with the 'maxdist' keyword.\n");
# endif
  return Action::OK;
}

// Action_PairDist::Setup()
Action::RetType Action_PairDist::Setup(ActionSetup& setup) {
  // Set up first mask
  if (setup.Top().SetupIntegerMask(mask1_) ) return Action::ERR;
  mprintf("\t");
  mask1_.BriefMaskInfo();
  if (mask1_.None()) {
    mprintf("Warning: Mask has no atoms.\n");
    return Action::SKIP;
  }
  // Set up second mask if specified
  if (!single_mask_) {
    if (setup.Top().SetupIntegerMask(mask2_) ) return Action::ERR;
    mask2_.BriefMaskInfo();
    mprintf("\n");
    if (mask2_.None()) {
      mprintf("Warning: Second mask has no atoms.\n");
      return Action::SKIP;
    }
    int numInCommon = mask1_.NumAtomsInCommon( mask2_ );
    if (numInCommon > 0) {
      mprinterr("Error: Masks must be non-overlapping (currently %i atoms in common.)\n",
                numInCommon);
      return Action::ERR;
    }
  }
  // See if imaging is possible
  imageOpt_.SetupImaging(setup.CoordInfo().TrajBox().HasBox() );
  if (imageOpt_.ImagingEnabled())
    mprintf("\tDistances will be imaged.\n");

  return Action::OK;
}

/** Do histogram binning */
void Action_PairDist::BinHist( std::vector<double>& tmp,
                               const double* XYZ1, const double* XYZ2, Box const& box)
{
  double Dist = DIST(imageOpt_.ImagingType(), XYZ1, XYZ2, box);

  unsigned long bin = (unsigned long) (Dist / delta_);
  // Resize on the fly if necessary
  if (bin > maxbin_) {
    maxbin_ = bin;
    tmp.resize(maxbin_ + 1);
    histogram_.resize(maxbin_ + 1);
  }
  //if (bin == 326) rprintf("DEBUG: Bin 326 frame %i idx1 %i idx2 %i %f\n", frame_, idx1_, idx2_, Dist);
  tmp[bin]++;
}

// Action_PairDist::DoAction()
Action::RetType Action_PairDist::DoAction(int frameNum, ActionFrame& frm) {
  std::vector<double> tmp(histogram_.size(), 0); // per frame histogram

  if (imageOpt_.ImagingEnabled())
    imageOpt_.SetImageType( frm.Frm().BoxCrd().Is_X_Aligned_Ortho() );

  if (single_mask_) {
    // Only 1 mask
    for (int idx1 = 0; idx1 != mask1_.Nselected(); idx1++) {
      for (int idx2 = idx1 + 1; idx2 != mask1_.Nselected(); idx2++) {
        //frame_ = frameNum; // DEBUG
        //idx1_ = idx1; // DEBUG
        //idx2_ = idx2; // DEBUG
        BinHist(tmp, frm.Frm().XYZ(mask1_[idx1]), frm.Frm().XYZ(mask1_[idx2]), frm.Frm().BoxCrd());
      }
    }
  } else {
    // 2 non-overlapping masks
    for (int idx1 = 0; idx1 != mask1_.Nselected(); idx1++) {
      for (int idx2 = 0; idx2 != mask2_.Nselected(); idx2++) {
        BinHist(tmp, frm.Frm().XYZ(mask1_[idx1]), frm.Frm().XYZ(mask2_[idx2]), frm.Frm().BoxCrd());
      }
    }
  }

  // Update the overall histogram 
  for (unsigned long i = 0; i < tmp.size(); i++)
    histogram_[i].accumulate(tmp[i]);

  nframes_++;

  return Action::OK;
}

/** Since the histogram can be resized on-the-fly, there may be bins in
  * histogram_ not present at the start of the calculation. These bins
  * need to have their n_ field updated to the actual number of frames
  * so it is as if they were present from the start.
  */
void Action_PairDist::UpdateHistogramFrames() {
  for (unsigned long i = 0; i < histogram_.size(); i++)
  {
    unsigned int ndata = (unsigned int)histogram_[i].nData();
    if (ndata != nframes_) {
      //rprintf("DEBUG: Hist bin %lu frames does not match (%u vs %u)\n", i, ndata, nframes_);
      for (unsigned int idx = ndata; idx < nframes_; idx++)
        histogram_[i].accumulate( 0 );
    }
  }
}

#ifdef MPI
/** Calculate overall average/stdev for histogram bins using the parallel
  * form of the Online algorithm ( Chan, T.F.; Golub, G.H.; LeVeque, R.J.
  * (1979), "Updating Formulae and a Pairwise Algorithm for Computing Sample
  * Variances.", Technical Report STAN-CS-79-773, Department of Computer
  * Science, Stanford University ).
  * NOTE: Due to the way pairwise currently works (histogram resized on the
  *       fly) this will not work correctly for a bin that is only sometimes
  *       present as the final bin counts will be different. The only way I
  *       could see around this is to have all the threads sync up any time
  *       a histogram is resized, which seems pretty bad performance-wise.
  *       Adding the 'maxdist' keyword seemed like a good compromise.
  */
int Action_PairDist::SyncAction() {
  if (trajComm_.Size() > 1) {
    //rprintf("DEBUG: Total number of frames: %u\n", nframes_);
    // Ensure histogram bins have the correct count on this rank 
    UpdateHistogramFrames();
    // Get the total number of frames on all ranks.
    unsigned int myframes = nframes_;
    trajComm_.AllReduce(&nframes_, &myframes, 1, MPI_UNSIGNED, MPI_SUM);
/*
    // DEBUG - Print out histogram on each rank.
    for (int rank = 0; rank < trajComm_.Size(); rank++) {
      if ( rank == trajComm_.Rank() ) {
        rprintf("DEBUG: Histogram for rank.\n");
        for (unsigned long i = 0; i < histogram_.size(); i++) {
          double dist = ((double) i + 0.5) * delta_;
          rprintf("DEBUG:\t%12lu %12.4f %12.4f %12.4f %12.4f\n", i, dist, histogram_[i].mean(), histogram_[i].M2(), histogram_[i].nData());
        }
      }
      trajComm_.Barrier();
    }
    trajComm_.Barrier();
    // END DEBUG
*/
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
  }
  return 0;
}
#endif

// Action_PairDist::print()
void Action_PairDist::Print()
{
  // Ensure bins not present at the beginning are properly updated
  UpdateHistogramFrames();
  // Set DataSets X dim
  Dimension Xdim( 0.5 * delta_, delta_, "Distance" );
  Pr_->SetDim(Dimension::X, Xdim);
  std_->SetDim(Dimension::X, Xdim);

  //mprintf("DEBUG: Final result:\n");
  for (unsigned long i = 0; i < histogram_.size(); i++) {
    //mprintf("DEBUG:\t%12lu %12.4f %12.4f %12.4f %12.4f\n", i, ((double) i + 0.5) * delta_, histogram_[i].mean(), histogram_[i].M2(), histogram_[i].nData());
    double Pr = histogram_[i].mean() / delta_;

    if (Pr > 0.0) {
      double dist = ((double) i + 0.5) * delta_;
      double sd = sqrt(histogram_[i].variance() );
      ((DataSet_Mesh*) Pr_)->AddXY(dist, Pr);
      ((DataSet_Mesh*) std_)->AddXY(dist, sd);
    }
  }
}
