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
  same_mask_(false),
  ub1_(0),
  ub2_(0)
{}

void Action_PairDist::Help() const {
  mprintf("\t[out <filename>] mask <mask> [mask2 <mask>] [delta <resolution>]\n"
          "\t[maxdist <distance>]\n"
          "  Calculate pair distribution function P(r) between two masks.\n"
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
  InitImaging(true);
  delta_ = actionArgs.getKeyDouble("delta", 0.01);
  double maxDist = actionArgs.getKeyDouble("maxdist", -1.0);
  if (maxDist > 0.0) {
    maxbin_ = (unsigned long)(maxDist / delta_);
    histogram_.resize( maxbin_ + 1 );
  }

  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );

  std::string mask1 = actionArgs.GetStringKey("mask");

  if (mask1.empty()) {
    mprinterr("Error: pairdist: No mask1 specified.\n");
    return Action::ERR;
  }

  mask1_.SetMaskString(mask1);

  std::string mask2 = actionArgs.GetStringKey("mask2");

  if (mask2.empty()) {
    same_mask_ = true;
    mask2_.SetMaskString(mask1);
  } else {
    mask2_.SetMaskString(mask2);

    if (mask1_.MaskExpression() != mask2_.MaskExpression() )
      same_mask_ = false;
    else
      same_mask_ = true;
  }

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
  if (!same_mask_)
    mprintf(" between atoms selected by '%s' and '%s'\n", mask1_.MaskString(), mask2_.MaskString());
  else
    mprintf(" for atoms selected by '%s'\n", mask1_.MaskString());
  if (outfile != 0)
    mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());
  mprintf("\tResolution is %f Ang.\n", delta_);
  if (maxDist > 0.0)
    mprintf("\tInitial histogram max distance= %f Ang\n", maxDist);
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
  if (setup.Top().SetupIntegerMask(mask1_) ) return Action::ERR;

  mprintf("\t");
  mask1_.BriefMaskInfo();

  if (mask1_.None()) {
    mprintf("Warning: Mask has no atoms.\n");
    return Action::SKIP;
  }

  if (setup.Top().SetupIntegerMask(mask2_) ) return Action::ERR;

  mask2_.BriefMaskInfo();
  mprintf("\n");

  if (mask2_.None()) {
    mprintf("Warning: PairDist::setup: Mask2 has no atoms.\n");
    return Action::SKIP;
  }

  if (mask1_.MaskExpression() != mask2_.MaskExpression() &&
      mask1_.NumAtomsInCommon(mask2_) > 0) {
    mprinterr("Error: mask expressions must be either "
	      "exactly the same\n\t(equivalent to mask2 omitted) or masks must "
	      "be non-overlapping.\n");
    return Action::ERR;
  }

  if (same_mask_) {
    ub1_ = mask1_.Nselected() - 1;
    ub2_ = mask1_.Nselected();
  } else {
    ub1_ = mask1_.Nselected();
    ub2_ = mask2_.Nselected();
  }

  SetupImaging(setup.CoordInfo().TrajBox().Type() );

  return Action::OK;
}


// Action_PairDist::action()
Action::RetType Action_PairDist::DoAction(int frameNum, ActionFrame& frm) {
  unsigned long bin, j;
  double Dist = 0.0;
  Matrix_3x3 ucell, recip;
  Vec3 a1, a2;
  std::vector<double> tmp;	// per frame histogram


  tmp.resize(histogram_.size() );

  for (unsigned long i = 0; i < ub1_; i++) {
    for (same_mask_ ? j = i + 1 : j = 0; j < ub2_; j++) {
      a1 = frm.Frm().XYZ(mask1_[i]);
      a2 = frm.Frm().XYZ(mask2_[j]);

      switch (ImageType() ) {
      case NONORTHO:
	frm.Frm().BoxCrd().ToRecip(ucell, recip);
	Dist = DIST2_ImageNonOrtho(a1, a2, ucell, recip);
	break;
      case ORTHO:
	Dist = DIST2_ImageOrtho(a1, a2, frm.Frm().BoxCrd());
	break;
      case NOIMAGE:
	Dist = DIST2_NoImage(a1, a2);
	break;
      }

      bin = (unsigned long) (sqrt(Dist) / delta_);

      if (bin > maxbin_) {
	maxbin_ = bin;
	tmp.resize(maxbin_ + 1);
	histogram_.resize(maxbin_ + 1);
      }

      tmp[bin]++;
    }
  }

  // FIXME: may be inefficient to call accumulate() on every data point
  // -> pass "array" to accumulate()?
  for (unsigned long i = 0; i < tmp.size(); i++) {
    histogram_[i].accumulate(tmp[i]);
  }

  return Action::OK;
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
  // Set DataSets X dim
  Dimension Xdim( 0.5 * delta_, delta_, "Distance" );
  Pr_->SetDim(Dimension::X, Xdim);
  std_->SetDim(Dimension::X, Xdim);

  double dist, Pr, sd;

  for (unsigned long i = 0; i < histogram_.size(); i++) {
    Pr = histogram_[i].mean() / delta_;

    if (Pr > 0.0) {
      dist = ((double) i + 0.5) * delta_;
      sd = sqrt(histogram_[i].variance() );
      ((DataSet_Mesh*) Pr_)->AddXY(dist, Pr);
      ((DataSet_Mesh*) std_)->AddXY(dist, sd);
    }
  }
}
