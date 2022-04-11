#include <cmath>
#include "Analysis_CrdFluct.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI
#include "StringRoutines.h" // integerToString

Analysis_CrdFluct::Analysis_CrdFluct() : 
  coords_(0),
  bfactor_(false),
  windowSize_(-1)
{}

void Analysis_CrdFluct::Help() const {
  mprintf("\t[crdset <crd set>] [<mask>] [out <filename>] [window <size>] [bfactor]\n"
          "  Calculate atomic positional fluctuations for atoms in <mask>\n"
          "  over windows of specified size.\n"
          "  <crd set> can be created with the 'createcrd' command.\n");
}

// Analysis_CrdFluct::Setup()
Analysis::RetType Analysis_CrdFluct::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  bfactor_ = analyzeArgs.hasKey("bfactor");
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: crdfluct: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    return Analysis::ERR;
  }
  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  windowSize_ = analyzeArgs.getKeyInt("window", -1);
  // Get mask
  if (mask_.SetMaskString( analyzeArgs.GetMaskNext() )) return Analysis::ERR;

  mprintf("    CRDFLUCT: Atomic fluctuations will be calcd for set %s, mask [%s]\n", 
          coords_->legend(), mask_.MaskString());
  if (windowSize_ != -1) mprintf("\tWindow size = %i\n", windowSize_);
  if (outfile != 0) mprintf("\tOutput to %s\n", outfile->DataFilename().base());

  // Set up data sets
  setname = analyzeArgs.GetStringNext();
  if (windowSize_ < 1) {
    // Only one data set for total B-factors
    DataSet* ds = setup.DSL().AddSet( DataSet::DOUBLE, setname, "fluct" );
    if (ds == 0) return Analysis::ERR;
    outSets_.push_back( ds );
    if (outfile != 0) outfile->AddDataSet( ds );
  } else {
    if (coords_->Size() == 0) {
      mprinterr("Error: window size > 0 and COORDS data set %s is empty.\n", 
                 coords_->legend());
      mprinterr("Error: Cannot predict how many window data sets will be needed.\n");
      return Analysis::ERR;
    }
    if (setname.empty()) setname = setup.DSL().GenerateDefaultName("fluct");
    // Determine how many windows will be needed
    int nwindows = coords_->Size() / windowSize_;
    for (int win = 0; win < nwindows; ++win) {
      int frame = (win + 1) * windowSize_;
      DataSet* ds = setup.DSL().AddSet( DataSet::DOUBLE, MetaData(setname, frame) );
      if (ds == 0) return Analysis::ERR;
      ds->SetLegend( "F_" + integerToString( frame ) );
      ds->SetDim( Dimension::X, Dimension(1.0, 1.0, "Atom") );
      outSets_.push_back( ds );
      if (outfile != 0) outfile->AddDataSet( ds );
    }
    if ( (coords_->Size() % windowSize_) != 0 ) {
      DataSet* ds = setup.DSL().AddSet( DataSet::DOUBLE, MetaData(setname, coords_->Size()) );
      ds->SetLegend("Final");
      outSets_.push_back( ds );
      if (outfile != 0) outfile->AddDataSet( ds );
    }
    for (SetList::iterator out = outSets_.begin(); out != outSets_.end(); ++out)
      mprintf("\t%s\n", (*out)->legend());
  }
  // Setup output file
/*  if (bfactor_)
    outfile->Dim(Dimension::Y).SetLabel("B-factors");
  outfile->Dim(Dimension::X).SetLabel("Atom");*/

  return Analysis::OK;
}

// Analysis_CrdFluct::CalcBfactors()
void Analysis_CrdFluct::CalcBfactors( Frame const& SumCoords, Frame const& SumCoords2, double Nsets,
                                      DataSet& outset )
{
  for (int ix = 0; ix != SumCoords.size(); ix++) {
    sum_[ix]  = SumCoords[ix]  / Nsets;
    sum2_[ix] = SumCoords2[ix] / Nsets;
  }
  //sum2_ = sum2_ - (sum_ * sum_);
  sum_ *= sum_;
  sum2_ -= sum_;
  AtomMask::const_iterator maskat = mask_.begin();
  if (bfactor_) {
    // Set up b factor normalization
    // B-factors are (8/3)*PI*PI * <r>**2 hence we do not sqrt the fluctuations
    double bfac = (8.0/3.0)*Constants::PI*Constants::PI;
    for (int i = 0; i < sum2_.size(); i+=3) {
      double fluct = (sum2_[i] + sum2_[i+1] + sum2_[i+2]) * bfac;
      outset.Add( *maskat, &fluct );
      ++maskat; 
    }
  } else {
    // Atomic fluctuations
    for (int i = 0; i < sum2_.size(); i+=3) {
      double fluct = sum2_[i] + sum2_[i+1] + sum2_[i+2];
      if (fluct > 0)
        outset.Add( *maskat, &fluct );
      ++maskat;
    }
  }
}

// Analysis_CrdFluct::Analyze()
Analysis::RetType Analysis_CrdFluct::Analyze() {
  // Set up mask
  if ( coords_->Top().SetupIntegerMask( mask_ )) return Analysis::ERR;
  mask_.MaskInfo();
  if (mask_.None()) return Analysis::ERR;
  int end = coords_->Size();
  mprintf("\tFluctuation analysis for %i frames (%i atoms each).\n", end, 
          coords_->Top().Natom());
  Frame currentFrame( mask_.Nselected() );
  Frame SumCoords( mask_.Nselected() );
  SumCoords.ZeroCoords();
  Frame SumCoords2( mask_.Nselected() );
  SumCoords2.ZeroCoords();
  sum_ = SumCoords;
  sum2_ = SumCoords2;
  int w_count = 0;
  SetList::iterator out = outSets_.begin();
  for (int frame = 0; frame < end; frame++) {
    coords_->GetFrame(frame, currentFrame, mask_);
    SumCoords += currentFrame;
    SumCoords2 += ( currentFrame * currentFrame );
    ++w_count;
    if (w_count == windowSize_) {
      CalcBfactors( SumCoords, SumCoords2, (double)frame, *(*out) );
      ++out;
      w_count = 0;
    }
  }

  if (windowSize_ < 1 || w_count != 0) {
    // For windowSize < 1 this is the only b-factor calc
    CalcBfactors( SumCoords, SumCoords2, (double)end, *(*out) );
    if (w_count != 0) 
      mprintf("Warning: Number of frames (%i) was not evenly divisible by window size.\n",
               end);
  }

  return Analysis::OK;
}
