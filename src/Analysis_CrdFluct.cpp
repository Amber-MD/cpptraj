#include <cmath>
#include "Analysis_CrdFluct.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI
#include "StringRoutines.h" // integerToString

Analysis_CrdFluct::Analysis_CrdFluct() : 
  coords_(0),
  bfactor_(true),
  windowSize_(-1)
{}

void Analysis_CrdFluct::Help() {
  mprintf("crdfluct <crd set name> [out <filename>] [window <size>]\n");
}

Analysis::RetType Analysis_CrdFluct::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, int debugIn)
{
  std::string setname = analyzeArgs.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: crdfluct: Specify set name.\n");
    Help();
    return Analysis::ERR;
  }
  coords_ = (DataSet_Coords*)datasetlist->FindSetOfType( setname, DataSet::COORDS );
  if (coords_ == 0) {
    mprinterr("Error: crdfluct: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    return Analysis::ERR;
  }
  outfilename_ = analyzeArgs.GetStringKey("out");
  windowSize_ = analyzeArgs.getKeyInt("window", -1);

  //outset_ = datasetlist->AddSet( DataSet::DOUBLE, analyzeArgs.GetStringNext(), "fluct" );
  setname_ = analyzeArgs.GetStringNext();
  if (setname_.empty()) setname_ = datasetlist->GenerateDefaultName("fluct");

  mprintf("    CRDFLUCT: Atomic fluctuations will be calcd for set %s\n", 
          coords_->Legend().c_str());
  return Analysis::OK;
}

void Analysis_CrdFluct::CalcBfactors( Frame SumCoords, Frame SumCoords2, double Nsets,
                                      DataSet& outset )
{
  SumCoords.Divide(Nsets);
  SumCoords2.Divide(Nsets);
  //SumCoords2 = SumCoords2 - (SumCoords * SumCoords);
  SumCoords *= SumCoords;
  SumCoords2 -= SumCoords;
  AtomMask::const_iterator maskat = coords_->Mask().begin();
  if (bfactor_) {
    // Set up b factor normalization
    // B-factors are (8/3)*PI*PI * <r>**2 hence we do not sqrt the fluctuations
    double bfac = (8.0/3.0)*PI*PI;
    for (int i = 0; i < SumCoords2.size(); i+=3) {
      double fluct = (SumCoords2[i] + SumCoords2[i+1] + SumCoords2[i+2]) * bfac;
      outset.Add( *maskat, &fluct );
      ++maskat; 
    }
  } else {
    // Atomic fluctuations
    for (int i = 0; i < SumCoords2.size(); i+=3) {
      double fluct = SumCoords2[i] + SumCoords2[i+1] + SumCoords2[i+2];
      if (fluct > 0)
        outset.Add( *maskat, &fluct );
      ++maskat;
    }
  }
}

Analysis::RetType Analysis_CrdFluct::Analyze() {
  int end = coords_->Size();
  mprintf("\tFluctuation analysis for %i frames (%i atoms each).\n", end, 
          coords_->Natom());
  Frame currentFrame( coords_->Natom() );
  Frame SumCoords( coords_->Natom() );
  SumCoords.ZeroCoords();
  Frame SumCoords2( coords_->Natom() );
  SumCoords2.ZeroCoords();
  int w_count = 0;
  for (int frame = 0; frame < end; frame++) {
    currentFrame = (*coords_)[ frame ];
    SumCoords += currentFrame;
    SumCoords2 += ( currentFrame * currentFrame );
    ++w_count;
    if (w_count == windowSize_) {
      DataSet* ds = outSets_.AddSetIdx( DataSet::DOUBLE, setname_, frame + 1 );
      ds->SetLegend("F_" + integerToString( frame+1 ) );
      CalcBfactors( SumCoords, SumCoords2, (double)frame, *ds );
      w_count = 0;
    }
  }

  if (windowSize_ < 1) {
    DataSet* ds = outSets_.AddSet( DataSet::DOUBLE, setname_, "fluct" );
    CalcBfactors( SumCoords, SumCoords2, (double)end, *ds );
  } else if (w_count != 0) {
    mprintf("Warning: Number of frames (%i) was not evenly divisible by window size.\n",
             end);
    DataSet* ds = outSets_.AddSetIdx( DataSet::DOUBLE, setname_, end );
    ds->SetLegend("F_" + integerToString( end ) );
    CalcBfactors( SumCoords, SumCoords2, (double)end, *ds );
  }

  return Analysis::OK;
}

void Analysis_CrdFluct::Print( DataFileList* datafilelist ) {
  if (outfilename_.empty()) return;
  DataFile* outfile = datafilelist->AddDataFile(outfilename_);
  if (outfile == 0) return;
  for (DataSetList::const_iterator set = outSets_.begin();
                                   set != outSets_.end(); ++set)
    outfile->AddSet( *set );
  if (bfactor_)
    outfile->ProcessArgs("ylabel B-factors");
  outfile->ProcessArgs("xlabel Atom");
}
