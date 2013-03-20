#include <cmath> // sqrt
#include "Action_Projection.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString

// CONSTRUCTOR
Action_Projection::Action_Projection() :
  beg_(0),
  end_(0)
{}

void Action_Projection::Help() {
  mprintf("\tmodes <modesfile> out <outfile> [beg <beg>] [end <end>] [<mask>]\n");
  mprintf("\t%s\n", ActionFrameCounter::HelpText);
  mprintf("\tCalculate projection of coordinates along given eigenmodes.\n");
}

Action::RetType Action_Projection::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get ibeg, iend, start, stop, offset
  // NOTE: Must get 'end' before InitFrameCounter since the latter checks for 'end'
  beg_ = actionArgs.getKeyInt("beg", 1);
  end_ = actionArgs.getKeyInt("end", 2);
  if (InitFrameCounter(actionArgs)) return Action::ERR;

  // Get modes from file
  std::string modesname = actionArgs.GetStringKey("modes");
  if (modesname.empty()) {
    mprinterr("Error: projection: no modes file specified ('modes <filename>')\n");
    return Action::ERR;
  }
  if (modinfo_.ReadEvecFile( modesname, beg_, end_ )) return Action::ERR;

  // Check modes type
  if (modinfo_.Type() != DataSet_Matrix::COVAR &&
      modinfo_.Type() != DataSet_Matrix::MWCOVAR &&
      modinfo_.Type() != DataSet_Matrix::IDEA)
  {
    mprinterr("Error: projection: evecs type is not COVAR, MWCOVAR, or IDEA.\n");
    return Action::ERR;
  }

  // Output Filename
  std::string filename_ = actionArgs.GetStringKey("out");
  if (filename_.empty()) {
    mprinterr("Error: projection: No output file specified ('out <filename>')\n");
    return Action::ERR;
  }

  // Get mask
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  // Set up data sets
  std::string setname = actionArgs.GetStringNext();
  if (setname.empty())
    setname = DSL->GenerateDefaultName("Proj");
  int imode = beg_;
  for (int mode = 0; mode < modinfo_.Nmodes(); ++mode) {
    if (modinfo_.Type() != DataSet_Matrix::IDEA) { // COVAR, MWCOVAR
      DataSet* dout = DSL->AddSetIdx( DataSet::FLOAT, setname, imode );
      if (dout == 0) {
        mprinterr("Error creating output dataset for mode %i\n",imode);
        return Action::ERR;
      }
      dout->SetLegend("Mode"+integerToString(imode));
      project_.push_back( dout );
      DFL->AddSetToFile( filename_, dout );
      imode++;
    } else { // IDEA TODO: Error check
      project_.push_back( DSL->AddSetIdxAspect( DataSet::FLOAT, setname, imode, "X") );
      DFL->AddSetToFile( filename_, project_.back() );
      project_.push_back( DSL->AddSetIdxAspect( DataSet::FLOAT, setname, imode, "Y") );
      DFL->AddSetToFile( filename_, project_.back() );
      project_.push_back( DSL->AddSetIdxAspect( DataSet::FLOAT, setname, imode, "Z") );
      DFL->AddSetToFile( filename_, project_.back() );
      project_.push_back( DSL->AddSetIdxAspect( DataSet::FLOAT, setname, imode++, "R") );
      DFL->AddSetToFile( filename_, project_.back() );
    }
  }
  // Set datafile args
  if (!filename_.empty()) {
    DataFile* df = DFL->GetDataFile( filename_ );
    if (df == 0) {
      mprinterr("Internal Error: File %s was not set up.\n", filename_.c_str());
      return Action::ERR;
    }
    df->ProcessArgs("noemptyframes");
  }

  mprintf("    PROJECTION: Calculating projection using modes %i to %i of file %s\n",
          beg_, end_, modinfo_.Legend().c_str());
  mprintf("                Results are written to %s\n", filename_.c_str());
  FrameCounterInfo();
  mprintf("                Atom Mask: [%s]\n", mask_.MaskString());

  return Action::OK;
}

// Action_Projection::setup()
Action::RetType Action_Projection::Setup(Topology* currentParm, Topology** parmAddress) {
  // Setup mask
  if (currentParm->SetupIntegerMask( mask_ )) return Action::ERR;
  if (mask_.None()) {
    mprinterr("Error: projection: No atoms selected.\n");
    return Action::ERR;
  }
  mask_.MaskInfo();
  // Check # of selected atoms against modes info
  if ( modinfo_.Type() == DataSet_Matrix::COVAR || 
       modinfo_.Type() == DataSet_Matrix::MWCOVAR)
  {
    // Check if (3 * number of atoms in mask) and nvectelem agree
    int natom3 = mask_.Nselected() * 3;
    if ( natom3 != modinfo_.NavgCrd() ) {
      mprinterr("Error: projection: # selected coords (%i) != # avg coords (%i) in %s\n",
                natom3, modinfo_.NavgCrd(), modinfo_.Legend().c_str());
      return Action::ERR;
    }
    if ( natom3 != modinfo_.VectorSize() ) {
      mprinterr("Error: projection: # selected coords (%i) != eigenvector size (%i)\n",
                natom3, modinfo_.VectorSize() );
      return Action::ERR;
    }
  } else if ( modinfo_.Type() == DataSet_Matrix::IDEA ) {
    // Check if (number of atoms in mask) and nvectelem agree
    if (//mask_.Nselected() != modinfo_.Navgelem() ||
        mask_.Nselected() != modinfo_.VectorSize()) 
    {
      mprinterr("Error: projection: # selected atoms (%i) != eigenvector size (%i)\n",
                mask_.Nselected(), modinfo_.VectorSize() );
      return Action::ERR;
    }
  }

  // Precalc sqrt of mass for each coordinate
  sqrtmasses_.clear();
  if ( modinfo_.Type() == DataSet_Matrix::MWCOVAR ) {
    sqrtmasses_.reserve( mask_.Nselected() );
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      sqrtmasses_.push_back( sqrt( (*currentParm)[*atom].Mass() ) );
  } else {
    // If not MWCOVAR no mass-weighting necessary
    sqrtmasses_.resize( mask_.Nselected(), 1.0 );
  }

  return Action::OK;
}

// Action_Projection::action()
Action::RetType Action_Projection::DoAction(int frameNum, Frame* currentFrame, 
                                            Frame** frameAddress)
{
  if ( CheckFrameCounter( frameNum ) ) return Action::OK;
  // Always start at first eigenvector element of first mode.
  const double* Vec = modinfo_.Eigenvector(0);
  // Project snapshots on modes
  if ( modinfo_.Type() == DataSet_Matrix::COVAR || 
       modinfo_.Type() == DataSet_Matrix::MWCOVAR ) 
  {
    for (int mode = 0; mode < modinfo_.Nmodes(); ++mode) {
      const double* Avg = modinfo_.AvgCrd();
      double proj = 0;
      std::vector<double>::const_iterator sqrtmass = sqrtmasses_.begin();
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      {
        const double* XYZ = currentFrame->XYZ( *atom );
        double mass = *(sqrtmass++);
        proj += (XYZ[0] - Avg[0]) * mass * Vec[0]; 
        proj += (XYZ[1] - Avg[1]) * mass * Vec[1]; 
        proj += (XYZ[2] - Avg[2]) * mass * Vec[2]; 
        Avg += 3;
        Vec += 3;
      }
      float fproj = (float)proj;
      project_[mode]->Add( frameNum, &fproj );
    }
  } else { // if modinfo_.Type() == IDEA
    int ip = 0;
    for (int mode = 0; mode < modinfo_.Nmodes(); ++mode) {
      double proj1 = 0;
      double proj2 = 0;
      double proj3 = 0;
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      {
        const double* XYZ = currentFrame->XYZ(*atom);
        proj1 += XYZ[0] * *(Vec  );
        proj2 += XYZ[1] * *(Vec  );
        proj3 += XYZ[2] * *(Vec++);
      }
      float fproj1 = (float)proj1;
      float fproj2 = (float)proj2;
      float fproj3 = (float)proj3;
      double proj4 = sqrt(proj1*proj1 + proj2*proj2 + proj3*proj3);
      float fproj4 = (float)proj4;
      project_[ip++]->Add( frameNum, &fproj1 );
      project_[ip++]->Add( frameNum, &fproj2 );
      project_[ip++]->Add( frameNum, &fproj3 );
      project_[ip++]->Add( frameNum, &fproj4 );
    }
  }
  return Action::OK;
}
