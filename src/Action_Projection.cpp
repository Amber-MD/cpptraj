#include <cmath> // sqrt
#include "Action_Projection.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "Constants.h" // DEGRAD

// CONSTRUCTOR
Action_Projection::Action_Projection() :
  modinfo_(0),
  beg_(0),
  end_(0)
{}

void Action_Projection::Help() const {
  mprintf("\t[<name>] evecs <evecs dataset> [out <outfile>] [beg <beg>] [end <end>]\n"
          "\t[<mask>] [dihedrals <dataset arg>]\n\t%s\n"
          "  Calculate projection along given eigenvectors.\n", ActionFrameCounter::HelpText);
}

// Action_Projection::Init()
Action::RetType Action_Projection::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get ibeg, iend, start, stop, offset
  // NOTE: Must get 'end' before InitFrameCounter since the latter checks for 'end'
  beg_ = actionArgs.getKeyInt("beg", 1) - 1;
  end_ = actionArgs.getKeyInt("end", 2);
  if (InitFrameCounter(actionArgs)) return Action::ERR;

  std::string modesname = actionArgs.GetStringKey("modes"); // For backwards compat.
  if (modesname.empty()) modesname = actionArgs.GetStringKey("evecs");
  if (modesname.empty()) {
    mprinterr("Error: No eigenvectors data set specified ('evecs <name>'). To load\n"
              "Error:   eigenvectors from a file use 'readdata <file>' prior to this command.\n");
    return Action::ERR;
  }
  // Check if DataSet exists
  modinfo_ = (DataSet_Modes*)init.DSL().FindSetOfType( modesname, DataSet::MODES );
  if (modinfo_ == 0) {
    // To preserve backwards compat., if no modes data set specified try to
    // load the data file.
    DataFile dataIn;
    dataIn.SetDebug( debugIn );
    if (dataIn.ReadDataOfType( modesname, DataFile::EVECS, init.DSL() ))
      return Action::ERR;
    modinfo_ = (DataSet_Modes*)init.DSL().FindSetOfType( modesname, DataSet::MODES );
    if (modinfo_ == 0) return Action::ERR;
  }
  if (modinfo_->Nmodes() < 1) {
    rprinterr("Error: modes set '%s' is empty.\n", modinfo_->legend());
    return Action::ERR;
  }
  // Check if beg and end are in bounds.
  if (end_ > modinfo_->Nmodes()) {
    mprintf("Warning: 'end' %i is greater than # evecs (%i); setting end to %i\n",
            end_, modinfo_->Nmodes(), modinfo_->Nmodes());
    end_ = modinfo_->Nmodes();
  }
  if (beg_ < 0 || beg_ >= end_) {
    mprinterr("Error: 'beg' %i out of bounds.\n", beg_+1);
    return Action::ERR;
  }

  // Check modes type
  if (modinfo_->Meta().ScalarType() != MetaData::COVAR &&
      modinfo_->Meta().ScalarType() != MetaData::MWCOVAR &&
      modinfo_->Meta().ScalarType() != MetaData::DIHCOVAR &&
      modinfo_->Meta().ScalarType() != MetaData::IDEA)
  {
    mprinterr("Error: evecs type is not COVAR, MWCOVAR, DIHCOVAR, or IDEA.\n");
    return Action::ERR;
  }

  // Output Filename
  DataFile* DF = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );

  // Get dihedral data sets 
  if (modinfo_->Meta().ScalarType() == MetaData::DIHCOVAR) {
    DihedralSets_.clear();
    DihedralSets_.AddTorsionSets( init.DSL().GetMultipleSets(actionArgs.GetStringKey("dihedrals")) );
    if ( DihedralSets_.empty() ) {
      mprinterr("Error: No valid data sets found.\n");
      return Action::ERR;
    } else if ((int)DihedralSets_.size() * 2 != modinfo_->VectorSize()) {
      mprinterr("Error: Number of dihedral data sets %zu does not correspond to"
                " number of eigenvectors %i\n", DihedralSets_.size()*2, modinfo_->VectorSize());
      return Action::ERR;
    } else if ((int)DihedralSets_.size() * 2 != modinfo_->NavgCrd()) {
      mprinterr("Error: Number of dihedral data sets %zu does not correspond to"
                " number of average elements %i\n", DihedralSets_.size()*2, modinfo_->NavgCrd());
      return Action::ERR;
    }
  } else {
    // Get mask
    if (mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  }

  // Set up data sets
  std::string setname = actionArgs.GetStringNext();
  if (setname.empty())
    setname = init.DSL().GenerateDefaultName("Proj");
  for (int mode = beg_; mode < end_; ++mode) {
    int imode = mode + 1;
    if (modinfo_->Meta().ScalarType() != MetaData::IDEA) { // COVAR, MWCOVAR
      DataSet* dout = init.DSL().AddSet( DataSet::FLOAT, MetaData(setname, imode) );
      if (dout == 0) {
        mprinterr("Error: Could not create output dataset for mode %i\n", imode);
        return Action::ERR;
      }
      dout->SetLegend("Mode"+integerToString(imode));
      project_.push_back( dout );
      if (DF != 0) DF->AddDataSet( dout );
    } else { // IDEA TODO: Error check
      project_.push_back( init.DSL().AddSet( DataSet::FLOAT, MetaData(setname, "X", imode) ) );
      if (DF != 0) DF->AddDataSet( project_.back() );
      project_.push_back( init.DSL().AddSet( DataSet::FLOAT, MetaData(setname, "Y", imode) ) );
      if (DF != 0) DF->AddDataSet( project_.back() );
      project_.push_back( init.DSL().AddSet( DataSet::FLOAT, MetaData(setname, "Z", imode) ) );
      if (DF != 0) DF->AddDataSet( project_.back() );
      project_.push_back( init.DSL().AddSet( DataSet::FLOAT, MetaData(setname, "R", imode) ) );
      if (DF != 0) DF->AddDataSet( project_.back() );
    }
  }
  // Set datafile args
  mprintf("    PROJECTION: Calculating projection using eigenvectors %i to %i of %s\n",
          beg_+1, end_, modinfo_->legend());
  if (DF != 0)
    mprintf("\tResults are written to %s\n", DF->DataFilename().full());
  FrameCounterInfo();
  if (modinfo_->Meta().ScalarType() == MetaData::DIHCOVAR)
    mprintf("\t%zu dihedral data sets.\n", DihedralSets_.size());
  else
    mprintf("\tAtom Mask: [%s]\n", mask_.MaskString());

  return Action::OK;
}

// Action_Projection::Setup()
Action::RetType Action_Projection::Setup(ActionSetup& setup) {
  if (modinfo_->Meta().ScalarType() != MetaData::DIHCOVAR) {
    // Setup mask
    if (setup.Top().SetupIntegerMask( mask_ )) return Action::ERR;
    if (mask_.None()) {
      mprintf("Warning: No atoms selected.\n");
      return Action::SKIP;
    }
    mask_.MaskInfo();
    // Check # of selected atoms against modes info
    if ( modinfo_->Meta().ScalarType() == MetaData::COVAR || 
         modinfo_->Meta().ScalarType() == MetaData::MWCOVAR)
    {
      // Check if (3 * number of atoms in mask) and nvectelem agree
      int natom3 = mask_.Nselected() * 3;
      if ( natom3 != modinfo_->NavgCrd() ) {
        mprinterr("Error: number selected coords (%i) != number avg coords (%i) in %s\n",
                  natom3, modinfo_->NavgCrd(), modinfo_->legend());
        return Action::ERR;
      }
      if ( natom3 != modinfo_->VectorSize() ) {
        mprinterr("Error: number selected coords (%i) != eigenvector size (%i)\n",
                  natom3, modinfo_->VectorSize() );
        return Action::ERR;
      }
    } else if ( modinfo_->Meta().ScalarType() == MetaData::IDEA ) {
      // Check if (number of atoms in mask) and nvectelem agree
      if (//mask_.Nselected() != modinfo_.Navgelem() ||
          mask_.Nselected() != modinfo_->VectorSize()) 
      {
        mprinterr("Error: number selected atoms (%i) != eigenvector size (%i)\n",
                  mask_.Nselected(), modinfo_->VectorSize() );
        return Action::ERR;
      }
    }

    // Precalc sqrt of mass for each coordinate
    sqrtmasses_.clear();
    if ( modinfo_->Meta().ScalarType() == MetaData::MWCOVAR ) {
      sqrtmasses_.reserve( mask_.Nselected() );
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
        sqrtmasses_.push_back( sqrt( setup.Top()[*atom].Mass() ) );
    } else {
      // If not MWCOVAR no mass-weighting necessary
      sqrtmasses_.resize( mask_.Nselected(), 1.0 );
    }
  }
  return Action::OK;
}

// Action_Projection::DoAction()
Action::RetType Action_Projection::DoAction(int frameNum, ActionFrame& frm) {
  if ( CheckFrameCounter( frm.TrajoutNum() ) ) return Action::OK;
  // Always start at first eigenvector element of first mode.
  const double* Vec = modinfo_->Eigenvector(beg_);
  // Project snapshots on modes
  if ( modinfo_->Meta().ScalarType() == MetaData::COVAR || 
       modinfo_->Meta().ScalarType() == MetaData::MWCOVAR ) 
  {
    for (int mode = beg_; mode < end_; ++mode) {
      DataSet_Modes::AvgIt Avg = modinfo_->AvgBegin();
      double proj = 0;
      std::vector<double>::const_iterator sqrtmass = sqrtmasses_.begin();
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      {
        const double* XYZ = frm.Frm().XYZ( *atom );
        double mass = *(sqrtmass++);
        proj += (XYZ[0] - *(Avg++)) * mass * Vec[0]; 
        proj += (XYZ[1] - *(Avg++)) * mass * Vec[1]; 
        proj += (XYZ[2] - *(Avg++)) * mass * Vec[2]; 
        Vec += 3;
      }
      float fproj = (float)proj;
      project_[mode]->Add( frameNum, &fproj );
    }
  } else if (modinfo_->Meta().ScalarType() == MetaData::DIHCOVAR ) {
    for (int mode = beg_; mode < end_; ++mode) {
      DataSet_Modes::AvgIt Avg = modinfo_->AvgBegin();
      double proj = 0.0;
      for (Array1D::const_iterator dih = DihedralSets_.begin();
                                   dih != DihedralSets_.end(); ++dih)
      {
        double theta = (*dih)->Dval( frm.TrajoutNum() ) * Constants::DEGRAD;
        proj += (cos(theta) - *(Avg++)) * Vec[0];
        proj += (sin(theta) - *(Avg++)) * Vec[1];
        Vec += 2;
      }
      // TODO: Convert to degrees?
      float fproj = (float)proj;
      project_[mode]->Add( frameNum, &fproj );
    }
  } else { // if modinfo_.ScalarType() == IDEA
    int ip = 0;
    for (int mode = beg_; mode < end_; ++mode) {
      double proj1 = 0;
      double proj2 = 0;
      double proj3 = 0;
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      {
        const double* XYZ = frm.Frm().XYZ(*atom);
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
