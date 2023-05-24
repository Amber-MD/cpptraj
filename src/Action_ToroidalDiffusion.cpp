#include "Action_ToroidalDiffusion.h"
#include "CpptrajStdio.h"
#include <cmath> // floor

/** CONSTRUCTOR */
Action_ToroidalDiffusion::Action_ToroidalDiffusion() :
  useMass_(false),
  avg_x_(0),
  avg_y_(0),
  avg_z_(0),
  avg_r_(0),
  avg_a_(0),
  time_(0)
{}

// Action_ToroidalDiffusion::Help()
void Action_ToroidalDiffusion::Help() const {
  mprintf("\t[<set name>] [<mask>] [mass]\n");
}

// Action_ToroidalDiffusion::Init()
Action::RetType Action_ToroidalDiffusion::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'tordiff' action does not work with > 1 process (%i processes currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out") );
  useMass_ = actionArgs.hasKey("mass");
  time_ = actionArgs.getNextDouble(1.0);
  if (mask1_.SetMaskString( actionArgs.GetMaskNext() )) {
    mprinterr("Error: Invalid mask string.\n");
    return Action::ERR;
  }
  // Add DataSets
  std::string dsname_ = actionArgs.GetStringNext();
  if (dsname_.empty())
    dsname_ = init.DSL().GenerateDefaultName("Diff");
  avg_x_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "X"));
  avg_y_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Y"));
  avg_z_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Z"));
  avg_r_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "R"));
  avg_a_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "A"));
  if (avg_x_ == 0 || avg_y_ == 0 || avg_z_ == 0 || avg_r_ == 0 || avg_a_ == 0) {
    mprinterr("Error: Could not allocate one or more average toroidal diffusion sets.\n");
    return Action::ERR;
  }
  if (outfile != 0) {
    outfile->AddDataSet( avg_r_ );
    outfile->AddDataSet( avg_x_ );
    outfile->AddDataSet( avg_y_ );
    outfile->AddDataSet( avg_z_ );
    outfile->AddDataSet( avg_a_ );
  }
  // Set X dim
  Dimension Xdim_ = Dimension(0.0, time_, "Time");
  avg_x_->SetDim(Dimension::X, Xdim_);
  avg_y_->SetDim(Dimension::X, Xdim_);
  avg_z_->SetDim(Dimension::X, Xdim_);
  avg_r_->SetDim(Dimension::X, Xdim_);
  avg_a_->SetDim(Dimension::X, Xdim_);

  mprintf("    TORDIFF: Toroidal-view-preserving diffusion calculation.\n");
  mprintf("\tCalculating diffusion for molecules selected by mask '%s'\n", mask1_.MaskString());
  if (useMass_)
    mprintf("\tUsing center of mass.\n");
  else
    mprintf("\tUsing geometric center.\n");
  mprintf("\tThe time between frames is %g ps.\n", time_);
  mprintf("\tData set name: %s\n", dsname_.c_str());
  if (outfile != 0)
    mprintf("\tOutput to file '%s'\n", outfile->DataFilename().full());
  mprintf("# Citation: Bullerjahn, von Bulow, Heidari, Henin, and Hummer.\n"
          "#           \"Unwrapping NPT Simulations to Calculate Diffusion Coefficients.\n"
          "#           https://arxiv.org/abs/2303.09418\n");

  return Action::OK;
}

/** \return Array of masks selecting molecules to track. */
Action_ToroidalDiffusion::Marray Action_ToroidalDiffusion::setup_entities(Topology const& topIn)
const
{
  Marray entitiesOut;
  std::vector<bool> molSelected(topIn.Nmol(), false);

  unsigned int molIdx = 0;
  for (Topology::mol_iterator mol = topIn.MolStart(); mol != topIn.MolEnd(); ++mol, ++molIdx)
  {
    for (Unit::const_iterator seg = mol->MolUnit().segBegin(); seg != mol->MolUnit().segEnd(); ++seg)
    {
      for (int at = seg->Begin(); at != seg->End(); at++)
      {
        if (mask1_.AtomInCharMask( at )) {
          if (!molSelected[molIdx]) {
            entitiesOut.push_back( AtomMask() );
            molSelected[molIdx] = true;
          }
          entitiesOut.back().AddAtom( at );
        }
      } // END loop over atoms in segment
    } // END loop over segments in molecule
  } // END loop over molecules

  mprintf("DEBUG: Entities selected by '%s'\n", mask1_.MaskString());
  for (Marray::const_iterator it = entitiesOut.begin(); it != entitiesOut.end(); ++it) {
    mprintf("\t");
    for (AtomMask::const_iterator at = it->begin(); at != it->end(); ++at)
      mprintf(" %i", *at + 1);
    mprintf("\n");
  }

  return entitiesOut;
}

// Action_ToroidalDiffusion::Setup()
Action::RetType Action_ToroidalDiffusion::Setup(ActionSetup& setup)
{
  if (setup.Top().SetupCharMask( mask1_ )) {
    mprinterr("Error: Could not set up mask '%s'\n", mask1_.MaskString());
    return Action::ERR;
  }
  mask1_.MaskInfo();
  if (mask1_.None()) {
    mprintf("Warning: No atoms selected by '%s'\n", mask1_.MaskString());
    return Action::SKIP;
  }
  if (entities_.empty()) {
    // First time set up.
    entities_ = setup_entities( setup.Top() );
    if (entities_.empty()) {
      mprinterr("Error: No entities to calculate diffusion for.\n");
      return Action::ERR;
    }
  } else {
    // Require that the same number of entities be selected TODO check mol #s?
    Marray newEntities = setup_entities( setup.Top() );
    if (newEntities.size() != entities_.size()) {
      mprinterr("Error: Number of entities selected for topology '%s' (%zu)\n"
                "Error:   differs from original number of entities (%zu)\n",
                setup.Top().c_str(), newEntities.size(), entities_.size());
      return Action::ERR;
    }
  }
  // Check for box
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprintf("Error: Topology '%s' does not contain box information.\n",
            setup.Top().c_str());
    return Action::ERR;
  }
  if (!setup.CoordInfo().TrajBox().Is_X_Aligned_Ortho()) {
    mprinterr("Error: Toroidal-preserving-view diffusion calculation currently only works\n"
              "Error:   for X-aligned orthogonal cells.\n");
    return Action::ERR;
  }

  return Action::OK;
}

// Action_ToroidalDiffusion::DoAction()
Action::RetType Action_ToroidalDiffusion::DoAction(int frameNum, ActionFrame& frm)
{
  Box const& currentBox = frm.Frm().BoxCrd();
  if (!currentBox.HasBox() || !currentBox.Is_X_Aligned_Ortho()) {
    mprinterr("Error: Either no box or box not X-aligned ortho. for frame %i\n", frameNum+1);
    return Action::ERR;
  }
  Vec3 boxVec = currentBox.Lengths();

  double average2 = 0.0;
  double avgx = 0.0;
  double avgy = 0.0;
  double avgz = 0.0;

  if (torPositions_.empty()) {
    // ----- First frame ---------------
    torPositions_.reserve( entities_.size() );
    prevPositions_.reserve( entities_.size() );
    initialPositions_.reserve( entities_.size() );
    if (useMass_) {
      for (Marray::const_iterator mask = entities_.begin(); mask != entities_.end(); ++mask) {
        torPositions_.push_back( frm.Frm().VCenterOfMass( *mask ) );
        prevPositions_.push_back( torPositions_.back() );
        initialPositions_.push_back( torPositions_.back() );
      }
    } else {
      for (Marray::const_iterator mask = entities_.begin(); mask != entities_.end(); ++mask) {
        torPositions_.push_back( frm.Frm().VGeometricCenter( *mask ) );
        prevPositions_.push_back( torPositions_.back() );
        initialPositions_.push_back( torPositions_.back() );
      }
    }
  } else {
    // ----- Subsequent frames ---------
    int idx;
    int maxidx = (int)entities_.size();
#   ifdef _OPENMP
#   pragma omp parallel private(idx) reduction(+ : average2, avgx, avgy, avgz)
    {
#   pragma omp for
#   endif
    for (idx = 0; idx < maxidx; idx++)
    {
      // wi+1 - wi
      Vec3 Wi1;
      if (useMass_)
        Wi1 = frm.Frm().VCenterOfMass( entities_[idx] );
      else
        Wi1 = frm.Frm().VGeometricCenter( entities_[idx] );
      Vec3 deltaW = Wi1 - prevPositions_[idx];
      // Calculate translation for toroidal scheme (3rd term of eq. 2)
      Vec3 trans;
      trans[0] = floor( (deltaW[0] / boxVec[0]) + 0.5 ) * boxVec[0];
      trans[1] = floor( (deltaW[1] / boxVec[1]) + 0.5 ) * boxVec[1];
      trans[2] = floor( (deltaW[2] / boxVec[2]) + 0.5 ) * boxVec[2];
      // Calculate current position in toroidal scheme
      Vec3 Ui1 = torPositions_[idx] + deltaW - trans;
      // Calculate distance from current toroidal position to initial position
      Vec3 delta = Ui1 - initialPositions_[idx];
      // Calculate diffusion
      double distx = delta[0] * delta[0];
      double disty = delta[1] * delta[1];
      double distz = delta[2] * delta[2];
      double dist2 = distx + disty + distz;
      // Accumulate averages
      avgx += distx;
      avgy += disty;
      avgz += distz;
      average2 += dist2;

      // Update arrays
      torPositions_[idx] = Ui1;
      prevPositions_[idx] = Wi1;
    } // END loop over entities
#   ifdef _OPENMP
    } /* END pragma omp parallel */
#   endif
  }
  // Calc averages
  double dNselected = 1.0 / (double)entities_.size();
  avgx *= dNselected;
  avgy *= dNselected;
  avgz *= dNselected;
  average2 *= dNselected;
  // Save averages
  avg_x_->Add(frameNum, &avgx);
  avg_y_->Add(frameNum, &avgy);
  avg_z_->Add(frameNum, &avgz);
  avg_r_->Add(frameNum, &average2);
  average2 = sqrt(average2);
  avg_a_->Add(frameNum, &average2);

  return Action::OK;
}
