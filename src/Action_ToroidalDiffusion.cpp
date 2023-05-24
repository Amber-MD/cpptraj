#include "Action_ToroidalDiffusion.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Action_ToroidalDiffusion::Action_ToroidalDiffusion() :
  useMass_(false)
{}

// Action_ToroidalDiffusion::Help()
void Action_ToroidalDiffusion::Help() const {
  mprintf("\t[<mask>] [mass]\n");
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
  useMass_ = actionArgs.hasKey("mass");
  if (mask1_.SetMaskString( actionArgs.GetMaskNext() )) {
    mprinterr("Error: Invalid mask string.\n");
    return Action::ERR;
  }

  mprintf("    TORDIFF: Toroidal-view-preserving diffusion calculation.\n");
  mprintf("\tCalculating diffusion for molecules selected by mask '%s'\n", mask1_.MaskString());
  if (useMass_)
    mprintf("\tUsing center of mass.\n");
  else
    mprintf("\tUsing geometric center.\n");
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
  if (torPositions_.empty()) {
    torPositions_.reserve( entities_.size() );
    if (useMass_) {
      for (Marray::const_iterator mask = entities_.begin(); mask != entities_.end(); ++mask)
        torPositions_.push_back( frm.Frm().VCenterOfMass( *mask ) );
    } else {
      for (Marray::const_iterator mask = entities_.begin(); mask != entities_.end(); ++mask)
        torPositions_.push_back( frm.Frm().VGeometricCenter( *mask ) );
    }
  }
    
  return Action::OK;
}
