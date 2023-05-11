#include "Action_Unwrap.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h"
#include "Image_List.h"
#include "DataSet_Mat3x3.h"
#ifdef _OPENMP
# include <omp.h>
#endif

// CONSTRUCTOR
Action_Unwrap::Action_Unwrap() :
  imageList_(0),
  imageMode_(Image::BYATOM),
  RefParm_(0),
  center_(false),
  refNeedsCalc_(false),
  avgucell_(0)
{ }

/** DESTRUCTOR */
Action_Unwrap::~Action_Unwrap() {
  if (imageList_ != 0) delete imageList_;
}

void Action_Unwrap::Help() const {
  mprintf("\t[center] [{byatom | byres | bymol}] [avgucell <avg ucell set>]\n"
          "\t[ %s ] [<mask>]\n", DataSetList::RefArgs);
  mprintf("  Reverse of 'image'; unwrap coordinates in <mask> according\n"
          "  to the first frame, or optionally a reference structure. Can\n"
          " unwrap by atom (default), by residue, or by molecule.\n");
}

// Action_Unwrap::Init()
Action::RetType Action_Unwrap::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'unwrap' action does not work with > 1 process (%i processes currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  // Get Keywords
  center_ = actionArgs.hasKey("center");
  if (actionArgs.hasKey("bymol"))
    imageMode_ = Image::BYMOL;
  else if (actionArgs.hasKey("byres"))
    imageMode_ = Image::BYRES;
  else if (actionArgs.hasKey("byatom")) {
    imageMode_ = Image::BYATOM;
    // Unwrapping to center by atom makes no sense
    if (center_) center_ = false;
  } else
    imageMode_ = Image::BYATOM;
  // Get reference
  ReferenceFrame REF = init.DSL().GetReferenceFrame( actionArgs );
  if (REF.error()) return Action::ERR;
  if (!REF.empty()) {
    RefFrame_ = REF.Coord();
    // Get reference parm for frame
    RefParm_ = REF.ParmPtr();
    refNeedsCalc_ = true;
  } else
    refNeedsCalc_ = false;
  // Get average box
  std::string avgucellstr = actionArgs.GetStringKey("avgucell");
  if (!avgucellstr.empty()) {
    avgucell_ = init.DSL().GetDataSet( avgucellstr );
    if (avgucell_ == 0) {
      mprinterr("Error: No data set selected by '%s'\n", avgucellstr.c_str());
      return Action::ERR;
    }
    if (avgucell_->Type() != DataSet::MAT3X3) {
      mprinterr("Error: Average unit cell set '%s' is not a 3x3 matrix set.\n", avgucell_->legend());
      return Action::ERR;
    }
    if (avgucell_->Size() < 1) {
      mprinterr("Error: Average unit cell set '%s' is empty.\n", avgucell_->legend());
      return Action::ERR;
    }
    DataSet_Mat3x3 const& matset = static_cast<DataSet_Mat3x3 const&>( *avgucell_ );
    Matrix_3x3 const& ucell = matset[0];
    if (avgbox_.SetupFromUcell( ucell )) {
      mprinterr("Error: Could not set up box from unit cell parameters in '%s'\n", avgucell_->legend());
      return Action::ERR;
    }
  }

  // Get mask string
  maskExpression_ = actionArgs.GetMaskNext();

  mprintf("    UNWRAP: By %s", Image::ModeString(imageMode_));
  if (!maskExpression_.empty())
    mprintf(" using mask '%s'", maskExpression_.c_str());
  else
    mprintf(" using all atoms");
  if (imageMode_ != Image::BYATOM) {
    if (center_)
      mprintf(" based on center of mass.");
    else
      mprintf(" based on first atom position.");
  }
  mprintf("\n");
  if ( !REF.empty())
    mprintf("\tReference is %s", REF.refName());
  else
    mprintf("\tReference is first frame.");
  mprintf("\n");
  if (avgucell_ != 0) {
    mprintf("\tUsing average unit cell vectors from set '%s' to remove box fluctuations.\n", avgucell_->legend());
    avgbox_.PrintInfo();
  }
  # ifdef _OPENMP
# pragma omp parallel
  {
# pragma omp master
  {
  mprintf("\tParallelizing calculation with %i threads.\n", omp_get_num_threads());
  }
  }
# endif

  return Action::OK;
}

// Action_Unwrap::Setup()
Action::RetType Action_Unwrap::Setup(ActionSetup& setup) {
  // Ensure same number of atoms in current parm and ref parm
  if ( RefParm_!=0 ) {
    if ( setup.Top().Natom() != RefParm_->Natom() ) {
      mprinterr("Error: unwrap: # atoms in reference parm %s is not\n", RefParm_->c_str());
      mprinterr("Error:         equal to # atoms in parm %s\n", setup.Top().c_str());
      return Action::ERR;
    }
  }
  // Check box type
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprintf("Error: unwrap: Parm %s does not contain box information.\n",
            setup.Top().c_str());
    return Action::ERR;
  }

  // Setup atom pairs to be unwrapped. Always use CoM TODO why?
  if (imageList_ != 0) delete imageList_;
  imageList_ = Image::CreateImageList(setup.Top(), imageMode_, maskExpression_,
                                      true, center_);
  if (imageList_ == 0) {
    mprinterr("Internal Error: Could not allocate unwrap list.\n");
    return Action::ERR;
  }
  if (imageList_->nEntities() < 1) {
    mprintf("Warning: Mask selects no atoms for topology '%s'.\n", setup.Top().c_str());
    return Action::SKIP;
  }
  mprintf("\tNumber of %ss to be unwrapped is %u\n",
          Image::ModeString(imageMode_), imageList_->nEntities());
  // Get entities that need to be updated in reference
  allEntities_ = imageList_->AllEntities();

  // Use current parm as reference if not already set
  if (RefParm_ == 0)
    RefParm_ = setup.TopAddress();
  return Action::OK;
}

// Action_Unwrap::DoAction()
Action::RetType Action_Unwrap::DoAction(int frameNum, ActionFrame& frm) {
  Matrix_3x3 const* ucell;
  if (avgucell_ == 0)
    ucell = &(frm.Frm().BoxCrd().UnitCell());
  else
    ucell = &(avgbox_.UnitCell());
  if (refNeedsCalc_) {
    // Calculate initial fractional coords from reference frame.
    Image::UnwrapFrac(fracCoords_, RefFrame_, *imageList_, *ucell, frm.Frm().BoxCrd().FracCell());
    refNeedsCalc_ = false;
  }

  Image::UnwrapFrac(fracCoords_, frm.ModifyFrm(), *imageList_, *ucell, frm.Frm().BoxCrd().FracCell());

  return Action::MODIFY_COORDS;
}
