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
  avgucell_(0),
  scheme_(FRAC)
{ }

/** DESTRUCTOR */
Action_Unwrap::~Action_Unwrap() {
  if (imageList_ != 0) delete imageList_;
}

void Action_Unwrap::Help() const {
  mprintf("\t[center] [{byatom | byres | bymol}] [avgucell <avg ucell set>]\n"
          "\t[ %s ] [<mask>]\n"
          "\t[scheme {frac|tor}]\n",
          DataSetList::RefArgs);
  mprintf("  Reverse of 'image'; unwrap coordinates in <mask> according\n"
          "  to the first frame, or optionally a reference structure. Can\n"
          "  unwrap by atom (default), by residue, or by molecule.\n");
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

  std::string sarg = actionArgs.GetStringKey("scheme");
  if (sarg.empty())
    scheme_ = FRAC;
  else if (sarg == "frac")
    scheme_ = FRAC;
  else if (sarg == "tor")
    scheme_ = TOR;
  else {
    mprinterr("Error: Unrecognized 'scheme' %s\n", sarg.c_str());
    return Action::ERR;
  }

  Image::Mode defaultMode;
  switch (scheme_) {
    case FRAC : defaultMode = Image::BYATOM; break;
    case TOR  : defaultMode = Image::BYMOL; break;
  }

  if (actionArgs.hasKey("bymol"))
    imageMode_ = Image::BYMOL;
  else if (actionArgs.hasKey("byres"))
    imageMode_ = Image::BYRES;
  else if (actionArgs.hasKey("byatom")) {
    imageMode_ = Image::BYATOM;
    // Unwrapping to center by atom makes no sense
    if (center_) center_ = false;
  } else
    imageMode_ = defaultMode;

  if (scheme_ == TOR && imageMode_ != Image::BYMOL) { // TODO warning only?
    mprinterr("Error: Toroidal-view-preserving unwrap only works with 'bymol'.\n");
    return Action::ERR;
  }

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
  switch (scheme_) {
    case FRAC : mprintf("\tUnwrapping using fractional coordinates.\n"); break;
    case TOR  : mprintf("\tUnwrapping using toroidal-view-preserving scheme.\n"); break;
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
  if (scheme_ == TOR) {
    mprintf("# Citation: Bullerjahn, von Bulow, Heidari, Henin, and Hummer.\n"
            "#           \"Unwrapping NPT Simulations to Calculate Diffusion Coefficients.\n"
            "#           https://arxiv.org/abs/2303.09418\n");
  }

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
  if (scheme_ == TOR) {
    if (!setup.CoordInfo().TrajBox().Is_X_Aligned_Ortho()) {
      mprinterr("Error: Toroidal-preserving-view unwrap calculation currently only works\n"
                "Error:   for X-aligned orthogonal cells.\n");
      return Action::ERR;
    }
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
    switch (scheme_) {
      case FRAC : Image::UnwrapFrac(fracCoords_, RefFrame_, *imageList_, *ucell, RefFrame_.BoxCrd().FracCell()); break;
      case TOR  : Image::UnwrapToroidal(torPositions_, fracCoords_, RefFrame_, *imageList_, RefFrame_.BoxCrd().Lengths()); break;
    }
    refNeedsCalc_ = false;
  }
  switch (scheme_) {
    case FRAC : Image::UnwrapFrac(fracCoords_, frm.ModifyFrm(), *imageList_, *ucell, frm.Frm().BoxCrd().FracCell()); break;
    case TOR  : Image::UnwrapToroidal(torPositions_, fracCoords_, frm.ModifyFrm(), *imageList_, frm.Frm().BoxCrd().Lengths()); break;
  }

  return Action::MODIFY_COORDS;
}
