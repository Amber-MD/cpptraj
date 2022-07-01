#include "Action_Rotate.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "DataSet_Mat3x3.h"

/** CONSTRUCTOR */
Action_Rotate::Action_Rotate() :
  rmatrices_(0),
  delta_(0.0),
  mode_(ROTATE),
  inverse_(false),
  all_atoms_selected_(false),
  dsout_tx_(0),
  dsout_ty_(0),
  dsout_tz_(0),
  dsout_t_(0)
{ }

/** Action help. */
void Action_Rotate::Help() const {
  mprintf("\t[<mask>] { [x <xdeg>] [y <ydeg>] [z <zdeg>]  |\n"
          "\t           axis0 <mask0> axis1 <mask1> <deg> |\n"
          "\t           usedata <set name> [inverse] |\n"
          "\t           calcfrom <set name> [name <output set name>] [out <file>]\n"
          "\t         }\n"
          "  Rotate atoms in <mask> either around the x, y, and/or z axes, around the\n"
          "  the axis defined by <mask0> to <mask1>, or using rotation matrices in\n"
          "  the specified data set.\n"
          "  If 'calcfrom' is specified, instead of rotating the system, calcuate\n"
          "  rotations around the X Y and Z axes from previously calculated\n"
          "  rotation matrices in specified data set.\n");
}

/** \return 3x3 matrix DataSet. */
int Action_Rotate::Get3x3Set(DataSetList const& DSL, std::string const& dsname) {
  rmatrices_ = (DataSet_Mat3x3*)DSL.FindSetOfType( dsname, DataSet::MAT3X3 );
  if (rmatrices_ == 0) {
    mprinterr("Error: No 3x3 matrices data set '%s'\n", dsname.c_str());
    return 1;
  }
  return 0;
}

/** Print error message if set is null. */
static inline int CheckSet(DataSet* ds, std::string const& name, const char* aspect) {
  if (ds == 0) {
    mprinterr("Error: Could not set up output set %s[%s]\n", name.c_str(), aspect);
    return 1;
  }
  return 0;
}

/** Create output DataSets */
int Action_Rotate::SetupOutputSets(DataSetList& DSL, std::string const& dsname,
                                   DataFile* outfile)
{
  static const DataSet::DataType type = DataSet::DOUBLE;
  dsout_tx_ = DSL.AddSet(type, MetaData(dsname, "TX"));
  if (CheckSet(dsout_tx_, dsname, "TX")) return 1;
  dsout_ty_ = DSL.AddSet(type, MetaData(dsname, "TY"));
  if (CheckSet(dsout_ty_, dsname, "TY")) return 1;
  dsout_tz_ = DSL.AddSet(type, MetaData(dsname, "TZ"));
  if (CheckSet(dsout_tz_, dsname, "TZ")) return 1;
  dsout_t_ = DSL.AddSet(type, MetaData(dsname, "T"));
  if (CheckSet(dsout_t_, dsname, "T")) return 1;

  if (outfile != 0) {
    outfile->AddDataSet( dsout_tx_ );
    outfile->AddDataSet( dsout_ty_ );
    outfile->AddDataSet( dsout_tz_ );
    outfile->AddDataSet( dsout_t_ );
  }
  return 0;
}

/** Initialize action. */
Action::RetType Action_Rotate::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  double xrot = 0.0, yrot = 0.0, zrot = 0.0;
  DataFile* outfile = 0;
  std::string output_setname;
  std::string dsname = actionArgs.GetStringKey("usedata");
  std::string calcfrom = actionArgs.GetStringKey("calcfrom");
  std::string axis = actionArgs.GetStringKey("axis0");

  if (!dsname.empty() && !calcfrom.empty() && !axis.empty()) {
    mprinterr("Error: Cannot combine any of 'axis0', 'usedata', and 'calcfrom'\n");
    return Action::ERR;
  }

  if (!dsname.empty()) {
    inverse_ = actionArgs.hasKey("inverse");
    // Check if DataSet exists
    if ( Get3x3Set(init.DSL(), dsname) )
      return Action::ERR;
    mode_ = DATASET;
  } else if (!calcfrom.empty()) {
    output_setname = actionArgs.GetStringKey("name");
    outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
    // Check if DataSet exists
    if ( Get3x3Set(init.DSL(), calcfrom) )
      return Action::ERR;
    mode_ = CALC;
  } else if (!axis.empty()) {
    // Get axis definition
    if (axis0_.SetMaskString( axis )) return Action::ERR;
    axis = actionArgs.GetStringKey("axis1");
    if (axis.empty()) {
      mprinterr("Error: 'axis1' must be specified if 'axis0' is.\n");
      return Action::ERR;
    }
    if (axis1_.SetMaskString( axis )) return Action::ERR;
    delta_ = actionArgs.getNextDouble(0.0);
    if ( !(delta_ > 0.0) && !(delta_ < 0.0) ) {
      mprinterr("Error: Must specify non-zero rotation.\n");
      return Action::ERR;
    }
    mode_ = AXIS;
  } else {
    // Calc rotation matrix
    xrot = actionArgs.getKeyDouble("x",0.0);
    yrot = actionArgs.getKeyDouble("y",0.0);
    zrot = actionArgs.getKeyDouble("z",0.0);
    RotMatrix_.CalcRotationMatrix( xrot * Constants::DEGRAD, 
                                   yrot * Constants::DEGRAD, 
                                   zrot * Constants::DEGRAD );
    if (debugIn > 0)
      RotMatrix_.Print("Rotation matrix:");
  }
  // Get mask
  if (mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;

  // Set up output data sets
  if (mode_ == CALC) {
    if (output_setname.empty())
      output_setname = init.DSL().GenerateDefaultName("ROTATE");
    if (SetupOutputSets( init.DSL(), output_setname, outfile)) return Action::ERR;
  }

  mprintf("    ROTATE: Rotating atoms in mask %s\n", mask_.MaskString());
  switch (mode_) {
    case ROTATE:
      mprintf("\t%f degrees around X, %f degrees around Y, %f degrees around Z\n",
              xrot, yrot, zrot);
      break;
    case DATASET:
      mprintf("\tUsing rotation matrices from set '%s'\n", rmatrices_->legend());
      if (inverse_) mprintf("\tPerforming inverse rotation.\n");
      break;
    case CALC:
      mprintf("\tCalculating rotations (in degrees) from rotation matrices in set '%s'\n",
              rmatrices_->legend());
      mprintf("\tOutput sets name: %s\n", output_setname.c_str());
      if (outfile != 0)
        mprintf("\tSets written to '%s'\n", outfile->DataFilename().full());
      break;
    case AXIS:
      mprintf("\t%f degrees around axis defined by '%s' and '%s'\n",
              delta_, axis0_.MaskString(), axis1_.MaskString());
      delta_ *= Constants::DEGRAD;
      break;
  }
  return Action::OK;
}

/** Set up action. */
Action::RetType Action_Rotate::Setup(ActionSetup& setup) {
  if ( setup.Top().SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }
  all_atoms_selected_ = (mask_.Nselected() == setup.Top().Natom());
  if (setup.CoordInfo().HasBox()) {
    if (all_atoms_selected_)
      mprintf("\tAll atoms selected for rotation. Rotating unit cell vectors as well.\n");
    else
      mprintf("\tNot all atoms selected for rotation. Not rotating unit cell vectors.\n");
  }
  if (mode_ == AXIS) {
    if ( setup.Top().SetupIntegerMask( axis0_ ) ||
         setup.Top().SetupIntegerMask( axis1_ ) )
      return Action::ERR;
    axis0_.MaskInfo();
    axis1_.MaskInfo();
    if (axis0_.None() || axis1_.None()) {
      mprintf("Warning: Not enough atoms selected to define axis.\n");
      return Action::SKIP;
    }
  }
  return Action::OK;
}

/** Do action. */
Action::RetType Action_Rotate::DoAction(int frameNum, ActionFrame& frm) {
  Action::RetType ret = Action::MODIFY_COORDS;
  if (mode_ == ROTATE) {
    // Rotate coordinates around X Y and Z axes
    frm.ModifyFrm().Rotate(RotMatrix_, mask_);
    if (all_atoms_selected_)
      frm.ModifyFrm().ModifyBox().RotateUcell( RotMatrix_ );
  } else if (mode_ == DATASET) {
    // Rotate coordinates using rotation matrices in DataSet
    if (frm.TrajoutNum() >= (int)rmatrices_->Size()) {
      mprintf("Warning: Frame %i out of range for set '%s'\n",
              frm.TrajoutNum()+1, rmatrices_->legend());
      return Action::ERR;
    }
    if (inverse_) {
      frm.ModifyFrm().InverseRotate((*rmatrices_)[frm.TrajoutNum()], mask_);
      if (all_atoms_selected_)
        frm.ModifyFrm().ModifyBox().InverseRotateUcell( (*rmatrices_)[frm.TrajoutNum()] );
    } else {
      frm.ModifyFrm().Rotate((*rmatrices_)[frm.TrajoutNum()], mask_);
      if (all_atoms_selected_)
        frm.ModifyFrm().ModifyBox().RotateUcell( (*rmatrices_)[frm.TrajoutNum()] );
    }
  } else if (mode_ == CALC) {
    // Calculate rotations around X Y and Z axes from rotation matrices in DataSet
    if (frm.TrajoutNum() >= (int)rmatrices_->Size()) {
      mprintf("Warning: Frame %i out of range for set '%s'\n",
              frm.TrajoutNum()+1, rmatrices_->legend());
      return Action::ERR;
    }
    Matrix_3x3 const& RM = (*rmatrices_)[frm.TrajoutNum()];
    double tx, ty, tz;
    RM.RotationAngles(tx, ty, tz);
    tx *= Constants::RADDEG;
    ty *= Constants::RADDEG;
    tz *= Constants::RADDEG;
    mprintf("DEBUG: tx= %g, ty= %g, tz= %g\n", tx, ty, tz);
    dsout_tx_->Add( frameNum, &tx );
    dsout_ty_->Add( frameNum, &ty );
    dsout_tz_->Add( frameNum, &tz );
    double theta = RM.RotationAngle() * Constants::RADDEG;
    dsout_t_->Add( frameNum, &theta );
    ret = Action::OK;
  } else if (mode_ == AXIS) {
    // Rotate around a user defined axis
    Vec3 a0 = frm.Frm().VCenterOfMass(axis0_);
    Vec3 axisOfRotation = frm.ModifyFrm().SetAxisOfRotation( a0,
                                                             frm.Frm().VCenterOfMass(axis1_) );
    RotMatrix_.CalcRotationMatrix(axisOfRotation, delta_);
    frm.ModifyFrm().Rotate(RotMatrix_, mask_);
    if (all_atoms_selected_)
      frm.ModifyFrm().ModifyBox().RotateUcell( RotMatrix_ );
    // SetAxisOfRotation moves a0 to center; move back.
    frm.ModifyFrm().Translate( a0 ); 
  } else {
    mprinterr("Internal Error: Unhandled rotation mode.\n");
    return Action::ERR;
  }
  
  return ret;
}
