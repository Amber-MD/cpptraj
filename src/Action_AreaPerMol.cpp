#include "Action_AreaPerMol.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_AreaPerMol::Action_AreaPerMol() :
  area_per_mol_(0),
  Nmols_(-1.0),
  Nlayers_(1.0),
  areaType_(XY)
{}

void Action_AreaPerMol::Help() const {
  mprintf("\t[<name>] {[<mask1>] [nlayers <#>] | nmols <#>} [out <filename>] [{xy | xz | yz}]\n"
          "  Calculate the specified area per molecule for molecules in <mask1>.\n");
}

static const char* APMSTRING[] = {"XY", "XZ", "YZ"};

// Action_AreaPerMol::Init()
Action::RetType Action_AreaPerMol::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  if (actionArgs.hasKey("xy")) areaType_ = XY;
  else if (actionArgs.hasKey("xz")) areaType_ = XZ;
  else if (actionArgs.hasKey("yz")) areaType_ = YZ;
  else areaType_ = XY;

  Nmols_ = (double)actionArgs.getKeyInt("nmols", -1);

  // Get Masks
  if (Nmols_ < 0.0) {
    Nlayers_ = (double)actionArgs.getKeyInt("nlayers", 1);
    if (Nlayers_ < 1.0) {
      mprinterr("Error: Number of layers must be > 0\n");
      return Action::ERR;
    }
    if (Mask1_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  }

  // DataSet
  area_per_mol_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"APM");
  if (area_per_mol_==0) return Action::ERR;
  // Add DataSet to DataFileList
  if (outfile != 0) outfile->AddDataSet( area_per_mol_ );

  mprintf("    AREAPERMOL: Calculating %s area per molecule", APMSTRING[areaType_]);
  if (Mask1_.MaskStringSet())
    mprintf(" using mask '%s', %.0f layers.\n", Mask1_.MaskString(), Nlayers_);
  else
    mprintf(" for %.0f mols\n", Nmols_);

  return Action::OK;
}

// Action_AreaPerMol::Setup()
/** Set angle up for this parmtop. Get masks etc.
  */
Action::RetType Action_AreaPerMol::Setup(ActionSetup& setup) {
  // Needs box info
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprintf("Warning: No box information for '%s', cannot calculate area.\n",
            setup.Top().c_str());
    return Action::SKIP;
  }
  // Probably will not work for non-orthorhombic cells
  if (!setup.CoordInfo().TrajBox().Is_X_Aligned_Ortho())
    mprintf("Warning: Box is not X-aligned orthorhombic, calculated area may not be correct.\n");
  // Determine how many molecules are selected
  if (Mask1_.MaskStringSet()) {
    if (setup.Top().SetupCharMask(Mask1_)) return Action::ERR;
    if (Mask1_.None()) {
      mprinterr("Warning: Mask '%s' selects no atoms.\n", Mask1_.MaskString());
      return Action::SKIP;
    }
    Nmols_ = 0.0;
    for (Topology::mol_iterator mol = setup.Top().MolStart();
                                mol != setup.Top().MolEnd(); ++mol)
    {
      if (Mask1_.AtomsInCharMask(mol->MolUnit()))
        Nmols_ += 1.0;
    }
    mprintf("\tMask '%s' selects %.0f molecules.\n", Mask1_.MaskString(), Nmols_);
    if (Nmols_ < 1.0) return Action::SKIP;
    Nmols_ /= Nlayers_;
    mprintf("\tArea per %.0f molecules (%0.f layers) will be determined.\n", Nmols_, Nlayers_);
  } else
    mprintf("\tArea per %.0f molecules will be determined.\n", Nmols_);
  return Action::OK;
}

// Action_AreaPerMol::DoAction()
Action::RetType Action_AreaPerMol::DoAction(int frameNum, ActionFrame& frm) {
  double area;
  if (areaType_ == XY)
    area = frm.Frm().BoxCrd().Param(Box::X) * frm.Frm().BoxCrd().Param(Box::Y);
  else if (areaType_ == XZ) 
    area = frm.Frm().BoxCrd().Param(Box::X) * frm.Frm().BoxCrd().Param(Box::Z);
  else // if areaType_ == YZ
    area = frm.Frm().BoxCrd().Param(Box::Y) * frm.Frm().BoxCrd().Param(Box::Z);

  area = area / Nmols_;

  area_per_mol_->Add(frameNum, &area);

  return Action::OK;
}
