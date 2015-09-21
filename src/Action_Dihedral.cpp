// Action_Dihedral
#include "Action_Dihedral.h"
#include "DataSet_double.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "TorsionRoutines.h"

// CONSTRUCTOR
Action_Dihedral::Action_Dihedral() :
  dih_(0),
  minTorsion_(-180.0),
  useMass_(false)
{ } 

void Action_Dihedral::Help() {
  mprintf("\t[<name>] <mask1> <mask2> <mask3> <mask4> [out filename] [mass]\n"
          "\t[type {alpha|beta|gamma|delta|epsilon|zeta|chi|c2p|h1p|phi|psi|pchi}]\n"
          "\t[range360] [idx <index>]\n"
          "  Calculate dihedral angle for atoms in masks 1-4.\n");
}

// Action_Dihedral::Init()
Action::RetType Action_Dihedral::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  useMass_ = actionArgs.hasKey("mass");
  if (actionArgs.hasKey("range360"))
    minTorsion_ = 0.0;
  else
    minTorsion_ = -180.0;
  int dsidx = actionArgs.getKeyInt("idx", -1);
  MetaData::scalarType stype = MetaData::UNDEFINED;
  std::string stypename = actionArgs.GetStringKey("type");
  if (!stypename.empty()) {
    stype = MetaData::TypeFromKeyword( stypename, MetaData::M_TORSION );
    if (stype == MetaData::UNDEFINED) {
      mprinterr("Error: Invalid torsion type keyword '%s'\n", stypename.c_str());
      return Action::ERR;
    }
  }

  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  std::string mask2 = actionArgs.GetMaskNext();
  std::string mask3 = actionArgs.GetMaskNext();
  std::string mask4 = actionArgs.GetMaskNext();
  if (mask1.empty() || mask2.empty() || mask3.empty() || mask4.empty()) {
    mprinterr("Error: dihedral: Requires 4 masks\n");
    return Action::ERR;
  }
  M1_.SetMaskString(mask1);
  M2_.SetMaskString(mask2);
  M3_.SetMaskString(mask3);
  M4_.SetMaskString(mask4);

  // Setup dataset
  dih_ = DSL->AddSet(DataSet::DOUBLE, MetaData(actionArgs.GetStringNext(), dsidx,
                                               MetaData::M_TORSION, stype),"Dih");
  if (dih_==0) return Action::ERR;
  // Add dataset to datafile list
  if (outfile != 0) outfile->AddDataSet( dih_ );

  mprintf("    DIHEDRAL: [%s]-[%s]-[%s]-[%s]\n", M1_.MaskString(), 
          M2_.MaskString(), M3_.MaskString(), M4_.MaskString());
  if (useMass_)
    mprintf("              Using center of mass of atoms in masks.\n");
  if (minTorsion_ > -180.0)
    mprintf("              Output range is 0 to 360 degrees.\n");
  else
    mprintf("              Output range is -180 to 180 degrees.\n");

  return Action::OK;
}

// Action_Dihedral::Setup
Action::RetType Action_Dihedral::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupIntegerMask(M1_)) return Action::ERR;
  if (currentParm->SetupIntegerMask(M2_)) return Action::ERR;
  if (currentParm->SetupIntegerMask(M3_)) return Action::ERR;
  if (currentParm->SetupIntegerMask(M4_)) return Action::ERR;
  mprintf("\t");
  M1_.BriefMaskInfo();
  M2_.BriefMaskInfo();
  M3_.BriefMaskInfo();
  M4_.BriefMaskInfo();
  mprintf("\n");
  if ( M1_.None() || M2_.None() || M3_.None() || M4_.None() ) {
    mprintf("Warning: dihedral: One or more masks have no atoms.\n");
    return Action::ERR;
  }

  return Action::OK;  
}

// Action_Dihedral::DoAction()
Action::RetType Action_Dihedral::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Vec3 a1, a2, a3, a4;

  if (useMass_) {
    a1 = currentFrame->VCenterOfMass( M1_ );
    a2 = currentFrame->VCenterOfMass( M2_ );
    a3 = currentFrame->VCenterOfMass( M3_ );
    a4 = currentFrame->VCenterOfMass( M4_ );
  } else {
    a1 = currentFrame->VGeometricCenter( M1_ );
    a2 = currentFrame->VGeometricCenter( M2_ );
    a3 = currentFrame->VGeometricCenter( M3_ );
    a4 = currentFrame->VGeometricCenter( M4_ );
  }
  double torsion = Torsion(a1.Dptr(), a2.Dptr(), a3.Dptr(), a4.Dptr());

  torsion *= Constants::RADDEG;

  if (torsion < minTorsion_)
    torsion += 360.0; 

  dih_->Add(frameNum, &torsion);

  return Action::OK;
}
