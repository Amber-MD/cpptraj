// Action_Pucker
#include "Action_Pucker.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "TorsionRoutines.h"
#include "DataSet_double.h"

// CONSTRUCTOR
Action_Pucker::Action_Pucker() :
  puck_(0),
  puckerMethod_(ALTONA),
  amplitude_(false),
  useMass_(true),
  range360_(false),
  offset_(0.0)
{ } 

void Action_Pucker::Help() {
  mprintf("\t[<name>] <mask1> <mask2> <mask3> <mask4> <mask5> out <filename>\n");
  mprintf("\t[range360] [amplitude] [altona | cremer] [offset <offset>]\n");
  mprintf("\tCalculate pucker of atoms in masks 1-5.\n");
}

// Action_Pucker::init()
Action::RetType Action_Pucker::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  if      (actionArgs.hasKey("altona")) puckerMethod_=ALTONA;
  else if (actionArgs.hasKey("cremer")) puckerMethod_=CREMER;
  amplitude_ = actionArgs.hasKey("amplitude");
  offset_ = actionArgs.getKeyDouble("offset",0.0);
  range360_ = actionArgs.hasKey("range360");
  DataSet::scalarType stype = DataSet::UNDEFINED;
  std::string stypename = actionArgs.GetStringKey("type");
  if ( stypename == "pucker" ) stype = DataSet::PUCKER;

  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  std::string mask2 = actionArgs.GetMaskNext();
  std::string mask3 = actionArgs.GetMaskNext();
  std::string mask4 = actionArgs.GetMaskNext();
  std::string mask5 = actionArgs.GetMaskNext();
  if (mask1.empty() || mask2.empty() || mask3.empty() || mask4.empty() || mask5.empty()) {
    mprinterr("Error: pucker: Requires 5 masks\n");
    return Action::ERR;
  }
  M1_.SetMaskString(mask1);
  M2_.SetMaskString(mask2);
  M3_.SetMaskString(mask3);
  M4_.SetMaskString(mask4);
  M5_.SetMaskString(mask5);

  // Setup dataset
  puck_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"Pucker");
  if (puck_ == 0) return Action::ERR;
  puck_->SetScalar( DataSet::M_PUCKER, stype );
  // Add dataset to datafile list
  if (outfile != 0) outfile->AddSet( puck_ );

  mprintf("    PUCKER: [%s]-[%s]-[%s]-[%s]-[%s]\n", M1_.MaskString(),M2_.MaskString(),
          M3_.MaskString(), M4_.MaskString(), M5_.MaskString());
  if (puckerMethod_==ALTONA) 
    mprintf("            Using Altona & Sundaralingam method.\n");
  else if (puckerMethod_==CREMER)
    mprintf("            Using Cremer & Pople method.\n");
  if (outfile != 0) 
    mprintf("            Data will be written to %s\n", outfile->DataFilename().base());
  if (amplitude_)
    mprintf("            Amplitudes will be stored instead of psuedorotation.\n");
  if (offset_!=0)
    mprintf("            Offset: %lf will be added to values.\n");
  if (range360_)
    mprintf("              Output range is 0 to 360 degrees.\n");
  else
    mprintf("              Output range is -180 to 180 degrees.\n");

  return Action::OK;
}

// Action_Pucker::setup
Action::RetType Action_Pucker::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask( M1_ ) ) return Action::ERR;
  if ( currentParm->SetupIntegerMask( M2_ ) ) return Action::ERR;
  if ( currentParm->SetupIntegerMask( M3_ ) ) return Action::ERR;
  if ( currentParm->SetupIntegerMask( M4_ ) ) return Action::ERR;
  if ( currentParm->SetupIntegerMask( M5_ ) ) return Action::ERR;
  mprintf("\t");
  M1_.BriefMaskInfo();
  M2_.BriefMaskInfo();
  M3_.BriefMaskInfo();
  M4_.BriefMaskInfo();
  M5_.BriefMaskInfo();
  mprintf("\n");

  if ( M1_.None() || M2_.None() || M3_.None() || M4_.None() || M5_.None() ) {
    mprintf("Warning: pucker: One or more masks have no atoms.\n");
    return Action::ERR;
  }

  return Action::OK;  
}

// Action_Pucker::action()
Action::RetType Action_Pucker::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Vec3 a1, a2, a3, a4, a5;
  double pval, amp;

  if (useMass_) {
    a1 = currentFrame->VCenterOfMass( M1_ );  
    a2 = currentFrame->VCenterOfMass( M2_ );  
    a3 = currentFrame->VCenterOfMass( M3_ );  
    a4 = currentFrame->VCenterOfMass( M4_ );  
    a5 = currentFrame->VCenterOfMass( M5_ );
  } else {
    a1 = currentFrame->VGeometricCenter( M1_ );
    a2 = currentFrame->VGeometricCenter( M2_ );
    a3 = currentFrame->VGeometricCenter( M3_ );
    a4 = currentFrame->VGeometricCenter( M4_ );
    a5 = currentFrame->VGeometricCenter( M5_ );
  }

  switch (puckerMethod_) {
    case ALTONA: 
      pval = Pucker_AS( a1.Dptr(), a2.Dptr(), a3.Dptr(), a4.Dptr(), a5.Dptr(), &amp );
      break;
    case CREMER:
      pval = Pucker_CP( a1.Dptr(), a2.Dptr(), a3.Dptr(), a4.Dptr(), a5.Dptr(), &amp );
      break;
  }
  if ( amplitude_ )
    pval = amp;
  
  pval *= RADDEG;

  puck_->Add(frameNum, &pval);

  return Action::OK;
} 

void Action_Pucker::Print() {
  double puckermin, puckermax;
  if (range360_) {
    puckermax =  360.0;
    puckermin =    0.0;
  } else {
    puckermax =  180.0;
    puckermin = -180.0;
  }
  // Deal with offset and wrap values
  DataSet_double* ds = (DataSet_double*)puck_;
  for (int i = 0; i < ds->Size(); i++) {
    (*ds)[i] += offset_;
    if ( (*ds)[i] > puckermax )
      (*ds)[i] -= 360.0;
    else if ( (*ds)[i] < puckermin )
      (*ds)[i] += 360.0;
  }
}
