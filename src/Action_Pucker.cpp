// Action_Pucker
#include "Action_Pucker.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "TorsionRoutines.h"

// CONSTRUCTOR
Action_Pucker::Action_Pucker() :
  puck_(NULL),
  puckerMethod_(ALTONA),
  amplitude_(false),
  useMass_(true),
  offset_(0),
  puckermin_( -180.0),
  puckermax_( 180.0)
{ } 

void Action_Pucker::Help() {
  mprintf("pucker [<name>] <mask1> <mask2> <mask3> <mask4> <mask5> out <filename>\n");
  mprintf("       [range360] [amplitude] [altona | cremer] [offset <offset>]\n");
}

// Action_Pucker::init()
Action::RetType Action_Pucker::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  std::string puckerFile = actionArgs.GetStringKey("out");
  if      (actionArgs.hasKey("altona")) puckerMethod_=ALTONA;
  else if (actionArgs.hasKey("cremer")) puckerMethod_=CREMER;
  amplitude_ = actionArgs.hasKey("amplitude");
  offset_ = actionArgs.getKeyDouble("offset",0.0);
  if (actionArgs.hasKey("range360")) {
    puckermax_=360.0;
    puckermin_=0.0;
  }
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
  if (puck_==NULL) return Action::ERR;
  puck_->SetScalar( DataSet::M_PUCKER, stype );
  // Add dataset to datafile list
  DFL->AddSetToFile(puckerFile, puck_, actionArgs);

  //dih->Info();
  mprintf("    PUCKER: [%s]-[%s]-[%s]-[%s]-[%s]\n", M1_.MaskString(),M2_.MaskString(),
          M3_.MaskString(), M4_.MaskString(), M5_.MaskString());
  if (puckerMethod_==ALTONA) 
    mprintf("            Using Altona & Sundaralingam method.\n");
  else if (puckerMethod_==CREMER)
    mprintf("            Using Cremer & Pople method.\n");
  if (!puckerFile.empty()) 
    mprintf("            Data will be written to %s\n",puckerFile.c_str());
  if (amplitude_)
    mprintf("            Amplitudes will be stored instead of psuedorotation.\n");
  if (offset_!=0)
    mprintf("            Offset: %lf will be added to values.\n");
  mprintf  ("            Values will range from %.1lf to %.1lf\n",puckermin_,puckermax_);

  return Action::OK;
}

// Action_Pucker::setup
Action::RetType Action_Pucker::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask( M1_ ) ) return Action::ERR;
  if ( currentParm->SetupIntegerMask( M2_ ) ) return Action::ERR;
  if ( currentParm->SetupIntegerMask( M3_ ) ) return Action::ERR;
  if ( currentParm->SetupIntegerMask( M4_ ) ) return Action::ERR;
  if ( currentParm->SetupIntegerMask( M5_ ) ) return Action::ERR;
  M1_.MaskInfo();
  M2_.MaskInfo();
  M3_.MaskInfo();
  M4_.MaskInfo();
  M5_.MaskInfo();

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

  // Deal with offset
  pval += offset_;

  // Wrap values > puckermax or < puckermin
  if      (pval > puckermax_) pval -= 360.0;
  else if (pval < puckermin_) pval += 360.0;

  puck_->Add(frameNum, &pval);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,pval);
  
  return Action::OK;
} 

