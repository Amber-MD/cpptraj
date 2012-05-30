// Action_Pucker
#include "Action_Pucker.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG

// CONSTRUCTOR
Action_Pucker::Action_Pucker() :
  puck_(NULL),
  puckerMethod_(ALTONA),
  amplitude_(false),
  offset_(0),
  puckermin_( -180.0),
  puckermax_( 180.0)
{
  //fprintf(stderr,"Pucker Con\n");
  useMass_ = true;
} 

// Action_Pucker::init()
/** Expected call: pucker <name> <mask1> <mask2> <mask3> <mask4> <mask5> out <filename>
  *                [range360] [amplitude] [altona | cremer] [offset <offset>]
  */
int Action_Pucker::init() {
  // Get keywords
  char* puckerFile = actionArgs.getKeyString("out",NULL);
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
  char* mask1 = actionArgs.getNextMask();
  char* mask2 = actionArgs.getNextMask();
  char* mask3 = actionArgs.getNextMask();
  char* mask4 = actionArgs.getNextMask();
  char* mask5 = actionArgs.getNextMask();
  if (mask1==NULL || mask2==NULL || mask3==NULL || mask4==NULL || mask5==NULL) {
    mprinterr("Error: pucker: Requires 5 masks\n");
    return 1;
  }
  M1_.SetMaskString(mask1);
  M2_.SetMaskString(mask2);
  M3_.SetMaskString(mask3);
  M4_.SetMaskString(mask4);
  M5_.SetMaskString(mask5);

  // Setup dataset
  puck_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"Pucker");
  if (puck_==NULL) return 1;
  puck_->SetScalar( DataSet::M_PUCKER, stype );
  // Add dataset to datafile list
  DFL->Add(puckerFile, puck_);

  //dih->Info();
  mprintf("    PUCKER: [%s]-[%s]-[%s]-[%s]-[%s]\n", M1_.MaskString(),M2_.MaskString(),
          M3_.MaskString(), M4_.MaskString(), M5_.MaskString());
  if (puckerMethod_==ALTONA) 
    mprintf("            Using Altona & Sundaralingam method.\n");
  else if (puckerMethod_==CREMER)
    mprintf("            Using Cremer & Pople method.\n");
  if (puckerFile!=NULL) 
    mprintf("            Data will be written to %s\n",puckerFile);
  if (amplitude_)
    mprintf("            Amplitudes will be stored instead of psuedorotation.\n");
  if (offset_!=0)
    mprintf("            Offset: %lf will be added to values.\n");
  mprintf  ("            Values will range from %.1lf to %.1lf\n",puckermin_,puckermax_);

  return 0;
}

// Action_Pucker::setup
int Action_Pucker::setup() {
  if ( currentParm->SetupIntegerMask( M1_ ) ) return 1;
  if ( currentParm->SetupIntegerMask( M2_ ) ) return 1;
  if ( currentParm->SetupIntegerMask( M3_ ) ) return 1;
  if ( currentParm->SetupIntegerMask( M4_ ) ) return 1;
  if ( currentParm->SetupIntegerMask( M5_ ) ) return 1;
  mprintf("\t%s (%i atoms)\n",M1_.MaskString(),M1_.Nselected());
  mprintf("\t%s (%i atoms)\n",M2_.MaskString(),M2_.Nselected());
  mprintf("\t%s (%i atoms)\n",M3_.MaskString(),M3_.Nselected());
  mprintf("\t%s (%i atoms)\n",M4_.MaskString(),M4_.Nselected());
  mprintf("\t%s (%i atoms)\n",M5_.MaskString(),M5_.Nselected());

  if ( M1_.None() || M2_.None() || M3_.None() || M4_.None() || M5_.None() ) {
    mprintf("Warning: pucker: One or more masks have no atoms.\n");
    return 1;
  }

  return 0;  
}

// Action_Pucker::action()
int Action_Pucker::action() {
  double pval = currentFrame->PUCKER(&M1_,&M2_,&M3_,&M4_,&M5_,puckerMethod_,amplitude_,useMass_);
  pval *= RADDEG;

  // Deal with offset
  pval += offset_;

  // Wrap values > puckermax or < puckermin
  if      (pval > puckermax_) pval -= 360.0;
  else if (pval < puckermin_) pval += 360.0;

  puck_->Add(frameNum, &pval);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,pval);
  
  return 0;
} 

