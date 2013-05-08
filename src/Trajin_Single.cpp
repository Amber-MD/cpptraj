#include "Trajin_Single.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists

// CONSTRUCTOR
Trajin_Single::Trajin_Single() :
  trajio_(0),
  velio_(0),
  trajIsOpen_(false)
{}

// DESTRUCTOR
Trajin_Single::~Trajin_Single() {
  if (trajio_!=0) {
    if (trajIsOpen_) EndTraj();
    delete trajio_;
  }
  if (velio_!=0) delete velio_;
}

int Trajin_Single::SetupTrajRead(std::string const& tnameIn, ArgList *argIn, Topology *tparmIn)
{
  return SetupTrajRead(tnameIn, argIn, tparmIn, true);
}

// Trajin_Single::SetupTrajRead()
int Trajin_Single::SetupTrajRead(std::string const& tnameIn, ArgList *argIn, 
                                 Topology *tparmIn, bool checkBox) 
{
  // Require a filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: Trajin_Single: No filename given.\n");
    return 1;
  }
  // Check and set associated parm file
  if ( SetTrajParm( tparmIn ) ) return 1;
  // Check that file exists
  // FIXME: Is this necessary? Just let cpptrajfile check this?
  if (!fileExists(tnameIn.c_str())) {
    mprinterr("Error: File %s does not exist.\n",tnameIn.c_str());
    return 1;
  }
  // Detect file format
  if ( (trajio_ = DetectFormat( tnameIn )) == 0 ) {
    mprinterr("Error: Could not determine trajectory %s format.\n", tnameIn.c_str());
    return 1;
  }
  trajio_->SetDebug( debug_ );
  // Set trajectory filename
  SetTrajFileName( tnameIn, true );
  // Set up the format for reading and get the number of frames.
  if (SetupTrajIO( tnameIn, *trajio_, argIn )) return 1;
  // Check how many frames will actually be read
  if (setupFrameInfo() == 0) return 1;
  // Check traj box info against parm box info
  if (checkBox) { 
    Box parmBox = tparmIn->ParmBox();
    if (CheckBoxInfo(tparmIn->c_str(), parmBox, trajio_->TrajBox())) return 1;
    tparmIn->SetBox( parmBox );
  }
  // Check if a separate mdvel file will be read
  if (argIn!=0 && argIn->Contains("mdvel")) {
    std::string mdvelname = argIn->GetStringKey("mdvel");
    if (mdvelname.empty()) {
      mprinterr("Error: mdvel: Usage 'mdvel <velocity filename>'\n");
      return 1;
    }
    // Detect mdvel format
    if ( (velio_ = DetectFormat( mdvelname )) == 0 ) {
      mprinterr("Error: Could not set up velocity file %s for reading.\n",mdvelname.c_str());
      return 1;
    }
    velio_->SetDebug( debug_ );
    // Set up the format for reading mdvel, get # of mdvel frames
    int vel_frames = velio_->setupTrajin(mdvelname, TrajParm());
    if (vel_frames != TotalFrames()) {
      mprinterr("Error: velocity file %s frames (%i) != traj file frames (%i)\n",
                mdvelname.c_str(), vel_frames, TotalFrames());
      return 1;
    }
    if (trajio_->IsSeekable() && !velio_->IsSeekable()) {
      mprinterr("Error: traj %s is seekable but velocity file %s is not.\n",
                TrajFilename().full(), mdvelname.c_str());
      return 1;
    }
  }
  return 0;
}

// Trajin_Single::BeginTraj()
int Trajin_Single::BeginTraj(bool showProgress) {
  // Open the trajectory
  if (trajio_->openTrajin()) {
    mprinterr("Error: Trajin_Single::BeginTraj: Could not open %s\n",TrajFilename().base());
    return 1;
  }
  // Open mdvel file if present
  if (velio_!=0 && velio_->openTrajin()) {
    mprinterr("Error: Could not open mdvel file.\n");
    return 1;
  }
  // Set progress bar, start and offset.
  PrepareForRead( showProgress, trajio_->IsSeekable() );
  trajIsOpen_ = true;
  return 0;
}

// Trajin_Single::EndTraj()
void Trajin_Single::EndTraj() {
  if (trajIsOpen_) {
    trajio_->closeTraj();
    if (velio_!=0) velio_->closeTraj();
    trajIsOpen_ = false;
  }
}

// Trajin_Single::GetNextFrame()
int Trajin_Single::GetNextFrame( Frame& frameIn ) {
  // If the current frame is out of range, exit
  if ( CheckFinished() ) return 0;
  bool tgtFrameFound = false;
  while ( !tgtFrameFound ) {
    if (trajio_->readFrame(CurrentFrame(), frameIn.xAddress(), frameIn.vAddress(),
                           frameIn.bAddress(), frameIn.tAddress()))
      return 0;
    if (velio_ != 0 && velio_->readVelocity(CurrentFrame(), frameIn.vAddress()))
      return 0;
    // Check if coords in frame are valid.
    if (frameIn.CheckCoordsInvalid())
      mprintf("Warning: Frame %i coords 1 & 2 overlap at origin; may be corrupt.\n",
              CurrentFrame()+1);
    //printf("DEBUG:\t%s:  current=%i  target=%i\n",trajName,currentFrame,targetSet);
    tgtFrameFound = ProcessFrame();
  }
  return 1;
}

// Trajin_Single::PrintInfo()
void Trajin_Single::PrintInfo(int showExtended) {
  mprintf("[%s] ",TrajFilename().base());
  trajio_->Info();
  mprintf(", Parm %s",TrajParm()->c_str());
  if (trajio_->HasBox()) mprintf(" (%s box)", trajio_->TrajBox().TypeName());
  if (showExtended==1) PrintFrameInfo(); 
  if (debug_>0)
    mprintf(", %i atoms, Box %i, seekable %i",TrajParm()->Natom(),(int)trajio_->HasBox(),
            (int)trajio_->IsSeekable());
  mprintf("\n");
  if (velio_!=0) {
    mprintf("\tMDVEL: ");
    velio_->Info();
    mprintf("\n");
  }
}

// Trajin_Single::HasVelocity()
bool Trajin_Single::HasVelocity() {
  if (trajio_!=0) {
    if (velio_ == 0)
      return trajio_->HasV();
    else
      return true;
  }
  return false;
}

