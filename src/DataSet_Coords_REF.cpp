#include "DataSet_Coords_REF.h"
#include "CpptrajStdio.h"
#include "Trajin_Single.h"
#include "ArgList.h"

/** Print info about REF coords. */
void DataSet_Coords_REF::Info() const {
  //if (!tag_.empty())
  //  mprintf(" %s", tag_.c_str());
  if (!Meta().Fname().empty() && Meta().Name() != Meta().Fname().Full())
    mprintf(" '%s'", Meta().Fname().full());
  CommonInfo();
}

/** Set up REF coords with given Topology and coordinate info. */
int DataSet_Coords_REF::CoordsSetup(Topology const& topIn, CoordinateInfo const& cInfoIn) {
  top_ = topIn;
  cInfo_ = cInfoIn;

  if (frame_.SetupFrameV( topIn.Atoms(), cInfoIn )) {
    mprinterr("Internal Error: DataSet_Coords_REF: Could not allocate Frame.\n");
    return 1;
  }
  return 0;
}

/** Set REF from incoming frame. */
void DataSet_Coords_REF::AddFrame(Frame const& fIn) {
  if (frame_.empty())
    frame_.SetupFrame( fIn ); // TODO always do?
  frame_.SetFrame( fIn );
}

/** Set REF from incoming frame. */
void DataSet_Coords_REF::SetCRD(int idx, Frame const& fIn) {
  // TODO warn if idx not zero?
  if (idx == 0)
    mprintf("Warning: In reference set '%s', attempting to set index which is not 1 (%i)\n", legend(), idx+1);
  if (frame_.empty())
    frame_.SetupFrame( fIn );
  frame_.SetFrame( fIn );
}

/** Load reference from file, no arguments or name. */
int DataSet_Coords_REF::LoadRefFromFile(FileName const& fname, Topology const& parmIn, int dbg)
{
  ArgList blank;
  return LoadRefFromFile(fname, "", parmIn, blank, dbg);
}

/** Load reference from file. */
int DataSet_Coords_REF::LoadRefFromFile(FileName const& fname, std::string const& nameIn,
                                        Topology const& parmIn, ArgList& argIn, int dbg)
{
  // Set up input trajectory
  Trajin_Single traj;
  traj.SetDebug( dbg );
  if ( traj.SetupTrajRead( fname, argIn, (Topology*)&parmIn ) ) { // TODO: should there be a version of Setup that takes const&? 
    mprinterr("Error: DataSet_Coords_REF: Could not set up trajectory.\n");
    return 1;
  }
  // Check number of frames to be read
  int trajFrames = traj.Traj().Counter().TotalReadFrames();
  if (trajFrames < 1) {
    mprinterr("Error: No frames could be read for reference '%s'\n", traj.Traj().Filename().full());
    return 1;
  } else if (trajFrames > 1)
    mprintf("Warning: Reference has %i frames, only reading frame %i\n",
            trajFrames, traj.Traj().Counter().Start()+1);
  // Set up REF set
  if (CoordsSetup( parmIn, traj.TrajCoordInfo() )) {
    mprinterr("Error: Could not set up reference DataSet for %s\n", traj.Traj().Filename().full());
    return 1;
  }
  // Start trajectory read
  if ( traj.BeginTraj() ) {
    mprinterr("Error: Could not open reference '%s'\n.", traj.Traj().Filename().full());
    return 1;
  }
  // Read reference frame
  traj.ReadTrajFrame( traj.Traj().Counter().Start(), frame_ );
  // Close trajectory
  traj.EndTraj();
  // Set up DataSet metadata
  MetaData md(fname, nameIn, traj.Traj().Counter().Start()+1);
  if (!traj.Title().empty())
    md.SetLegend( traj.Title() );
  if (SetMeta( md )) return 1;
  return 0;
}

/** Set reference from given COORDS. */ // TODO deprecate
int DataSet_Coords_REF::SetRefFromCoords(DataSet_Coords* CRD, std::string const& nameIn, int fnum)
{
  if (CRD==0) return 1;
  frame_ = CRD->AllocateFrame();
  CRD->GetFrame( fnum, frame_ );
  CoordsSetup( CRD->Top(), CRD->CoordsInfo() );
  std::string setname;
  if (nameIn.empty())
    setname = CRD->Meta().Name();
  else
    setname = nameIn;
  if (SetMeta( MetaData(setname, fnum+1) )) return 1;
  return 0;
}

/** Strip down reference according to the given atom mask. */
int DataSet_Coords_REF::StripRef(std::string const& maskexpr) {
  if (maskexpr.empty()) return 1;
  AtomMask stripMask( maskexpr );
  if (Top().SetupIntegerMask( stripMask )) return 1;
  return StripRef( stripMask );
}

/** Currently used by 'reference' and 'atommap' */
int DataSet_Coords_REF::StripRef(AtomMask const& stripMask) {
    Frame stripFrame( frame_, stripMask );
    Topology* stripParm = Top().modifyStateByMask( stripMask );
    if (stripParm == 0) {
      mprinterr("Error: Could not create stripped reference topology.\n");
      return 1;
    }
    stripParm->Brief("Stripped ref parm:");
    frame_ = stripFrame;
    CoordsSetup( *stripParm, cInfo_ );
    delete stripParm; // OK to free, parm has been copied by CoordsSetup.
    return 0;
}
