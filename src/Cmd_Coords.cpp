#include "Cmd_Coords.h"
#include "CpptrajStdio.h"
#include "Timer.h"
#include "ProgressBar.h"
#include "Trajin_Single.h"
#include "DataSet_Coords_TRJ.h"

void Help_CrdAction() {
  mprintf("\t<crd set> <actioncmd> [<action args>] [crdframes <start>,<stop>,<offset>]\n"
          "  Perform action <actioncmd> on COORDS data set <crd set>.\n");
}

Cmd::RetType CrdAction(CpptrajState& State, ArgList& actionargs, DataSet_Coords* CRD,
                       Action* act, TrajFrameCounter const& frameCount)
{
  Timer total_time;
  total_time.Start();
  ActionInit state(*State.DSL(), *State.DFL());
  if ( act->Init( actionargs, state, State.Debug() ) != Action::OK ) {
    delete act;
    return Cmd::ERR;
  }
  actionargs.CheckForMoreArgs();
  // Set up frame and parm for COORDS.
  ActionSetup originalSetup( CRD->TopPtr(), CRD->CoordsInfo(), CRD->Size() );
  Frame originalFrame = CRD->AllocateFrame();
  ActionFrame frm( &originalFrame );
  // Set up for this topology
  Action::RetType setup_ret = act->Setup( originalSetup );
  if ( setup_ret == Action::ERR ) {
    delete act;
    return Cmd::ERR;
  }
  // Loop over all frames in COORDS.
  ProgressBar progress( frameCount.TotalReadFrames() );
  int set = 0;
  for (int frame = frameCount.Start(); frame < frameCount.Stop();
           frame += frameCount.Offset(), ++set)
  {
    progress.Update( set );
    CRD->GetFrame( frame, originalFrame );
    Action::RetType ret = act->DoAction( set, frm );
    if (ret == Action::ERR) {
      mprinterr("Error: crdaction: Frame %i, set %i\n", frame + 1, set + 1);
      break;
    }
    // Check if frame was modified. If so, update COORDS.
    if ( ret == Action::MODIFY_COORDS ) 
      CRD->SetCRD( frame, frm.Frm() );
  }
  // Check if parm was modified. If so, update COORDS.
  if ( setup_ret == Action::MODIFY_TOPOLOGY ) {
    mprintf("Info: crdaction: Parm for %s was modified by action %s\n",
            CRD->legend(), actionargs.Command());
    CRD->CoordsSetup( originalSetup.Top(), originalSetup.CoordInfo() );
  }
  act->Print();
  State.MasterDataFileWrite();
  delete act;
  total_time.Stop();
  mprintf("TIME: Total action execution time: %.4f seconds.\n", total_time.Total());
  return Cmd::OK;
}

// -----------------------------------------------------------------------------
void Help_CrdOut() {
  mprintf("\t<crd set> <filename> [<trajout args>] [crdframes <start>,<stop>,<offset>]\n"
          "  Write COORDS data set <crd set> to trajectory file <filename>\n");
}

Cmd::RetType CrdOut(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: crdout: Specify COORDS dataset name.\n");
    return Cmd::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL()->FindCoordsSet( setname );
  if (CRD == 0) {
    mprinterr("Error: crdout: No COORDS set with name %s found.\n", setname.c_str());
    return Cmd::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  setname = argIn.GetStringNext(); // Output traj file name
  // Start, stop, offset
  TrajFrameCounter frameCount;
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  if (frameCount.CheckFrameArgs( CRD->Size(), crdarg )) return Cmd::ERR;
  frameCount.PrintInfoLine( CRD->legend() );
  Trajout_Single outtraj;
  if (outtraj.PrepareTrajWrite( setname, argIn, CRD->TopPtr(), CRD->CoordsInfo(),
                                CRD->Size(), TrajectoryFile::UNKNOWN_TRAJ))
  {
    mprinterr("Error: crdout: Could not set up output trajectory.\n");
    return Cmd::ERR;
  }
  outtraj.PrintInfo(0);
  Frame currentFrame = CRD->AllocateFrame(); 
  ProgressBar progress( frameCount.TotalReadFrames() );
  int set = 0;
  for (int frame = frameCount.Start(); frame < frameCount.Stop();
           frame += frameCount.Offset(), ++set)
  {
    progress.Update( set );
    CRD->GetFrame( frame, currentFrame );
    if ( outtraj.WriteSingle( frame, currentFrame ) ) {
      mprinterr("Error writing %s to output trajectory, frame %i.\n",
                CRD->legend(), frame + 1);
      break;
    }
  }
  return Cmd::OK;
}

// -----------------------------------------------------------------------------
void Help_LoadCrd() {
  mprintf("\t<filename> [%s] [<trajin args>] [name <name>]\n", DataSetList::TopArgs);
  mprintf("  Load trajectory <filename> as a COORDS data set named <name> (default <filename>).\n");
}

Cmd::RetType LoadCrd(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  // Get parm
  Topology* parm = State.DSL()->GetTopology( argIn );
  if (parm == 0) {
    mprinterr("Error: loadcrd: No parm files loaded.\n");
    return Cmd::ERR;
  }
  // Load trajectory
  Trajin_Single trajin;
  trajin.SetDebug( State.Debug() );
  if (trajin.SetupTrajRead(argIn.GetStringNext(), argIn, parm)) {
    mprinterr("Error: loadcrd: Could not set up input trajectory.\n");
    return Cmd::ERR;
  }
  // Create input frame
  Frame frameIn;
  frameIn.SetupFrameV(parm->Atoms(), trajin.TrajCoordInfo());
  // Set up metadata with file name and output set name
  std::string setname = argIn.GetStringKey("name");
  // NOTE: For backwards compat.
  if (setname.empty()) setname = argIn.GetStringNext();
  MetaData md( trajin.Traj().Filename(), setname, -1 );
  // Check if set already present
  DataSet_Coords* coords = 0;
  DataSet* ds = State.DSL()->FindSetOfType( setname, DataSet::COORDS );
  if (ds == 0) {
    // Create Set 
    coords = (DataSet_Coords*)State.DSL()->AddSet(DataSet::COORDS, md);
    if (coords == 0) {
      mprinterr("Error: loadcrd: Could not set up COORDS data set.\n");
      return Cmd::ERR;
    }
    coords->CoordsSetup( *parm, trajin.TrajCoordInfo() );
    mprintf("\tLoading trajectory '%s' as '%s'\n", trajin.Traj().Filename().full(),
            coords->legend());
  } else {
    // Check that set is actually coords.
    if (ds->Type() != DataSet::COORDS) {
      mprinterr("Error: Set %s present but is not of type COORDS.\n", ds->legend());
      return Cmd::ERR;
    }
    coords = (DataSet_Coords*)ds;
    // Check that topology matches. For now just check # atoms.
    if (parm->Natom() != coords->Top().Natom()) {
      mprinterr("Error: Trajectory '%s' # atoms %i does not match COORDS data set '%s' (%i)\n",
                trajin.Traj().Filename().full(), parm->Natom(),
                coords->legend(), coords->Top().Natom());
      return Cmd::ERR;
    }
    mprintf("\tAppending trajectory '%s' to COORDS data set '%s'\n", 
            trajin.Traj().Filename().full(), coords->legend());
  }
  // Read trajectory TODO progress bar
  trajin.BeginTraj();
  trajin.Traj().PrintInfoLine();
  while (trajin.GetNextFrame( frameIn ))
    coords->AddFrame( frameIn );
  trajin.EndTraj();
  return Cmd::OK;
}

// -----------------------------------------------------------------------------
void Help_LoadTraj() {
  mprintf("\tname <setname> [<filename>]\n"
          "  Create/add to TRAJ data set named <setname>. If no <filename> given, convert\n"
          "  currently loaded input trajectories to TRAJ data set; otherwise add <filename>\n"
          "  to TRAJ data set <setname>\n");
}

Cmd::RetType LoadTraj(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  // Get Keywords
  std::string setname = argIn.GetStringKey("name");
  if (setname.empty()) {
    mprinterr("Error: Must provide data set name ('name <setname>')\n");
    return Cmd::ERR;
  }
  DataSet_Coords_TRJ* trj = (DataSet_Coords_TRJ*)
                            State.DSL()->FindSetOfType(setname, DataSet::TRAJ);
  if (trj == 0)
    trj = (DataSet_Coords_TRJ*)
          State.DSL()->AddSet(DataSet::TRAJ, setname, "__DTRJ__");
  if (trj == 0) {
    mprinterr("Error: Could not set up TRAJ data set.\n");
    return Cmd::ERR;
  }
  std::string trajname = argIn.GetStringNext();
  if (trajname.empty()) {
    // Add all existing input trajectories
    if (State.InputTrajList().empty()) {
      mprinterr("Error: No input trajectories loaded.\n");
      return Cmd::ERR;
    }
    if (State.InputTrajList().Mode() != TrajinList::NORMAL) {
      mprinterr("Error: Cannot convert ensemble input trajectories to data.\n");
      return Cmd::ERR;
    }
    mprintf("\tSaving currently loaded input trajectories as data set with name '%s'\n",
            setname.c_str());
    for (TrajinList::trajin_it Trajin = State.InputTrajList().trajin_begin();
                               Trajin != State.InputTrajList().trajin_end(); ++Trajin)
      if (trj->AddInputTraj( *Trajin )) return Cmd::ERR;
    // TODO: Clear input trajectories from trajinList?
  } else {
    // Add the named trajectory
    if (trj->AddSingleTrajin( trajname, argIn, State.DSL()->GetTopology(argIn) ))
      return Cmd::ERR;
  }
  return Cmd::OK;
}

// -----------------------------------------------------------------------------
void Help_CombineCoords() {
  mprintf("\t<crd1> <crd2> ... [parmname <topname>] [crdname <crdname>]\n"
          "  Combine two or more COORDS data sets.\n");
}

Cmd::RetType CombineCoords(CpptrajState& State, ArgList& argIn, Cmd::AllocType Alloc)
{
  std::string parmname = argIn.GetStringKey("parmname");
  std::string crdname  = argIn.GetStringKey("crdname");
  // Get COORDS DataSets.
  std::vector<DataSet_Coords*> CRD;
  std::string setname = argIn.GetStringNext();
  while (!setname.empty()) {
    DataSet_Coords* ds = (DataSet_Coords*)State.DSL()->FindCoordsSet( setname );
    if (ds == 0) {
      mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
      return Cmd::ERR;
    }
    CRD.push_back( ds );
    setname = argIn.GetStringNext();
  }
  if (CRD.size() < 2) {
    mprinterr("Error: %s: Must specify at least 2 COORDS data sets\n", argIn.Command());
    return Cmd::ERR;
  }
  // Only add the topology to the list if parmname specified
  bool addTop = true;
  Topology CombinedTop;
  if (parmname.empty()) {
    parmname = CRD[0]->Top().ParmName() + "_" + CRD[1]->Top().ParmName();
    addTop = false;
  }
  CombinedTop.SetParmName( parmname, FileName() );
  // TODO: Check Parm box info.
  size_t minSize = CRD[0]->Size();
  for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum) {
    if (CRD[setnum]->Size() < minSize)
      minSize = CRD[setnum]->Size();
    CombinedTop.AppendTop( CRD[setnum]->Top() );
  }
  CombinedTop.Brief("Combined parm:");
  if (addTop) {
    if (State.AddTopology( CombinedTop, parmname )) return Cmd::ERR;
  }
  // Combine coordinates
  if (crdname.empty())
    crdname = CRD[0]->Meta().Legend() + "_" + CRD[1]->Meta().Legend();
  mprintf("\tCombining %zu frames from each set into %s\n", minSize, crdname.c_str());
  DataSet_Coords* CombinedCrd = (DataSet_Coords*)State.DSL()->AddSet(DataSet::COORDS, crdname, "CRD");
  if (CombinedCrd == 0) {
    mprinterr("Error: Could not create COORDS data set.\n");
    return Cmd::ERR;
  }
  // FIXME: Only copying coords for now
  CombinedCrd->CoordsSetup( CombinedTop, CoordinateInfo() );
  Frame CombinedFrame( CombinedTop.Natom() * 3 );
  std::vector<Frame> InputFrames;
  for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum)
    InputFrames.push_back( CRD[setnum]->AllocateFrame() );
  for (size_t nf = 0; nf != minSize; ++nf) {
    CombinedFrame.ClearAtoms();
    for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum)
    {
      CRD[setnum]->GetFrame( nf, InputFrames[setnum] );
      for (int atnum = 0; atnum < CRD[setnum]->Top().Natom(); atnum++)
        CombinedFrame.AddXYZ( InputFrames[setnum].XYZ(atnum) );
    }
    CombinedCrd->AddFrame( CombinedFrame );
  }
/* FIXME: This code is fast but only works for DataSet_Coords_CRD
  Frame::CRDtype CombinedFrame( CombinedTop->Natom() * 3 );
  for (size_t nf = 0; nf != minSize; ++nf) {
    size_t offset = 0;
    for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum)
    {
      size_t crd_offset = (size_t)CRD[setnum]->Top().Natom() * 3;
      std::copy( CRD[setnum]->CRD(nf).begin(), CRD[setnum]->CRD(nf).begin() + crd_offset,
                 CombinedFrame.begin() + offset );
      offset += crd_offset;
    }
    CombinedCrd->AddCRD( CombinedFrame );
  }
*/
  return Cmd::OK;
}
