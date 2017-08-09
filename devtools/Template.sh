#!/bin/bash
# This script can be used to generate a template for common cpptraj components.
# Will create both the .h and .cpp files.
# Daniel R. Roe 2016

Help() {
  echo "Usage: $0 <name> [<type>]"
  echo "  <type>: Action Analysis Exec Traj DataIO"
  echo ""
}

echo "CPPTRAJ file template generator."

NAME=$1
if [ -z "$NAME" ] ; then
  echo "Enter name."
  Help
  exit 1
fi

TYPE=$2
if [ -z "$TYPE" ] ; then
  echo "Enter type."
  Help
  exit 1 
fi
if [ "$TYPE" != 'Action' -a "$TYPE" != 'Analysis' -a "$TYPE" != 'Exec' -a "$TYPE" != 'Traj' -a "$TYPE" != 'DataIO' ] ; then
  echo "Type $TYPE not recognized."
  Help
  exit 1
fi

CLASS=$TYPE"_"$NAME
H_FILE=$CLASS".h"
C_FILE=$CLASS".cpp"
if [ "$TYPE" = 'Traj' ] ; then
  TYPEH='TrajectoryIO'
else
  TYPEH="$TYPE"
fi
echo "Creating class $CLASS in $H_FILE and $C_FILE"

# Check files
if [ -e "$H_FILE" -o -e "$C_FILE" ] ; then
  echo Already there.
  exit 1
fi

# Header protect
cat > $H_FILE <<EOF
#ifndef INC_${TYPE^^}_${NAME^^}_H
#define INC_${TYPE^^}_${NAME^^}_H
#include "$TYPEH.h"
/// <Enter description of $CLASS here>
class $CLASS : public $TYPEH {
EOF

# Class-specific header/implementation
# ----- Exec -------------------------------------
if [ "$TYPE" = 'Exec' ] ; then
  cat >> $H_FILE <<EOF
  public:
    $CLASS() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new $CLASS(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
EOF
  cat > $C_FILE <<EOF
#include "$H_FILE"
#include "CpptrajStdio.h"

// $CLASS::Help()
void $CLASS::Help() const
{

}

// $CLASS::Execute()
Exec::RetType $CLASS::Execute(CpptrajState& State, ArgList& argIn)
{

}
EOF
# ----- Action -----------------------------------
elif [ "$TYPE" = 'Action' ] ; then
  cat >> $H_FILE <<EOF
  public:
    $CLASS() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new $CLASS(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

};
#endif
EOF
  cat > $C_FILE <<EOF
#include "$H_FILE"
#include "CpptrajStdio.h"

// $CLASS::Help()
void $CLASS::Help() const {

}

// $CLASS::Init()
Action::RetType $CLASS::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{

}

// $CLASS::Setup()
Action::RetType $CLASS::Setup(ActionSetup& setup)
{

}

// $CLASS::DoAction()
Action::RetType $CLASS::DoAction(int frameNum, ActionFrame& frm)
{

}
EOF
# ----- Analysis ---------------------------------
elif [ "$TYPE" = 'Analysis' ] ; then
  cat >> $H_FILE <<EOF
  public:
    $CLASS() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new $CLASS(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:

};
#endif
EOF
  cat > $C_FILE <<EOF
#include "$H_FILE"
#include "CpptrajStdio.h"

// $CLASS::Help()
void $CLASS::Help() const {

}

// $CLASS::Setup()
Analysis::RetType $CLASS::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{

}

// $CLASS::Analyze()
Analysis::RetType $CLASS::Analyze() {

}
EOF
# ----- Traj -------------------------------------
elif [ "$TYPE" = 'Traj' ] ; then
  cat >> $H_FILE <<EOF
  public:
    $CLASS();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new $CLASS(); }
    static void WriteHelp();
    static void ReadHelp();
  private:
    // ----- Inherited functions -----------------
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int readVelocity(int, Frame&);
    int readForce(int, Frame&);
    int processWriteArgs(ArgList&);
    int processReadArgs(ArgList&);
    // -------------------------------------------
#   ifdef MPI
    // ----- Parallel functions ------------------
    int parallelOpenTrajin(Parallel::Comm const&);
    int parallelOpenTrajout(Parallel::Comm const&);
    int parallelSetupTrajout(FileName const&, Topology*, CoordinateInfo const&,
                             int, bool, Parallel::Comm const&);
    int parallelReadFrame(int, Frame&);
    int parallelWriteFrame(int, Frame const&);
    void parallelCloseTraj();
    // -------------------------------------------
#   endif
};
#endif
EOF
  cat > $C_FILE <<EOF
#include "$H_FILE"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
$CLASS::$CLASS() {}

/** Identify trajectory format. File should be setup for READ */
bool $CLASS::ID_TrajFormat(CpptrajFile& fileIn) {

  return false;
}

/** Print trajectory info to stdout. */
void $CLASS::Info() {
  mprintf("is a <type>");
}

/** Close file. */
void $CLASS::closeTraj() {

}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
int $CLASS::openTrajin() {

  return 0;
}

/** Read help */
void $CLASS::ReadHelp() {

}

/** Process read arguments. */
int $CLASS::processReadArgs(ArgList& argIn) {

  return 0;
}

/** Set up trajectory for reading.
  * \return Number of frames in trajectory.
  */
int $CLASS::setupTrajin(FileName const& fname, Topology* trajParm)
{

  return TRAJIN_ERR;
}

/** Read specified trajectory frame. */
int $CLASS::readFrame(int set, Frame& frameIn) {

  return 0;
}

/** Read velocities from specified frame. */
int $CLASS::readVelocity(int set, Frame& frameIn) {

  return 0;
}

/** Read forces from specified frame. */
int $CLASS::readForce(int set, Frame& frameIn) {

  return 0;
}

// -----------------------------------------------------------------------------
/** Write help. */
void $CLASS::WriteHelp() {

}

/** Process write arguments. */
int $CLASS::processWriteArgs(ArgList& argIn) {

  return 0;
}

/** Set up trajectory for write. */
int $CLASS::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{

  return 1;
}

/** Write specified trajectory frame. */
int $CLASS::writeFrame(int set, Frame const& frameOut) {

  return 0;
}

// =============================================================================
#ifdef MPI
/** Open trajectory for reading in parallel. */
int $CLASS::parallelOpenTrajin(Parallel::Comm const& commIn) {
  return 1;
}

/** Open trajectory for writing in parallel. */
int $CLASS::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return 1;
}

/** Set up trajectory for write in parallel. */
int $CLASS::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{

  return 1;
}

/** Read frame in parallel. */
int $CLASS::parallelReadFrame(int set, Frame& frameIn) {

  return 1;
}

/** Write frame in parallel. */
int $CLASS::parallelWriteFrame(int set, Frame const& frameOut) {

  return 1;
}

/** Close trajectory in parallel. */
void $CLASS::parallelCloseTraj() {

}
#endif
EOF
# ----- DataIO -----------------------------------
elif [ "$TYPE" = 'DataIO' ] ; then
  cat >> $H_FILE <<EOF
  public:
    $CLASS();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new $CLASS(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
};
#endif
EOF
  cat > $C_FILE <<EOF
#include "$H_FILE"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
$CLASS::$CLASS()
{

}

bool $CLASS::ID_DataFormat(CpptrajFile& infile)
{

  return false;
}

// $CLASS::ReadHelp()
void $CLASS::ReadHelp()
{

}

// $CLASS::WriteHelp()
void $CLASS::WriteHelp()
{

}

// $CLASS::processReadArgs()
int $CLASS::processReadArgs(ArgList& argIn)
{

  return 0;
}

// $CLASS::ReadData()
int $CLASS::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{

  return 1;
}

// $CLASS::processWriteArgs()
int $CLASS::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// $CLASS::WriteData()
int $CLASS::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
EOF

# ------------------------------------------------
else
  echo "Unrecognized type: $TYPE."
  exit 1
fi
