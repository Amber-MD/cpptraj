# Source this from all unit test scripts

# The variable UNITSOURCES will list all cpptraj source files needed by the unit test.

# ----- Environment Variables --------------------------------------------------
#   CPPTRAJ_TEST_SETUP   : 'yes' if setup is complete
#   CPPTRAJ_RM           : Command used to remove files
#   CPPTRAJ_TEST_RESULTS : File to record individual test results to.
#   CPPTRAJ_ERROR        : File to direct cpptraj STDERR to.
#   CPPTRAJ_TEST_CLEAN   : If 1, only cleaning tests; do not run them.

# ----- Variables local to single test ---------------------
PROGERROR=0              # Total number of program errors this test

# ==============================================================================
# TestHeader() <outfile>
#   Write test header (working directory) to specified file. If no file
#   specified will be written to STDOUT.
TestHeader() {
  if [ ! -z "$1" ] ; then
    echo "********************************************************************************" > $1
    echo "TEST: $TEST_WORKDIR" >> $1
  else
    echo "********************************************************************************"
    echo "TEST: $TEST_WORKDIR"
  fi
}

# --------------------------------------
# ErrMsg() <message>
# Write out error message to stderr prefaced with 'Error:'.
ErrMsg() {
  >&2 echo "Error: $*"
}

# --------------------------------------
# OUT() <message>
#   Write message to CPPTRAJ_TEST_RESULTS
OUT() {
  echo "$1" >> $CPPTRAJ_TEST_RESULTS
}

# ------------------------------------------------------------------------------
# CleanFiles() <file1> ... <fileN>
#   For every arg passed to the function, check for the file or directory and
#   remove it.
CleanFiles() {
  while [ ! -z "$1" ] ; do
    if [ -d "$1" ] ; then
      rmdir $1
    elif [ -f "$1" ] ; then
      $CPPTRAJ_RM $1
    fi
    shift
  done
  # If only cleaning requested no run needed, exit now
  if [ $CPPTRAJ_TEST_CLEAN -eq 1 ] ; then
    exit 0
  fi
}


# ----- Create Makefile ----------------
CreateMakefile() {
  cat > Makefile <<EOF
include $CPPTRAJ_CONFIGH

main.o : main.cpp
	\$(CXX) \$(CXXFLAGS) -I$CPPTRAJ_SRCDIR -c -o main.o main.cpp
EOF
  UNITOBJECTS='main.o'
  for CSRCFILE in $UNITSOURCES ; do
    # Create object name
    ONAME=${CSRCFILE/.cpp/.o}
    UNITOBJECTS="$UNITOBJECTS $ONAME"
    # Add the rule
  cat >> Makefile <<EOF

$ONAME: $CPPTRAJ_SRCDIR/$CSRCFILE
	\$(CXX) \$(DIRECTIVES) \$(CXXFLAGS) -I$CPPTRAJ_SRCDIR -c -o $ONAME $CPPTRAJ_SRCDIR/$CSRCFILE
EOF
  done
  cat >> Makefile <<EOF
all: a.out

a.out: $UNITOBJECTS
	\$(CXX) -o a.out $UNITOBJECTS

EOF
}

# --------------------------------------
# ProgramError() <message>
ProgramError() {
  ErrMsg " $1"
  OUT "Error: $1"
  ((PROGERROR++))
} 

# ----- Run make, check result ---------
RunMake() {
  ((PROGCOUNT++))
  echo ""
  echo "  CPPTRAJ: $1"
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    OUT "  CPPTRAJ: $1"
  fi
  make all
  if [ $? -ne 0 -o ! -f 'a.out' ] ; then
    ProgramError "$1"
  else
    if [ -z "$CPPTRAJ_ERROR" ] ; then
      ./a.out
    else
      $VALGRIND ./a.out 2>> $CPPTRAJ_ERROR
    fi
    if [ $? -ne 0 ] ; then
      ProgramError "$1"
    fi
  fi
}

# ------------------------------------------------------------------------------
# EndTest()
#   Print a summary of the current tests results. Should be called at the end of
#   every test script.
EndTest() {
  #echo "DEBUG: EndTest"
  echo ""
  echo "  $PROGERROR out of $PROGCOUNT executions exited with an error."
  echo "  $PROGERROR out of $PROGCOUNT executions exited with an error." >> $CPPTRAJ_TEST_RESULTS
}

# ------------------------------------------------------------------------------
CmdLineOpts() {
  CPPTRAJ_TEST_CLEAN=0 # Will be exported
  while [ ! -z "$1" ] ; do
    case "$1" in
      "clean"     ) CPPTRAJ_TEST_CLEAN=1 ;;
    esac
    shift
  done
  export CPPTRAJ_TEST_CLEAN
}

# ==============================================================================
# M A I N
# ==============================================================================
if [ -z "$CPPTRAJ_TEST_SETUP" ] ; then
  if [ -z "$CPPTRAJHOME" -o ! -d "$CPPTRAJHOME" ] ; then
    echo "Cannot find CPPTRAJ home directory, set CPPTRAJHOME."
    exit 1
  fi

  # Determine cpptraj source directory, CPPTRAJ_SRCDIR
  if [ -z "$CPPTRAJ_SRCDIR" ] ; then

    CPPTRAJ_SRCDIR=$CPPTRAJHOME/src
    echo "  CPPTRAJ source directory: $CPPTRAJ_SRCDIR"
    if [ ! -d "$CPPTRAJ_SRCDIR" ] ; then
      echo "CPPTRAJ source directory not found."
      exit 1
    fi

    if [ ! -f "$CPPTRAJ_SRCDIR/main.cpp" ] ; then
      echo "CPPTRAJ source directory $CPPTRAJ_SRCDIR appears empty."
      exit 1
    fi
    export CPPTRAJ_SRCDIR
  fi

  # Determine the config.h file, CPPTRAJ_CONFIGH
  if [ -z "$CPPTRAJ_CONFIGH" ] ; then
    CPPTRAJ_CONFIGH=$CPPTRAJHOME/config.h
    echo "  CPPTRAJ config.h file: $CPPTRAJ_CONFIGH"
    if [ ! -f "$CPPTRAJ_CONFIGH" ] ; then
      echo "CPPTRAJ config.h file not found."
      exit 1
    fi
    export CPPTRAJ_CONFIGH
  fi

  # Ensure required binaries are set up
  if [ -z "$CPPTRAJ_RM" ] ; then
    # TODO is this being too paranoid?
    if [ ! -f '/bin/rm' ] ; then
      ErrMsg "Required binary '/bin/rm' not found."
      exit 1
    fi
    export CPPTRAJ_RM='/bin/rm -f'
  fi

  # Set some defaults
  export CPPTRAJ_TEST_RESULTS='Test_Results.dat'
  export CPPTRAJ_ERROR=''

  # Process command line options
  CmdLineOpts $*

  export CPPTRAJ_TEST_SETUP='yes'
fi # END initial setup

# Clean any existing test results files.
if [ -f "$CPPTRAJ_TEST_RESULTS" ] ; then
  $CPPTRAJ_RM $CPPTRAJ_TEST_RESULTS
fi

# Determine mode of execution: individual test or multiple tests.
if [ "$CPPTRAJ_TEST_MODE" = 'master' ] ; then
  # Executing multiple tests
  echo Master Make
else
  # Executing single test
  echo non-master make
  if [ $CPPTRAJ_TEST_CLEAN -eq 0 ] ; then
    TEST_WORKDIR=`pwd`
    TestHeader
    TestHeader "$CPPTRAJ_TEST_RESULTS"
  fi
fi
