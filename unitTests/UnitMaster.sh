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
# ----- Local setup variables ------------------------------
SUMMARY=0                # If 1 print summary of CPPTRAJ_TEST_RESULTS only
TARGET=""                # Make target if multiple tests being run

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
#   Unit tests either pass or fail.
EndTest() {
  #echo "DEBUG: EndTest"
  # Sanity check
  if [ $PROGCOUNT -gt 1 ] ; then
    ErrMsg "Too many executions; got $PROGCOUNT, expected 1."
    exit 1
  fi
  echo ""
  if [ $PROGERROR -eq 0 ] ; then
    echo "  Unit test passed."
    echo "  Unit test passed." >> $CPPTRAJ_TEST_RESULTS
  else
    echo "  Unit test failed."
    echo "  Unit test failed." >> $CPPTRAJ_TEST_RESULTS
  fi
  #echo "  $PROGERROR out of $PROGCOUNT executions exited with an error."
  #echo "  $PROGERROR out of $PROGCOUNT executions exited with an error." >> $CPPTRAJ_TEST_RESULTS
  echo ""
  exit $PROGERROR
}

# ------------------------------------------------------------------------------
# Summary()
#  Print a summary of results in all CPPTRAJ_TEST_RESULTS files and exit.
#  Optionally print an error summary and/or valgrind summary.
Summary() {
  ERR_STATUS=0
  if [ ! -z "$CPPTRAJ_TEST_RESULTS" ] ; then
    RESULTSFILES=`ls */$CPPTRAJ_TEST_RESULTS 2> tmp.cpptrajtest.devnull`
    rm tmp.cpptrajtest.devnull
    if [ ! -z "$RESULTSFILES" ] ; then
      cat $RESULTSFILES > $CPPTRAJ_TEST_RESULTS
      echo "===================== TEST SUMMARY ======================"
      awk 'BEGIN{
        program_exe  = 0; # Number of unit test executions
        program_pass = 0; # Number of unit test pass
        program_err  = 0; # Number of unit test failures
      }{
        if ($1 == "CPPTRAJ:")
          program_exe++;
        else if ($1 == "Unit" && $2 == "test") {
          if ($3 == "passed.")
            program_pass++;
          else if ($3 == "failed.")
            program_err++;
        }
      }END{
        if (program_exe > 0)
          printf("  %i out of %i unit tests passed (%i failed).\n",
                 program_pass, program_exe, program_err);
        exit program_err;
      }' $CPPTRAJ_TEST_RESULTS
      echo "========================================================="
      ERR_STATUS=$?
    fi
  fi
  exit $ERR_STATUS
}

# ------------------------------------------------------------------------------
CmdLineOpts() {
  CPPTRAJ_TEST_CLEAN=0 # Will be exported
  while [ ! -z "$1" ] ; do
    case "$1" in
      "clean"     ) CPPTRAJ_TEST_CLEAN=1 ;;
      "summary"   ) SUMMARY=1 ;;
      "--target"  ) shift ; TARGET=$1 ;;
      * ) ErrMsg "Unknown option: $1" ; exit 1 ;;
    esac
    shift
  done
  export CPPTRAJ_TEST_CLEAN
}

# ==============================================================================
# M A I N
# ==============================================================================
if [ -z "$CPPTRAJ_TEST_SETUP" ] ; then
  #if [ -z "$CPPTRAJHOME" -o ! -d "$CPPTRAJHOME" ] ; then
  #  echo "Cannot find CPPTRAJ home directory, set CPPTRAJHOME."
  #  exit 1
  #fi
  # Assume we are in the cpptraj source directory, <dir>/unitTests
  # or <dir>/unitTests/<testname>
  CPPTRAJDIR=''
  CURRENTDIR=`pwd`
  #echo "DEBUG : current dir $CURRENTDIR"
  if [ "`basename $CURRENTDIR`" = 'unitTests' ] ; then
    # This is <dir>/unitTests
    CPPTRAJDIR=`dirname $CURRENTDIR`
  else
    ONE_UP=`dirname $CURRENTDIR`
    if [ "`basename $ONE_UP`" = 'unitTests' ] ; then
      # This is <dir>/unitTests/<dir>
      CPPTRAJDIR=`dirname $ONE_UP`
    else
      ErrMsg "Could not determine cpptraj source directory location."
      exit 1
    fi
  fi
  if [ -z "$CPPTRAJDIR" -o ! -d "$CPPTRAJDIR" ] ; then
    ErrMsg "Cpptraj source directory $CPPTRAJDIR empty or is not a directory."
    exit 1
  fi

  # Determine cpptraj source file directory, CPPTRAJ_SRCDIR
  if [ -z "$CPPTRAJ_SRCDIR" ] ; then
    CPPTRAJ_SRCDIR=$CPPTRAJDIR/src
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
    CPPTRAJ_CONFIGH=$CPPTRAJDIR/config.h
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
  # Probably executed via make.
  # If summary requested just do that and exit.
   if [ $SUMMARY -ne 0 ] ; then
    Summary
  fi
  # Need a Makefile.
  if [ ! -f 'Makefile' ] ; then
    ErrMsg "test Makefile not found."
    exit 1
  fi
  #Required "make"
  if [ -z "$TARGET" ] ; then
    # If no target specified, probably not executed via make.
    ErrMsg "No test directories specified."
    exit 1
  fi
  make $TARGET
  if [ $? -ne 0 ] ; then
    exit 1
  fi
else
  # Executing single test
  if [ $CPPTRAJ_TEST_CLEAN -eq 0 ] ; then
    TEST_WORKDIR=`pwd`
    TestHeader
    TestHeader "$CPPTRAJ_TEST_RESULTS"
  fi
fi
