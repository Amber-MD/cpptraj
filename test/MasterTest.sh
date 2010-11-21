# This should be sourced from Test run scripts

# If arg is Key=Value, separate into Key and Value
ParseArg() {
  KEY=`echo "$1" | awk 'BEGIN{FS = "=";}{print $1;}'`
  VALUE=`echo "$1" | awk 'BEGIN{FS = "=";}{print $2;}'`
  if [[ $VALUE = $KEY ]] ; then
    VALUE=""
  fi
}

# DoTest(): Compare File1 to File2, print an error if they differ
DoTest() {
  if [[ $NOTEST -eq 0 ]] ; then
    ((NUMTEST++))
    if [[ ! -e $2 ]] ; then
      echo "  $2 not found." >> $TEST_RESULTS
      ((ERR++))
    elif [[ `diff $1 $2 | wc -l` -gt 0 ]] ; then
      #echo "  $1 $2 differ."
      echo "  $1 $2 are different." >> $TEST_RESULTS
      diff $1 $2 >> $TEST_RESULTS 2>&1
      #echo "------------------------------------------------------------" >> $TEST_RESULTS
      ((ERR++))
    else
      #echo "  $2 OK."
      echo "  $2 OK." >> $TEST_RESULTS
    fi
  fi
}

# CheckTest(): Report if the error counter is greater than 0 and exit if so
CheckTest() {
  if [[ $ERR -gt 0 ]] ; then
    echo "  $ERR out of $NUMTEST comparisons failed."
    echo "  $ERR out of $NUMTEST comparisons failed." >> $TEST_RESULTS
    #echo "---------------------------------------------------------"
    #echo "---------------------------------------------------------" >> $TEST_RESULTS
    exit 1
  else
   echo "  $NUMTEST comparisons ok."
   echo "  $NUMTEST comparisons ok." >> $TEST_RESULTS
  fi
}

# RunCpptraj(): Run cpptraj with the given options. Start and stop MPI if requested.
RunCpptraj() {
  echo ""
  # If only cleaning requested no run needed, exit now
  if [[ $CLEAN -eq 1 ]] ; then
    exit 0
  fi
  echo "  CPPTRAJ: $1"
  echo "  CPPTRAJ: $1" >> $TEST_RESULTS
  $STARTMPI
  if [[ ! -z $DEBUG ]] ; then
    echo "$TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $TOP $INPUT $DEBUG >> $OUTPUT 2>>$ERROR"
  fi
  $TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $TOP $INPUT $DEBUG >> $OUTPUT 2>>$ERROR
  $STOPMPI
}

# EndTest(): Called at the end of every test script if no errors found.
EndTest() {
  echo "All $NUMTEST comparisons passed." 
  echo "All $NUMTEST comparisons passed." >> $TEST_RESULTS 
  echo ""
  if [[ ! -z $VALGRIND ]] ; then
    echo "Valgrind summary:"
    grep ERROR $ERROR
    grep heap $ERROR
    grep LEAK $ERROR
    echo ""
  fi
  #echo "---------------------------------------------------------"
  #echo "---------------------------------------------------------" >> $TEST_RESULTS
}

# CleanFiles(): For every arg passed to the function, check for the file and rm it
CleanFiles() {
  while [[ ! -z $1 ]] ; do
    #for RMFILE in `find . -name "$1"` ; do
    if [[ -e $1 ]] ; then
      #echo "  Cleaning $1"
      rm $1
    fi
    #done
    shift
  done
  # If only cleaning requested no run needed, exit now
  if [[ $CLEAN -eq 1 ]] ; then
    exit 0
  fi
}

# Library Checks - Tests that depend on certain libraries like Zlib can run
# these to make sure cpptraj was compiled with that library - exit gracefully
# if not.
CheckZlib() {
  if [[ -z $ZLIB ]] ; then
    echo "This test requires zlib. Cpptraj was compiled without zlib support."
    echo "Skipping test."
    #echo "---------------------------------------------------------"
    exit 0
  fi
}

CheckBzlib() {
  if [[ -z $BZLIB ]] ; then
    echo "This test requires bzlib. Cpptraj was compiled without bzlib support."
    echo "Skipping test."
    #echo "---------------------------------------------------------"
    exit 0
  fi
}

CheckNetcdf() {
  if [[ -z $NETCDFLIB ]] ; then
    echo "This test requires Netcdf. Cpptraj was compiled without Netcdf support."
    echo "Skipping test."
    #echo "---------------------------------------------------------"
    exit 0
  fi
}
#==============================================================================
# CPPTRAJHOME should be defined
if [[ -z $CPPTRAJHOME ]] ; then
  echo "CPPTRAJHOME not defined."
  # Try AMBERHOME
  if [[ -z $AMBERHOME ]] ; then
    echo "Tests require CPPTRAJHOME or AMBERHOME to be defined."
    echo ""
    exit 0
  fi
  echo "Checking for cpptraj directory in AMBERHOME ($AMBERHOME)"
  CPPTRAJHOME=$AMBERHOME/AmberTools/src/cpptraj
fi
if [[ ! -e $CPPTRAJHOME ]] ; then
  echo "Cpptraj directory $CPPTRAJHOME does not exist."
  echo ""
  exit 1
fi
  
# Check for config.h in CPPTRAJHOME
CONFIGH=$CPPTRAJHOME/config.h
if [[ ! -e $CONFIGH ]] ; then
  echo "$CONFIGH not found."
  exit 1
fi
KEY=""
VALUE=""
# Check for libraries in config.h
ParseArg `grep BZLIB $CONFIGH`
BZLIB=$VALUE
ParseArg `grep ZLIB $CONFIGH`
ZLIB=$VALUE
ParseArg `grep NETCDFLIB $CONFIGH` 
NETCDFLIB=$VALUE
# Check for binary location in config.h
ParseArg `grep CPPTRAJBIN $CONFIGH`
if [[ -z $VALUE ]] ; then
  echo "CPPTRAJBIN not found in config.h"
  echo "Trying ../../bin/cpptraj"
  CPPTRAJ=../../bin/cpptraj
else
  CPPTRAJ=$VALUE/cpptraj
fi

TEST_RESULTS=Test_Results.dat
CleanFiles $TEST_RESULTS

# Option defaults
TIME=""
VALGRIND=""
DO_PARALLEL=""
STARTMPI=""
STOPMPI=""
NP=1
MPI=0
TOP=
INPUT=
NOTEST=0
OUTPUT="test.out"
CleanFiles $OUTPUT
#if [[ -e $OUTPUT ]] ; then
#  rm $OUTPUT
#fi
#OUTPUT="/dev/stdout"
ERROR="/dev/stderr"
DEBUG=""
CLEAN=0
while [[ ! -z $1 ]] ; do
  case "$1" in
    "stdout" ) OUTPUT="/dev/stdout" ;;
    "vg"  ) 
      echo "  Using Valgrind."
      VALGRIND="valgrind --tool=memcheck --leak-check=yes --show-reachable=yes" 
      ERROR="valgrind.out"
      CleanFiles $ERROR
      #if [[ -e $ERROR ]] ; then
      #  rm $ERROR
      #fi
    ;;
    "mpi" ) MPI=1 ;;
    "time") TIME="time" ;;
    "np"  ) 
      shift
      NP=$1
    ;;
    "-i" )
      shift
      INPUT=$1
      echo "Using input file: $INPUT"
      NOTEST=1
    ;;
    "-p" )
      shift
      TOP="$TOP -p $1"
      echo "Using top file: $1"
      NOTEST=1
    ;;
    "notest" )
      echo "End of run tests will be skipped."
      NOTEST=1
    ;;
    "-d"    ) DEBUG="-debug 2" ;;
    "clean" ) CLEAN=1 ;;
    * ) echo "Unknown opt: $1" ;;
  esac
  shift
done

# If not cleaning, check for binary
if [[ $CLEAN -eq 0 ]] ; then
  # If not defined, attempt default location
  if [[ -z $CPPTRAJ ]] ; then
    CPPTRAJ=../../bin/cpptraj
  fi
  if [[ ! -e $CPPTRAJ ]] ; then
    echo "CPPTRAJ not found ($CPPTRAJ)."
    exit 1
  fi
  if [[ ! -z $DEBUG ]] ; then
    ls -l -t $CPPTRAJ
  fi
  # Start test results file
  echo "---------------------------------------------------------"
  echo "TEST: `pwd`" 
  echo "---------------------------------------------------------" > $TEST_RESULTS
  echo "TEST: `pwd`" >> $TEST_RESULTS
fi

# Set up MPI environment if specified or if >1 processor requested.
if [[ $MPI -eq 1 || $NP -gt 1 ]] ; then
  echo "  Using MPI with $NP processors."
  STARTMPI="mpdboot -n 1"
  STOPMPI="mpdallexit"
  DO_PARALLEL="mpiexec -n $NP"
fi

NUMTEST=0
ERR=0

