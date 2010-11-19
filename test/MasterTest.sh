# This should be sourced from Test run scripts

# DoTest(): Compare File1 to File2, print an error if they differ
DoTest() {
  if [[ $NOTEST -eq 0 ]] ; then
    ((NUMTEST++))
    if [[ ! -e $2 ]] ; then
      echo "  $2 not found."
      ((ERR++))
    elif [[ `diff $1 $2 | wc -l` -gt 0 ]] ; then
      #echo "  $1 $2 differ."
      echo "  $1 $2 differ." >> $TEST_RESULTS
      diff $1 $2 >> $TEST_RESULTS 2>&1
      echo "------------------------------------------------------------" >> $TEST_RESULTS
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
    echo ""
    exit 1
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
  echo "---------------------------------------------------------"
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

#==============================================================================
TEST_RESULTS=Test_Results.dat
CleanFiles $TEST_RESULTS

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


