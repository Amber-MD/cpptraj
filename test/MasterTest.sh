# This should be sourced at the top of CPPTRAJ test run scripts.
#
# The purpose of this script is to provide all common functionality required
# by all tests, including environment setup and collecting test results.
# Requires awk, grep, sed, diff (if not AmberTools), rm, make (multiple tests)
#
# -----===== Environment variables =====-----
# Binary locations
#   CPPTRAJ              : Set to cpptraj binary being tested.
#   AMBPDB               : Set to ambpdb binary being tested.
#   VALGRIND             : Set to 'valgrind' command if memory check requested.
#   CPPTRAJ_DIFF         : Command used to check for test differences.
#   CPPTRAJ_DACDIF       : Set if using 'dacdif' in AmberTools to checkout for test differences.
#   NDIFF                : Set to Nelson H. F. Beebe's ndiff.awk for numerical diff calc.
#   CPPTRAJ_NCDUMP       : Set to ncdump command; needed for NcTest()
#   REMOVE               : Command used to remove files
#   TIME                 : Set to the 'time' command if timing requested.
# Test output locations
#   CPPTRAJ_TEST_RESULTS : File to record test results to.
#   CPPTRAJ_TEST_ERROR   : File to record test errors/diffs to.
#   CPPTRAJ_OUTPUT       : File to direct cpptraj STDOUT to.
#   CPPTRAJ_ERROR        : File to direct cpptraj STDERR to.
# Other variables
#   CPPTRAJ_TEST_ROOT    : Test root directory.
#   CPPTRAJ_TEST_MODE    : 'single' or 'multiple'.
#   CPPTRAJ_STANDALONE   : If 0, part of AmberTools. If 1, stand-alone (e.g. from GitHub).
#   CPPTRAJ_TEST_CLEAN   : If 1, only cleaning tests; do not run them.
#   TEST_OS              : Operating system on which tests are being run. If blank assume linux.
#   N_THREADS            : Set to number of MPI threads if parallel.
#   OMP_NUM_THREADS      : Set to max number of OpenMP threads.
#   DO_PARALLEL          : Set to the MPI run command (e.g. 'mpirun -n 11')
#   CPPTRAJ_DEBUG        : Can be set to pass global debug flag to cpptraj.
# Cpptraj binary characteristics 
#   CPPTRAJ_DEFINES      : Set to the output of 'CPPTRAJ --defines'
#   CPPTRAJ_ZLIB         : If set CPPTRAJ has zlib support.
#   CPPTRAJ_BZLIB
#   CPPTRAJ_NETCDFLIB
#   CPPTRAJ_MPILIB
#   CPPTRAJ_NOMATHLIB
#   CPPTRAJ_OPENMP
#   CPPTRAJ_PNETCDFLIB
#   CPPTRAJ_SANDERLIB

# FIXME Variables to check
SUMMARY=0           # If 1, only summary of results needs to be performed.
SHOWERRORS=0        # If 1, print test errors to STDOUT after summary.
CPPTRAJ_PROFILE=0           # If 1, end of test profiling with gprof performed #FIXME
USE_DACDIF=1         # If 0 do not use dacdif even if in AmberTools
SFX=""              # CPPTRAJ binary suffix
#NPROC=""            # nproc binary for counting threads in parallel tests.

# Variables local to single test.
NUMTEST=0           # Total number of times DoTest has been called this test.
ERRCOUNT=0          # Total number of errors detected by DoTest this test.
WARNCOUNT=0         # Total number of warnings detected by DoTest this test.
PROGERROR=0         # Total number of program errors this test

# ==============================================================================
# Output() <Message>
#   Send <Message> to CPPTRAJ_TEST_RESULTS and CPPTRAJ_TEST_ERROR
Output() {
  echo "$1" >> $CPPTRAJ_TEST_RESULTS
  echo "$1" >> $CPPTRAJ_TEST_ERROR
}

# ------------------------------------------------------------------------------
# DoTest() <File1> <File2> [-r <relative err>] [-a <absolute err>] [<arg1>] ... [<argN>]
#   Compare File1 (the 'save' file) to File2 (test output), print an error if
#   they differ. If '-r' or '-a' are specified the test will pass as long as the
#   maximum relative or absolulte errors are below the given value (using Nelson
#   H. F. Beebe's ndiff.awk script. The remaining args can be used to pass
#   options to CPPTRAJ_DIFF.
DoTest() {
  echo "DEBUG: DoTest $1 $2"
  if [ ! -z "$CPPTRAJ_DACDIF" ] ; then
    # AmberTools - use dacdif. Use any '-r <X>' or '-a <X>' args found.
    # Ignore the rest.
    DIFFARGS="$1 $2"
    shift # Save file
    shift # Test file
    # Process remaining args
    while [ ! -z "$1" ] ; do
      case "$1" in
        "-r" ) shift ; DIFFARGS=" -r $1 "$DIFFARGS ;;
        "-a" ) shift ; DIFFARGS=" -a $1 "$DIFFARGS ;;
      esac
      shift
    done
    $CPPTRAJ_DACDIF $DIFFARGS
  else
    # Standalone - will use diff, or ndiff where '-r' or '-a' specified.
    ((NUMTEST++))
    DIFFARGS='--strip-trailing-cr'
    NDIFFARGS=""
    # First two arguments are files to compare.
    F1=$1 ; shift
    F2=$1 ; shift
    # Process remaining arguments.
    USE_NDIFF=0
    while [ ! -z "$1" ] ; do
      case "$1" in
        "-r"           ) USE_NDIFF=1; shift; NDIFFARGS="$NDIFFARGS -v RELERR=$1" ;;
        "-a"           ) USE_NDIFF=1; shift; NDIFFARGS="$NDIFFARGS -v ABSERR=$1" ;;
        *              ) DIFFARGS=$DIFFARGS" $1" ;;
      esac
      shift
    done
    if [ ! -f "$F1" ] ; then
      echo "  $F1 not found." >> $CPPTRAJ_TEST_RESULTS
      echo "  $F1 not found." >> $CPPTRAJ_TEST_ERROR
      ((ERRCOUNT++))
    elif [ ! -f "$F2" ] ; then
      echo "  $F2 not found." >> $CPPTRAJ_TEST_RESULTS
      echo "  $F2 not found." >> $CPPTRAJ_TEST_ERROR
      ((ERRCOUNT++))
    else
      if [ $USE_NDIFF -eq 0 ] ; then
        $CPPTRAJ_DIFF $DIFFARGS $DIFFOPTS $F1 $F2 > temp.diff 2>&1
      else
        $NDIFF $NDIFFARGS $F1 $F2 > temp.diff 2>&1
      fi
      if [ -s 'temp.diff' ] ; then
        echo "  $F1 $F2 are different." >> $CPPTRAJ_TEST_RESULTS
        echo "  $F1 $F2 are different." >> $CPPTRAJ_TEST_ERROR
        cat temp.diff >> $CPPTRAJ_TEST_ERROR
        ((ERRCOUNT++))
      else
        echo "  $F2 OK." >> $CPPTRAJ_TEST_RESULTS
      fi
      $REMOVE temp.diff
    fi
  fi
}

# ------------------------------------------------------------------------------
# NcTest() <1> <2>
#   Compare NetCDF files <1> and <2>. Use CPPTRAJ_NCDUMP to convert to ASCII
#   first, removing ==> line and :programVersion attribute.
NcTest() {
  echo "DEBUG: NcTest $1 $2"
  if [ -z "$1" -o -z "$2" ] ; then
    echo "Error: NcTest(): One or both files not specified." > /dev/stderr
    exit 1
  fi
  if [ -z "$CPPTRAJ_NCDUMP" -o ! -e "$CPPTRAJ_NCDUMP" ] ; then
    echo "ncdump missing." > /dev/stderr
    exit 1
  fi
  # Save remaining args for DoTest
  F1=$1
  F2=$2
  shift
  shift
  DIFFARGS="nc0.save nc0"
  CALC_NUM_ERR=0
  while [ ! -z "$1" ] ; do
    if [ "$1" = '-r' -o "$1" = '-a' ] ; then
      CALC_NUM_ERR=1
    fi
    DIFFARGS=$DIFFARGS" $1"
    shift
  done
  # Prepare files.
  if [ ! -e "$F1" ] ; then
    echo "Error: $F1 missing." >> $CPPTRAJ_TEST_ERROR
  elif [ ! -e "$F2" ] ; then
    echo "Error: $F2 missing." >> $CPPTRAJ_TEST_ERROR
  else
    if [ $CALC_NUM_ERR -eq 1 ] ; then
      # FIXME: Must remove commas here because I cannot figure out how to pass
      # the regular expression to ndiff.awk FS without the interpreter giving
      # this error for FS='[ \t,()]':
      # awk: fatal: Invalid regular expression: /'[/
      $CPPTRAJ_NCDUMP -n nctest $F1 | grep -v "==>\|:programVersion" | sed 's/,/ /g' > nc0.save
      $CPPTRAJ_NCDUMP -n nctest $F2 | grep -v "==>\|:programVersion" | sed 's/,/ /g' > nc0
    else
      $CPPTRAJ_NCDUMP -n nctest $F1 | grep -v "==>\|:programVersion" > nc0.save
      $CPPTRAJ_NCDUMP -n nctest $F2 | grep -v "==>\|:programVersion" > nc0
    fi
    DoTest $DIFFARGS 
    $REMOVE nc0.save nc0
  fi
}

# ------------------------------------------------------------------------------
# RunCpptraj() <title>
#   Run cpptraj test with given title.
RunCpptraj() {
  # If only cleaning requested no run needed, exit now
  if [ $CPPTRAJ_TEST_CLEAN -eq 1 ] ; then
    exit 0
  fi
  echo ""
  echo "  CPPTRAJ: $1"
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    echo "  CPPTRAJ: $1" >> $CPPTRAJ_TEST_RESULTS
  fi
#  if [[ ! -z $DEBUG ]] ; then
    echo "$TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $DEBUG $TOP $INPUT >> $CPPTRAJ_OUTPUT 2>>$CPPTRAJ_ERROR"
#  fi
  $TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $DEBUG $TOP $INPUT >> $CPPTRAJ_OUTPUT 2>>$CPPTRAJ_ERROR
  STATUS=$?
  echo "DEBUG: Cpptraj exited with status $STATUS"
  if [ $STATUS -ne 0 ] ; then
    echo "Error: cpptraj exited with status $STATUS" 2> /dev/stderr
    echo "Error: cpptraj exited with status $STATUS" > $CPPTRAJ_TEST_RESULTS
  fi
}

# ------------------------------------------------------------------------------
# EndTest()
#   Called at the end of every test script if no errors found.
EndTest() {
  echo "DEBUG: EndTest"
  # Report only when not using dacdif 
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    if [ $ERRCOUNT -gt 0 ] ; then
      echo "  $ERRCOUNT out of $NUMTEST comparisons failed."
      echo "  $ERRCOUNT out of $NUMTEST comparisons failed." >> $CPPTRAJ_TEST_RESULTS
      echo "  $ERRCOUNT out of $NUMTEST comparisons failed." >> $CPPTRAJ_TEST_ERROR
    elif [ $WARNCOUNT -gt 0 ] ; then
      ((PASSCOUNT = $NUMTEST - $WARNCOUNT))
      echo "  $PASSCOUNT out of $NUMTEST passing comparisons. $WARNCOUNT warnings."
      echo "  $PASSCOUNT out of $NUMTEST passing comparisons. $WARNCOUNT warnings." >> $CPPTRAJ_TEST_RESULTS
    else 
      echo "All $NUMTEST comparisons passed." 
      echo "All $NUMTEST comparisons passed." >> $CPPTRAJ_TEST_RESULTS 
    fi
    echo ""
    if [ ! -z "$VALGRIND" ] ; then
      echo "Valgrind summary:"
      grep ERROR $CPPTRAJ_ERROR
      grep heap $CPPTRAJ_ERROR
      grep LEAK $CPPTRAJ_ERROR
      echo ""
      echo "Valgrind summary:" >> $CPPTRAJ_TEST_RESULTS
      grep ERROR $CPPTRAJ_ERROR >> $CPPTRAJ_TEST_RESULTS
      grep heap $CPPTRAJ_ERROR >> $CPPTRAJ_TEST_RESULTS
      grep LEAK $CPPTRAJ_ERROR >> $CPPTRAJ_TEST_RESULTS
    fi
    if [ $CPPTRAJ_PROFILE -eq 1 ] ; then
      if [ -e 'gmon.out' ] ; then
        gprof $CPPTRAJ > profiledata.txt
      fi
    fi
  fi
}

# ------------------------------------------------------------------------------
# CleanFiles() <file1> ... <fileN>
#   For every arg passed to the function, check for the file and remove it.
CleanFiles() {
  while [ ! -z "$1" ] ; do
    if [ -e "$1" ] ; then
      $REMOVE $1
    fi
    shift
  done
  # If only cleaning requested no run needed, exit now
  if [ $CPPTRAJ_TEST_CLEAN -eq 1 ] ; then
    exit 0
  fi
}

#-------------------------------------------------------------------------------
# NotParallel() <Test title>
NotParallel() {
  if [ ! -z "$DO_PARALLEL" ] ; then
    echo ""
    echo "  CPPTRAJ: $1"
    echo "  This test cannot be run in parallel. Skipping test."
    return 1
 fi
 return 0
}

#-------------------------------------------------------------------------------
# SetNthreads()
# Use NPROC to set N_THREADS if not already set.
SetNthreads() {
  if [ -z "$N_THREADS" ] ; then
    if [ ! -f "$NPROC" ] ; then
      return 1
    fi
    N_THREADS=`$DO_PARALLEL $NPROC`
  fi
  return 0
}

#-------------------------------------------------------------------------------
# RequiresThreads() <# threads> <Test title>
RequiresThreads() {
  if [ ! -z "$DO_PARALLEL" ] ; then
    SetNthreads
    if [ $? -ne 0 ] ; then
      echo "Error: Program to find # threads not found ($NPROC)" > /dev/stderr
      echo "Error: Test requires $1 parallel threads. Attempting to run test anyway." > /dev/stderr
      return 0
    fi
    REMAINDER=`echo "$N_THREADS % $1" | bc`
    if [ -z "$REMAINDER" -o $REMAINDER -ne 0 ] ; then
      echo ""
      if [ ! -z "$2" ] ; then
        echo "  CPPTRAJ: $2"
      fi
      echo "  Warning: Test requires a multiple of $1 parallel threads. Skipping."
      return 1
    fi
  fi
  return 0
}

#-------------------------------------------------------------------------------
# MaxThreads() <# threads> <Test title>
MaxThreads() {
  if [ ! -z "$DO_PARALLEL" ] ; then
    SetNthreads
    if [ $? -ne 0 ] ; then
      echo "Error: Program to find # threads not found ($NPROC)" > /dev/stderr
      echo "Error: Test can only run with $1 or fewer threads. Attempting to run test anyway." > /dev/stderr
      return 0
    fi
    if [ $N_THREADS -gt $1 ] ; then
      echo ""
      if [ ! -z "$2" ] ; then
        echo "  CPPTRAJ: $2"
      fi
      echo "  Warning: Test can only run with $1 or fewer parallel threads. Skipping."
      return 1
    fi
  fi
  return 0
}

#-------------------------------------------------------------------------------
# Help(): Print help
Help() {
  echo "Command line flags:"
  echo "  summary    : Print summary of test results only."
  echo "  showerrors : (summary only) Print all test errors to STDOUT after summary."
  echo "  stdout     : Print CPPTRAJ test output to STDOUT."
  echo "  mpi        : Use MPI version of CPPTRAJ (automatically triggered if DO_PARALLEL set)."
  echo "  openmp     : Use OpenMP version of CPPTRAJ."
  echo "  cuda       : Use CUDA version of CPPTRAJ."
  echo "  vg         : Run test with valgrind memcheck."
  echo "  vgh        : Run test with valgrind helgrind."
  echo "  time       : Time the test."
  echo "  clean      : Clean test output."
  #echo "  -at        : Force AmberTools tests."
  echo "  -nodacdif  : Do not use dacdif for test comparisons."
  echo "  -d         : Run CPPTRAJ with global debug level 4."
  echo "  -debug <#> : Run CPPTRAJ with global debug level #."
  echo "  -cpptraj <file> : Use CPPTRAJ binary <file>."
  echo "  -ambpdb <file>  : Use AMBPDB binary <file>."
  echo "  -profile        : Profile results with 'gprof' (requires special compile)."
  echo "Important environment variables:"
  echo "  DO_PARALLEL: MPI run command."
  echo "  N_THREADS  : Number of MPI threads (only needed if 'DO_PARALLEL nproc' fails)."
  echo ""
}

#-------------------------------------------------------------------------------
# CmdLineOpts()
#   Process test script command line options. Only executed if CPPTRAJ_TEST_MODE
#   is not already set.
CmdLineOpts() {
  CPPTRAJ_TEST_CLEAN=0 # Will be exported
  VGMODE=0 # Valgrind mode: 0 none, 1 memcheck, 2 helgrind
  SFX_OMP=0
  SFX_CUDA=0
  SFX_MPI=0
  GET_TIMING=0
  while [ ! -z "$1" ] ; do
    case "$1" in
      "clean"     ) CPPTRAJ_TEST_CLEAN=1 ; break ;;
#     "summary"   ) SUMMARY=1 ;;
#     "showerrors") SHOWERRORS=1 ;;
#     "stdout"    ) OUTPUT="/dev/stdout" ;;
     "openmp"    ) SFX_OMP=1 ;;
     "cuda"      ) SFX_CUDA=1 ;;
     "mpi"       ) SFX_MPI=1 ;;
#     "vg"        ) VGMODE=1 ;;
#     "vgh"       ) VGMODE=2 ;;
#     "time"      ) GET_TIMING=0 ;;
#      "-at"       ) FORCE_AMBERTOOLS=1 ;;
#     "-d"        ) DEBUG="-debug 4" ;;
#     "-debug"    ) shift ; DEBUG="-debug $1" ;;
     "-nodacdif" ) USE_DACDIF=0 ;;
#     "-cpptraj"  ) shift ; CPPTRAJ=$1 ; echo "Using cpptraj: $CPPTRAJ" ;;
#     "-ambpdb"   ) shift ; AMBPDB=$1  ; echo "Using ambpdb: $AMBPDB" ;;
#     "-profile"  ) PROFILE=1 ; echo "Performing gnu profiling during EndTest." ;;
      "-h" | "--help" ) Help ; exit 0 ;;
      *           ) echo "Error: Unknown opt: $1" > /dev/stderr ; exit 1 ;;
    esac
    shift
  done
  export CPPTRAJ_TEST_CLEAN
  # Set up SFX
  if [ "$SFX_MPI"  -eq 1 ] ; then SFX="$SFX.MPI"  ; fi
  if [ "$SFX_OMP"  -eq 1 ] ; then SFX="$SFX.OMP"  ; fi
  if [ "$SFX_CUDA" -eq 1 ] ; then SFX="$SFX.cuda" ; fi
}

#-------------------------------------------------------------------------------
# SetBinaries()
#   Set paths for all binaries if not already set.
SetBinaries() {
  # Determine location of CPPTRAJ if not specified.
  if [ -z "$CPPTRAJ" ] ; then
     if [ $CPPTRAJ_STANDALONE -eq 0 ] ; then
       DIRPREFIX=$AMBERHOME
     elif [ -z "$CPPTRAJHOME" ] ; then
       DIRPREFIX=../../
     else
       DIRPREFIX=$CPPTRAJHOME
     fi
     CPPTRAJ=$DIRPREFIX/bin/cpptraj$SFX
     if [ ! -f "$CPPTRAJ" ] ; then
       echo "Error: CPPTRAJ $CPPTRAJ not found." > /dev/stderr
       exit 1
     fi
     export CPPTRAJ
  fi
}

#-------------------------------------------------------------------------------
# Required() <binary>
#   Insure that specified binary is in the PATH
Required() {
  if [ -z "`which $1`" ] ; then
    echo "Error: Required binary '$1' not found." > /dev/stderr
    exit 1
  fi
}

# ==============================================================================
# M A I N
# ==============================================================================

if [ -z "$CPPTRAJ_TEST_MODE" ] ; then
  echo "DEBUG: Initial test setup."
  # MasterTest.sh has not been called yet; set up test environment.
  # Determine mode of execution: individual test or multiple tests.
  if [ -f 'RunTest.sh' ] ; then
    # Assume we are only executing a single test.
    CPPTRAJ_TEST_MODE='single'
  else
    # Assume we are executing multiple tests. Need a Makefile.
    if [ ! -f 'Makefile' ] ; then
      echo "Error: test Makefile not found." > /dev/stderr
      exit 1
    fi
    CPPTRAJ_TEST_MODE='multiple'
  fi
  export CPPTRAJ_TEST_MODE
  export CPPTRAJ_TEST_ROOT=`pwd`
  # Ensure required binaries are set up
  if [ -z "$REMOVE" ] ; then
    # TODO is this being too paranoid?
    if [ ! -f '/bin/rm' ] ; then
      echo "Error: Required binary '/bin/rm' not found." > /dev/stderr
      exit 1
    fi
    export REMOVE='/bin/rm -f'
  fi
  Required "grep"
  Required "sed"
  Required "awk"
  # Process command line options
  CmdLineOpts $*
  # If not cleaning see what else needs to be set up. 
  if [ $CPPTRAJ_TEST_CLEAN -eq 0 ] ; then
    # Determine standalone or AmberTools
    if [ -z "$CPPTRAJ_STANDALONE" ] ; then # FIXME check needed?
      if [ ! -z "`echo "$CPPTRAJ_TEST_ROOT" | grep AmberTools`" ] ; then
        # Assume AmberTools. Need AMBERHOME.
        CPPTRAJ_STANDALONE=0
        if [ -z "$AMBERHOME" ] ; then
          echo "Error: In AmberTools and AMBERHOME is not set. Required for tests." > /dev/stderr
          exit 1
        fi
      else
        # Standalone. Never use dacdif.
        CPPTRAJ_STANDALONE=1
        USE_DACDIF=0
      fi
      export CPPTRAJ_STANDALONE
    fi
    # Determine if diff or dacdif will be used.
    CPPTRAJ_DIFF=''
    CPPTRAJ_DACDIF=''
    if [ $USE_DACDIF -eq 1 ] ; then
      CPPTRAJ_DACDIF="$AMBERHOME/test/dacdif"
      if [ ! -f "$CPPTRAJ_DACDIF" ] ; then
        echo "$CPPTRAJ_DACDIF not found. Required for AmberTools tests."
        exit 1
      fi
    else
      Required "diff"
      CPPTRAJ_DIFF=`which diff`
    fi
    export CPPTRAJ_DACDIF
    export CPPTRAJ_DIFF
    # Determine test output and error files: standalone only
    if [ -z "$CPPTRAJ_DACDIF" ] ; then
      if [ -z "$CPPTRAJ_TEST_RESULTS"  -o -z "$CPPTRAJ_TEST_ERROR" ] ; then
        export CPPTRAJ_TEST_RESULTS='Test_Results.dat'
        export CPPTRAJ_TEST_ERROR='Test_Error.dat'
      fi
    fi

    # Determine binary locations
    SetBinaries
  fi # END if not cleaning
  echo "DEBUG: Initial test setup complete."
  # If running multiple tests execute now.
  if [ "$CPPTRAJ_TEST_MODE" = 'multiple' ] ; then
    make test.test
  fi
fi # END initial setup

# Always clean up test output and error files
if [ "$CPPTRAJ_OUTPUT" != '/dev/stdout' -a -f "$CPPTRAJ_OUTPUT" ] ; then
  $REMOVE $CPPTRAJ_OUTPUT
fi
if [ "$CPPTRAJ_ERROR" != '/dev/stderr' -a -f "$CPPTRAJ_ERROR" ] ; then
  $REMOVE $CPPTRAJ_ERROR
fi

echo "DEBUG: Test mode is $CPPTRAJ_TEST_MODE"
