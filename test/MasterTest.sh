# This should be sourced at the top of CPPTRAJ test run scripts.

# Environment variables
# TEST_OS: Operating system on which tests are being run. If blank assume linux
# N_THREADS: Set to number of test threads
# DIFFOPTS: Additional options to pass to DIFF (non-DACDIF only)

# MasterTest.sh command line options
CLEAN=0             # If 1, only file cleaning needs to be performed.
SUMMARY=0           # If 1, only summary of results needs to be performed.
SHOWERRORS=0        # If 1, print test errors to STDOUT after summary.
STANDALONE=0        # If 0, part of AmberTools. If 1, stand-alone (e.g. from GitHub).
PROFILE=0           # If 1, end of test profiling with gprof performed
FORCE_AMBERTOOLS=0  # FIXME: currently needed to get extended tests to work
USEDACDIF=1         # If 0 do not use dacdif even if in AmberTools
CPPTRAJ=""          # CPPTRAJ binary
SFX=""              # CPPTRAJ binary suffix
AMBPDB=""           # ambpdb binary
NPROC=""            # nproc binary for counting threads in parallel tests.
TIME=""             # Set to the 'time' command if timing requested.
VALGRIND=""         # Set to 'valgrind' command if memory check requested.
DIFFCMD=""          # Command used to check for test differences
DACDIF=""           # Set if using 'dacdif' to test differences
NDIFF=""            # Set to Nelson H. F. Beebe's ndiff.awk for numerical diff calc
REMOVE="/bin/rm -f" # Remove command
NCDUMP=""           # ncdump command; needed for NcTest()
OUTPUT="test.out"   # File to direct test STDOUT to.
ERROR="/dev/stderr" # File to direct test STDERR to.
TEST_RESULTS=""     # For standalone, file to record test results to.
TEST_ERROR=""       # For standalone, file to record test errors to.
DEBUG=""            # Can be set to pass global debug flag to CPPTRAJ.
NUMTEST=0           # Total number of times DoTest has been called this test.
ERRCOUNT=0          # Total number of errors detected by DoTest this test.
WARNCOUNT=0         # Total number of warnings detected by DoTest this test.
PROGERROR=0         # Total number of program errors this test

# Options used in tests
TOP=""   # CPPTRAJ topology file/command line arg
INPUT="" # CPPTRAJ input file/command line arg

# Variables that describe how CPPTRAJ was compiled
ZLIB=""
BZLIB=""
NETCDFLIB=""
MPILIB=""
NOMATHLIB=""
OPENMP=""
PNETCDFLIB=""
SANDERLIB=""

# ------------------------------------------------------------------------------
# DoTest() <File1> <File2> [-r <relative err>] [-a <absolute err>] [<arg1>] ... [<argN>]
#   Compare File1 (the 'save' file) to File2 (test output), print an error if
#   they differ. If '-r' or '-a' are specified the test will pass as long as the
#   maximum relative or absolulte errors are below the given value (using Nelson
#   H. F. Beebe's ndiff.awk script. The remaining args can be used to pass
#   options to DIFFCMD.
DoTest() {
  echo "DEBUG: DoTest $1 $2"
  if [[ ! -z $DACDIF ]] ; then
    # AmberTools - use dacdif. Use any '-r <X>' or '-a <X>' args found.
    # Ignore the rest.
    DIFFARGS="$1 $2"
    shift # Save file
    shift # Test file
    # Process remaining args
    while [[ ! -z $1 ]] ; do
      case "$1" in
        "-r" ) shift ; DIFFARGS=" -r $1 "$DIFFARGS ;;
        "-a" ) shift ; DIFFARGS=" -a $1 "$DIFFARGS ;;
      esac
      shift
    done
    $DACDIF $DIFFARGS
  else
    # Standalone - will use diff, or ndiff where '-r' or '-a' specified.
    ((NUMTEST++))
    DIFFARGS="--strip-trailing-cr"
    NDIFFARGS=""
    # First two arguments are files to compare.
    F1=$1 ; shift
    F2=$1 ; shift
    # Process remaining arguments.
    USE_NDIFF=0
    while [[ ! -z $1 ]] ; do
      case "$1" in
        "-r"           ) USE_NDIFF=1; shift; NDIFFARGS="$NDIFFARGS -v RELERR=$1" ;;
        "-a"           ) USE_NDIFF=1; shift; NDIFFARGS="$NDIFFARGS -v ABSERR=$1" ;;
        *              ) DIFFARGS=$DIFFARGS" $1" ;;
      esac
      shift
    done
    if [[ ! -f "$F1" ]] ; then
      echo "  $F1 not found." >> $TEST_RESULTS
      echo "  $F1 not found." >> $TEST_ERROR
      ((ERRCOUNT++))
    elif [[ ! -f "$F2" ]] ; then
      echo "  $F2 not found." >> $TEST_RESULTS
      echo "  $F2 not found." >> $TEST_ERROR
      ((ERRCOUNT++))
    else
      if [[ $USE_NDIFF -eq 0 ]] ; then
        $DIFFCMD $DIFFARGS $DIFFOPTS $F1 $F2 > temp.diff 2>&1
      else
        $NDIFF $NDIFFARGS $F1 $F2 > temp.diff 2>&1
      fi
      if [[ -s temp.diff ]] ; then
        echo "  $F1 $F2 are different." >> $TEST_RESULTS
        echo "  $F1 $F2 are different." >> $TEST_ERROR
        cat temp.diff >> $TEST_ERROR
        ((ERRCOUNT++))
      else
        echo "  $F2 OK." >> $TEST_RESULTS
      fi
      $REMOVE temp.diff
    fi
  fi
}

# ------------------------------------------------------------------------------
# NcTest(): Compare NetCDF files <1> and <2>. Use NCDUMP to convert to ASCII
# first, removing ==> line and :programVersion attribute.
NcTest() {
  echo "DEBUG: NcTest $1 $2"
  if [[ -z $1 || -z $2 ]] ; then
    echo "Error: NcTest(): One or both files not specified." > /dev/stderr
    exit 1
  fi
  if [[ -z $NCDUMP || ! -e $NCDUMP ]] ; then
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
  while [[ ! -z $1 ]] ; do
    if [[ $1 = '-r' || $1 = '-a' ]] ; then
      CALC_NUM_ERR=1
    fi
    DIFFARGS=$DIFFARGS" $1"
    shift
  done
  # Prepare files.
  if [[ ! -e $F1 ]] ; then
    echo "Error: $F1 missing." >> $TEST_ERROR
  elif [[ ! -e $F2 ]] ; then
    echo "Error: $F2 missing." >> $TEST_ERROR
  else
    if [[ $CALC_NUM_ERR -eq 1 ]] ; then
      # FIXME: Must remove commas here because I cannot figure out how to pass
      # the regular expression to ndiff.awk FS without the interpreter giving
      # this error for FS='[ \t,()]':
      # awk: fatal: Invalid regular expression: /'[/
      $NCDUMP -n nctest $F1 | grep -v "==>\|:programVersion" | sed 's/,/ /g' > nc0.save
      $NCDUMP -n nctest $F2 | grep -v "==>\|:programVersion" | sed 's/,/ /g' > nc0
    else
      $NCDUMP -n nctest $F1 | grep -v "==>\|:programVersion" > nc0.save
      $NCDUMP -n nctest $F2 | grep -v "==>\|:programVersion" > nc0
    fi
    DoTest $DIFFARGS 
    $REMOVE nc0.save nc0
  fi
}

# ------------------------------------------------------------------------------
# CheckTest(): Report if the error counter is greater than 0. TODO Remove
CheckTest() {
  # Only use when not using dacdif 
  if [[ -z $DACDIF ]] ; then
    if [[ $ERRCOUNT -gt 0 ]] ; then
      echo "  $ERRCOUNT comparisons failed so far."
    fi
  fi
}

# ------------------------------------------------------------------------------
# RunCpptraj(): Run cpptraj with the given options.
RunCpptraj() {
  # If only cleaning requested no run needed, exit now
  if [[ $CLEAN -eq 1 ]] ; then
    exit 0
  fi
  echo ""
  echo "  CPPTRAJ: $1"
  if [[ -z $DACDIF ]] ; then
    echo "  CPPTRAJ: $1" >> $TEST_RESULTS
  fi
#  if [[ ! -z $DEBUG ]] ; then
    echo "$TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $DEBUG $TOP $INPUT >> $OUTPUT 2>>$ERROR"
#  fi
  $TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $DEBUG $TOP $INPUT >> $OUTPUT 2>>$ERROR
  STATUS=$?
  echo "DEBUG: Cpptraj exited with status $STATUS"
  if [ "$STATUS" -ne 0 ] ; then
    echo "Error: cpptraj exited with status $STATUS" 2> /dev/stderr
    echo "Error: cpptraj exited with status $STATUS" > $TEST_RESULTS
  fi
}

# ------------------------------------------------------------------------------
# EndTest(): Called at the end of every test script if no errors found.
EndTest() {
  echo "DEBUG: EndTest"
  # Report only when not using dacdif 
  if [[ -z $DACDIF ]] ; then
    if [[ $ERRCOUNT -gt 0 ]] ; then
      echo "  $ERRCOUNT out of $NUMTEST comparisons failed."
      echo "  $ERRCOUNT out of $NUMTEST comparisons failed." >> $TEST_RESULTS
      echo "  $ERRCOUNT out of $NUMTEST comparisons failed." >> $TEST_ERROR
    elif [[ $WARNCOUNT -gt 0 ]] ; then
      ((PASSCOUNT = $NUMTEST - $WARNCOUNT))
      echo "  $PASSCOUNT out of $NUMTEST passing comparisons. $WARNCOUNT warnings."
      echo "  $PASSCOUNT out of $NUMTEST passing comparisons. $WARNCOUNT warnings." >> $TEST_RESULTS
    else 
      echo "All $NUMTEST comparisons passed." 
      echo "All $NUMTEST comparisons passed." >> $TEST_RESULTS 
    fi
    echo ""
    if [[ ! -z $VALGRIND ]] ; then
      echo "Valgrind summary:"
      grep ERROR $ERROR
      grep heap $ERROR
      grep LEAK $ERROR
      echo ""
      echo "Valgrind summary:" >> $TEST_RESULTS
      grep ERROR $ERROR >> $TEST_RESULTS
      grep heap $ERROR >> $TEST_RESULTS
      grep LEAK $ERROR >> $TEST_RESULTS

    fi
    if [[ $PROFILE -eq 1 ]] ; then
      if [[ -e gmon.out ]] ; then
        gprof $CPPTRAJ > profiledata.txt
      fi
    fi
  fi
}

# ------------------------------------------------------------------------------
# CleanFiles(): For every arg passed to the function, check for the file and rm it
CleanFiles() {
  while [[ ! -z $1 ]] ; do
    if [[ -e $1 ]] ; then
      $REMOVE $1
    fi
    shift
  done
  # If only cleaning requested no run needed, exit now
  if [[ $CLEAN -eq 1 ]] ; then
    exit 0
  fi
}

# ------------------------------------------------------------------------------
# Library Checks - Tests that depend on certain libraries like Zlib can run
# these to make sure cpptraj was compiled with that library - exit gracefully
# if not.
# Should not be called if CLEAN==1, CleanFiles should always be called first.
CheckZlib() {
  if [[ -z $ZLIB ]] ; then
    echo "This test requires zlib. Cpptraj was compiled without zlib support."
    echo "Skipping test."
    exit 0
  fi
}

CheckBzlib() {
  if [[ -z $BZLIB ]] ; then
    echo "This test requires bzlib. Cpptraj was compiled without bzlib support."
    echo "Skipping test."
    exit 0
  fi
}

CheckNetcdf() {
  if [[ -z $NETCDFLIB ]] ; then
    echo "This test requires Netcdf. Cpptraj was compiled without Netcdf support."
    echo "Skipping test."
    exit 0
  fi
}

CheckPtrajAnalyze() {
  if [[ ! -z $NOMATHLIB ]] ; then
    echo "This test requires LAPACK/ARPACK/BLAS routines from AmberTools."
    echo "Cpptraj was compiled with -DNO_MATHLIB. Skipping test."
    exit 0
  fi
}

CheckSanderlib() {
  if [[ -z $SANDERLIB ]] ; then
    if [[ ! -z $1 ]] ; then
      DESCRIP=$1
    else
      DESCRIP="This test"
    fi
    echo "$DESCRIP requires compilation with the Sander API from AmberTools."
    echo "Skipping test."
    return 1
  fi
  return 0
}

CheckPnetcdf() {
  if [[ ! -z $DO_PARALLEL ]] ; then
    DESCRIP="This test"
    if [[ ! -z $1 ]] ; then
      DESCRIP="Test '$1'"
    fi
    if [[ -z $PNETCDFLIB ]] ; then
      echo "$DESCRIP requires compilation with Pnetcdf."
      echo "Cpptraj was compiled without Pnetcdf support. Skipping test."
      return 1
    fi
  fi
  return 0
}

# NotParallel() <Test title>
NotParallel() {
  if [[ ! -z $DO_PARALLEL ]] ; then
    echo ""
    echo "  CPPTRAJ: $1"
    echo "  This test cannot be run in parallel. Skipping test."
    return 1
 fi
 return 0
}

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

# RequiresThreads() <# threads> <Test title>
RequiresThreads() {
  if [[ ! -z $DO_PARALLEL ]] ; then
    SetNthreads
    if [ $? -ne 0 ] ; then
      echo "Error: Program to find # threads not found ($NPROC)" > /dev/stderr
      echo "Error: Test requires $1 parallel threads. Attempting to run test anyway." > /dev/stderr
      return 0
    fi
    REMAINDER=`echo "$N_THREADS % $1" | bc`
    if [[ -z $REMAINDER || $REMAINDER -ne 0 ]] ; then
      echo ""
      if [[ ! -z $2 ]] ; then
        echo "  CPPTRAJ: $2"
      fi
      echo "  Warning: Test requires a multiple of $1 parallel threads. Skipping."
      return 1
    fi
  fi
  return 0
}

# MaxThreads() <# threads> <Test title>
MaxThreads() {
  if [[ ! -z $DO_PARALLEL ]] ; then
    SetNthreads
    if [ $? -ne 0 ] ; then
      echo "Error: Program to find # threads not found ($NPROC)" > /dev/stderr
      echo "Error: Test can only run with $1 or fewer threads. Attempting to run test anyway." > /dev/stderr
      return 0
    fi
    if [[ $N_THREADS -gt $1 ]] ; then
      echo ""
      if [[ ! -z $2 ]] ; then
        echo "  CPPTRAJ: $2"
      fi
      echo "  Warning: Test can only run with $1 or fewer parallel threads. Skipping."
      return 1
    fi
  fi
  return 0
}
#-------------------------------------------------------------------------------
# Summary(): Print a summary of the tests.
Summary() {
  RESULTFILES=""
  if [[ ! -z $TEST_RESULTS ]] ; then
    RESULTFILES=`ls */$TEST_RESULTS 2> /dev/null`
  else
    exit 0
  fi
  echo "===================== TEST SUMMARY ======================"
  if [[ ! -z $RESULTFILES ]] ; then
    cat $RESULTFILES > $TEST_RESULTS
    # DoTest - Number of comparisons OK
    OK=`grep OK $TEST_RESULTS | wc -l`
    # DoTest - Number of warnings
    WARN=`grep Warning $TEST_RESULTS | wc -l`
    # DoTest - Number of comparisons different
    ERR=`grep different $TEST_RESULTS | wc -l`
    NOTFOUND=`grep "not found" $TEST_RESULTS | wc -l`
    PERR=`grep "Error:" $TEST_RESULTS | wc -l`
    ((ERR = $ERR + $NOTFOUND + $PERR))
    # Number of tests run
    NTESTS=`grep "TEST:" $TEST_RESULTS | wc -l`
    # Number of tests successfully finished
    PASSED=`grep "comparisons passed" $TEST_RESULTS | wc -l`
    ((NCOMPS = $OK + $ERR + $WARN))
    echo "  $OK out of $NCOMPS comparisons OK ($ERR failed, $WARN warnings)."
    echo "  $PASSED out of $NTESTS tests completed with no issues."
    RESULTFILES=`ls */$TEST_ERROR 2> /dev/null`
    if [[ ! -z $RESULTFILES ]] ; then
      cat $RESULTFILES > $TEST_ERROR
    fi 
  else
    echo "No Test Results files (./*/$TEST_RESULTS) found."
  fi

  if [[ $SHOWERRORS -eq 1 && $ERR -gt 0 ]]; then
    echo "Obtained the following errors:"
    echo "---------------------------------------------------------"
    cat $TEST_ERROR
    echo "---------------------------------------------------------"
  fi

  if [[ ! -z $VALGRIND ]] ; then
    RESULTFILES=`ls */$ERROR 2> /dev/null`
    if [[ ! -z $RESULTFILES ]] ; then
      echo "---------------------------------------------------------"
      echo "Valgrind summary:"
      NUMVGERR=`grep ERROR $RESULTFILES | awk 'BEGIN{sum=0;}{sum+=$4;}END{print sum;}'`
      echo "    $NUMVGERR errors."
      NUMVGOK=`grep "All heap" $RESULTFILES | wc -l`
      echo "    $NUMVGOK memory leak checks OK."
      NUMVGLEAK=`grep LEAK $RESULTFILES | wc -l`
      echo "    $NUMVGLEAK memory leak reports."
    else
      echo "No valgrind test results found."
      exit $ERR
    fi
  fi
  echo "========================================================="
  exit $ERR
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
# CmdLineOpts(): Process test script command line options
CmdLineOpts() {
  VGMODE=0 # Valgrind mode: 0 none, 1 memcheck, 2 helgrind
  SFX_OMP=0
  SFX_CUDA=0
  SFX_MPI=0
  while [[ ! -z $1 ]] ; do
    case "$1" in
      "summary"   ) SUMMARY=1 ;;
      "showerrors") SHOWERRORS=1 ;;
      "stdout"    ) OUTPUT="/dev/stdout" ;;
      "openmp"    ) SFX_OMP=1 ;;
      "cuda"      ) SFX_CUDA=1 ;;
      "mpi"       ) SFX_MPI=1 ;;
      "vg"        ) VGMODE=1 ;;
      "vgh"       ) VGMODE=2 ;;
      "time"      ) TIME=`which time` ;;
      "-at"       ) FORCE_AMBERTOOLS=1 ;;
      "-d"        ) DEBUG="-debug 4" ;;
      "-debug"    ) shift ; DEBUG="-debug $1" ;;
      "-nodacdif" ) USEDACDIF=0 ;;
      "-cpptraj"  ) shift ; CPPTRAJ=$1 ; echo "Using cpptraj: $CPPTRAJ" ;;
      "-ambpdb"   ) shift ; AMBPDB=$1  ; echo "Using ambpdb: $AMBPDB" ;;
      "-profile"  ) PROFILE=1 ; echo "Performing gnu profiling during EndTest." ;;
      "-h" | "--help" ) Help ; exit 0 ;;
      *           ) echo "Error: Unknown opt: $1" > /dev/stderr ; exit 1 ;;
    esac
    shift
  done
  # If DO_PARALLEL has been set force MPI
  if [ ! -z "$DO_PARALLEL" ] ; then
    SFX_MPI=1
    MPI=1
  fi
  # Warn if using OpenMP but OMP_NUM_THREADS not set.
  if [ "$SFX_OMP" -eq 1 -a -z "$OMP_NUM_THREADS" ] ; then
    echo "Warning: Using OpenMP but OMP_NUM_THREADS is not set."
  fi
  # Set up SFX
  if [ "$SFX_OMP"  -eq 1 ] ; then SFX="$SFX.OMP"  ; fi
  if [ "$SFX_CUDA" -eq 1 ] ; then SFX="$SFX.cuda" ; fi
  if [ "$SFX_MPI"  -eq 1 ] ; then SFX="$SFX.MPI"  ; fi
  # Set up valgrind if necessary
  if [[ $VGMODE -ne 0 ]] ; then
    VG=`which valgrind`
    if [[ -z $VG ]] ; then
      echo "Error: Valgrind not found." > /dev/stderr
      echo "Error:    Make sure valgrind is installed and in your PATH" > /dev/stderr
      exit 1
    fi
    echo "  Using Valgrind."
    ERROR="valgrind.out"
    if [[ $VGMODE -eq 1 ]] ; then
      VALGRIND="valgrind --tool=memcheck --leak-check=yes --show-reachable=yes"
    elif [[ $VGMODE -eq 2 ]] ; then
      VALGRIND="valgrind --tool=helgrind"
    fi
  fi
  # Figure out if we are a part of AmberTools
  if [[ -z $CPPTRAJ ]] ; then
    if [[ ! -z `pwd | grep AmberTools` || $FORCE_AMBERTOOLS -eq 1 ]] ; then
      STANDALONE=0
    else
      STANDALONE=1
    fi
  else
    # CPPTRAJ was specified. Assume standalone.
    STANDALONE=1
  fi
}

#-------------------------------------------------------------------------------
# SetBinaries(): Set and check CPPTRAJ etc binaries
SetBinaries() {
  # Set default command locations
  DIFFCMD=`which diff`
  NCDUMP=`which ncdump`
  # Set locations for CPPTRAJ, AMBPDB, NPROC, NDIFF if not already set.
  if [ "$STANDALONE" -eq 0 ] ; then
    # AmberTools
    if [ -z "$AMBERHOME" ] ; then
      echo "Warning: AMBERHOME is not set."
      # Assume we are running in $AMBERHOME/AmberTools/src/test/Test_X
      DIRPREFIX=../../../../
    else
      DIRPREFIX=$AMBERHOME
    fi
    # If not using DACDIF, NDIFF needs to be set.
    if [ "$USEDACDIF" -eq 1 ] ; then
      DACDIF=$DIRPREFIX/test/dacdif
      if [[ ! -f "$DACDIF" ]] ; then
        echo "Error: dacdiff command '$DACDIFF' not found." > /dev/stderr
        echo "Error: dacdiff command '$DACDIFF' not found." > $TEST_RESULTS
        exit 1
      fi
    else
      NDIFF=$DIRPREFIX/test/ndiff.awk
    fi
    if [ -f "$DIRPREFIX/bin/ncdump" ] ; then
      NCDUMP=$DIRPREFIX/bin/ncdump
    fi
    if [ -z "$CPPTRAJ" ] ; then CPPTRAJ=$DIRPREFIX/bin/cpptraj$SFX ; fi
    if [ -z "$AMBPDB"  ] ; then AMBPDB=$DIRPREFIX/bin/ambpdb ; fi
    NPROC=$DIRPREFIX/AmberTools/test/numprocs
  else
    # Standalone: GitHub etc
    if [ -z "$CPPTRAJHOME" ] ; then
      echo "Warning: CPPTRAJHOME is not set."
      DIRPREFIX=../../
    else
      DIRPREFIX=$CPPTRAJHOME
    fi
    if [ -z "$CPPTRAJ" ] ; then CPPTRAJ=$DIRPREFIX/bin/cpptraj$SFX ; fi
    if [ -z "$AMBPDB"  ] ; then AMBPDB=$DIRPREFIX/bin/ambpdb ; fi
    NPROC=$DIRPREFIX/test/nproc
    NDIFF=$DIRPREFIX/util/ndiff/ndiff.awk
    if [ ! -f "$NDIFF" ] ; then
      echo "Error: 'ndiff.awk' not present in cpptraj: $NDIFF"
      exit 1
    fi
  fi
  if [ ! -z "$NDIFF" ] ; then
    NDIFF="awk -f $NDIFF"
  fi
  # Check binaries
  if [ ! -f "$NCDUMP" ] ; then
    echo "Warning: 'ncdump' not found; NetCDF file comparisons cannot be performed."
  fi
  if [ ! -f "$DIFFCMD" ] ; then
    echo "Error: diff command '$DIFFCMD' not found." > /dev/stderr
    echo "Error: diff command '$DIFFCMD' not found." > $TEST_RESULTS
    exit 1
  fi
  if [ ! -f "$CPPTRAJ" ] ; then
    echo "Error: cpptraj binary '$CPPTRAJ' not found." > /dev/stderr
    echo "Error: cpptraj binary '$CPPTRAJ' not found." > $TEST_RESULTS
    exit 1
  fi
  if [ ! -z "$DEBUG" -o -z "$DACDIF" ] ; then
    ls -l $CPPTRAJ
  fi
  if [ ! -f "$AMBPDB" ] ; then
    # Try to locate it based on the location of CPPTRAJ
    DIRPREFIX=`dirname $CPPTRAJ`
    AMBPDB=$DIRPREFIX/ambpdb
    if [ ! -f "$AMBPDB" ] ; then
      echo "Warning: AMBPDB not present."
      AMBPDB=""
    fi
  fi
  # Print DEBUG info
  if [ ! -z "$DEBUG" ] ; then
    if [ $STANDALONE -eq 1 ] ; then
      echo "DEBUG: Standalone mode."
    else
      echo "DEBUG: AmberTools mode."
    fi
    echo "DEBUG: CPPTRAJ: $CPPTRAJ"
    echo "DEBUG: AMBPDB:  $AMBPDB"
    echo "DEBUG: NPROC:   $NPROC"
    echo "DEBUG: NCDUMP:  $NCDUMP"
    echo "DEBUG: DIFFCMD: $DIFFCMD"
    echo "DEBUG: DACDIF:  $DACDIF"
    echo "DEBUG: NDIFF:   $NDIFF"
  fi
}

#-------------------------------------------------------------------------------
# CheckDefines(): Check how CPPTRAJ was compiled.
CheckDefines() {
  ZLIB=''
  BZLIB=''
  NETCDFLIB=''
  MPILIB=''
  NOMATHLIB=''
  OPENMP=''
  PNETCDFLIB=''
  SANDERLIB=''
  CUDA=''
  NO_XDRFILE=''
  for DEFINE in `$CPPTRAJ --defines` ; do
    case "$DEFINE" in
      '-DHASGZ'         ) ZLIB=$DEFINE ;;
      '-DHASBZ2'        ) BZLIB=$DEFINE ;;
      '-DBINTRAJ'       ) NETCDFLIB=$DEFINE ;;
      '-DMPI'           ) MPILIB=$DEFINE ;;
      '-DNO_MATHLIB'    ) NOMATHLIB=$DEFINE ;;
      '-D_OPENMP'       ) OPENMP=$DEFINE ;;
      '-DHAS_PNETCDF'   ) PNETCDFLIB=$DEFINE ;;
      '-DUSE_SANDERLIB' ) SANDERLIB=$DEFINE ;;
      '-DCUDA'          ) CUDA=$DEFINE ;;
      '-DNO_XDRFILE'    ) NO_XDRFILE=$DEFINE ;;
    esac
  done
  #echo "DEBUG: $ZLIB $BZLIB $NETCDFLIB $MPILIB $NOMATHLIB $OPENMP $PNETCDFLIB $SANDERLIB $CUDA $NO_XDRFILE"
}

#===============================================================================
# If the first argument is "clean" then no set-up is required. Script will
# exit when either CleanFiles or RunCpptraj is called from sourcing script.
if [[ $1 = "clean" ]] ; then
  CLEAN=1
else
  CmdLineOpts $*
  # Set results files if not using dacdif 
  if [[ -z $DACDIF ]] ; then
    TEST_RESULTS=Test_Results.dat
    TEST_ERROR=Test_Error.dat
    if [[ -f "$TEST_RESULTS" ]] ; then
      $REMOVE $TEST_RESULTS
    fi
    if [[ -f "$TEST_ERROR" ]] ; then
      $REMOVE $TEST_ERROR
    fi
    if [ ! -z "$DIFFOPTS" ] ; then
      echo "Warning: DIFFOPTS is set to '$DIFFOPTS'"
    fi
  fi
  # Only a summary of previous results has been requested.
  if [[ $SUMMARY -eq 1 ]] ; then
    Summary
  fi
  # If TEST_OS is not set, assume linux
  if [[ -z $TEST_OS ]] ; then
    TEST_OS="linux"
  fi
  # Set binary locations
  SetBinaries
  # Check how CPPTRAJ was compiled
  CheckDefines
  # If CUDA, use cuda-memcheck instead of valgrind
  if [[ ! -z $VALGRIND && ! -z $CUDA ]] ; then
    VALGRIND='cuda-memcheck'
  fi
  # Start test results file
  echo "**************************************************************"
  echo "TEST: `pwd`"
  if [[ -z $DACDIF ]] ; then
    echo "**************************************************************" > $TEST_RESULTS
    echo "TEST: `pwd`" >> $TEST_RESULTS
  fi
fi
# Always clean up OUTPUT and ERROR
if [[ $OUTPUT != "/dev/stdout" && -f "$OUTPUT" ]] ; then
  $REMOVE $OUTPUT
fi
if [[ $ERROR != "/dev/stderr" && -f "$ERROR" ]] ; then
  $REMOVE $ERROR
fi
