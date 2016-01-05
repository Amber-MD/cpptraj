# This should be sourced at the top of CPPTRAJ test run scripts.

# Environment variables
# TEST_OS: Operating system on which tests are being run. If blank assume linux

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

# ------------------------------------------------------------------------------
# DoTest() <File1> <File2> [allowfail <OS>] [<arg1>] ... [<argN>]
#   Compare File1 (the 'save' file) to File2 (test output), print an error if
#   they differ. The 'allowfail' keyword allows test to fail with just a
#   warning, and is currently intended for tests with known failures on given
#   <OS>. The remaining args can be used to pass options to DIFFCMD.
DoTest() {
  if [[ ! -z $DACDIF ]] ; then
    # AmberTools - will use dacdif.
    $DACDIF $1 $2
  else
    # Standalone - will use diff.
    ((NUMTEST++))
    DIFFARGS="--strip-trailing-cr"
    # First two arguments are files to compare.
    F1=$1 ; shift
    F2=$1 ; shift
    # Process remaining arguments.
    ALLOW_FAIL=0
    FAIL_OS=""
    while [[ ! -z $1 ]] ; do
      case "$1" in
        "allowfail" ) ALLOW_FAIL=1 ; shift ; FAIL_OS=$1 ;;
        *           ) DIFFARGS=$DIFFARGS" $1" ;;
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
      $DIFFCMD $DIFFARGS $F1 $F2 > temp.diff 2>&1
      if [[ -s temp.diff ]] ; then
        if [[ $ALLOW_FAIL -eq 1 && $TEST_OS = $FAIL_OS ]] ; then
          echo "  Warning: Differences between $F1 and $F2 detected."
          echo "  Warning: Differences between $F1 and $F2 detected." >> $TEST_RESULTS
          echo "           This test is known to fail on $TEST_OS."
          echo "           This test is known to fail on $TEST_OS." >> $TEST_RESULTS
          echo "           The differences below should be carefully inspected."
          echo "--------------------------------------------------------------------------------"
          cat temp.diff
          echo "--------------------------------------------------------------------------------"
          ((WARNCOUNT++))
        else
          echo "  $F1 $F2 are different." >> $TEST_RESULTS
          echo "  $F1 $F2 are different." >> $TEST_ERROR
          cat temp.diff >> $TEST_ERROR
          ((ERRCOUNT++))
        fi
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
  if [[ -z $1 || -z $2 ]] ; then
    echo "Error: NcTest(): One or both files not specified." > /dev/stderr
    exit 1
  fi
  if [[ -z $NCDUMP || ! -e $NCDUMP ]] ; then
    echo "ncdump missing." > /dev/stderr
    exit 1
  fi
  # Prepare files.
  if [[ ! -e $1 ]] ; then
    echo "Error: $1 missing." >> $TEST_ERROR
  elif [[ ! -e $2 ]] ; then
    echo "Error: $2 missing." >> $TEST_ERROR
  else
    $NCDUMP -n nctest $1 | grep -v "==>\|:programVersion" > nc0
    $NCDUMP -n nctest $2 | grep -v "==>\|:programVersion" > nc1
    DoTest nc0 nc1
    $REMOVE nc0 nc1
  fi
}

# ------------------------------------------------------------------------------
# CheckTest(): Report if the error counter is greater than 0. TODO Remove
CheckTest() {
  # Only use when not using dacdif 
  if [[ -z $DACDIF ]] ; then
    if [[ $ERR -gt 0 ]] ; then
      echo "  $ERR comparisons failed so far."
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
  if [[ ! -z $DEBUG ]] ; then
    echo "$TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $DEBUG $TOP $INPUT >> $OUTPUT 2>>$ERROR"
  fi
  $TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $DEBUG $TOP $INPUT >> $OUTPUT 2>>$ERROR
}

# ------------------------------------------------------------------------------
# EndTest(): Called at the end of every test script if no errors found.
EndTest() {
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

CheckPnetcdf() {
  DESCRIP="This test"
  if [[ ! -z $1 ]] ; then
    DESCRIP="Test '$1'"
  fi
  if [[ -z $PNETCDFLIB ]] ; then
    echo "$DESCRIP requires compilation with Pnetcdf."
    echo "Cpptraj was compiled without Pnetcdf support. Skipping test."
    return 1
  fi
  return 0
}

NotParallel() {
  if [[ ! -z $DO_PARALLEL ]] ; then
    echo ""
    echo "  CPPTRAJ: $1"
    echo "  This test cannot be run in parallel. Skipping test."
    return 1
 fi
 return 0
}

RequiresThreads() {
  if [[ ! -z $DO_PARALLEL ]] ; then
    if [[ ! -f "$NPROC" ]] ; then
      echo "Error: Program to find # threads not found ($NPROC)" > /dev/stderr
      echo "Error: Test requires $1 threads. Attempting to run test anyway." > /dev/stderr
      return 0
    fi
    N_THREADS=`$DO_PARALLEL $NPROC`
    if [[ $N_THREADS -ne $1 ]] ; then
      echo ""
      echo "Warning: Test requires $1 threads. Skipping."
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
    OK=`cat $TEST_RESULTS | grep OK | wc -l`
    # DoTest - Number of warnings
    WARN=`cat $TEST_RESULTS | grep Warning | wc -l`
    # DoTest - Number of comparisons different
    ERR=`cat $TEST_RESULTS | grep different | wc -l`
    NOTFOUND=`cat $TEST_RESULTS | grep "not found" | wc -l`
    ((ERR = $ERR + $NOTFOUND))
    # Number of tests run
    NTESTS=`cat $TEST_RESULTS | grep "TEST:" | wc -l`
    # Number of tests successfully finished
    PASSED=`cat $TEST_RESULTS | grep "comparisons passed" | wc -l`
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
      NUMVGERR=`cat $RESULTFILES | grep ERROR | awk 'BEGIN{sum=0;}{sum+=$4;}END{print sum;}'`
      echo "    $NUMVGERR errors."
      NUMVGOK=`cat $RESULTFILES | grep "All heap" | wc -l`
      echo "    $NUMVGOK memory leak checks OK."
      NUMVGLEAK=`cat $RESULTFILES | grep LEAK | wc -l`
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
  echo "  mpi        : Use MPI version of CPPTRAJ."
  echo "  openmp     : Use OpenMP version of CPPTRAJ."
  echo "  vg         : Run test with valgrind memcheck."
  echo "  vgh        : Run test with valgrind helgrind."
  echo "  time       : Time the test."
  echo "  -at        : Force AmberTools tests."
  echo "  -nodacdif  : Do not use dacdif for test comparisons."
  echo "  -d         : Run CPPTRAJ with global debug level 4."
  echo "  -debug <#> : Run CPPTRAJ with global debug level #."
  echo "  -cpptraj <file> : Use CPPTRAJ binary <file>."
  echo "  -ambpdb <file>  : Use AMBPDB binary <file>."
  echo "  -profile        : Profile results with 'gprof' (requires special compile)."
}

#-------------------------------------------------------------------------------
# CmdLineOpts(): Process test script command line options
CmdLineOpts() {
  VGMODE=0 # Valgrind mode: 0 none, 1 memcheck, 2 helgrind
  while [[ ! -z $1 ]] ; do
    case "$1" in
      "summary"   ) SUMMARY=1 ;;
      "showerrors") SHOWERRORS=1 ;;
      "stdout"    ) OUTPUT="/dev/stdout" ;;
      "mpi"       ) SFX=".MPI" ;;
      "openmp"    ) SFX=".OMP" ;;
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
  # If DO_PARALLEL has been set force MPI
  if [[ ! -z $DO_PARALLEL ]] ; then
    SFX=".MPI"
    MPI=1
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
  # Set CPPTRAJ binary location if not already set.
  if [[ -z $CPPTRAJ ]] ; then
    if [[ $STANDALONE -eq 0 ]] ; then
      # AmberTools
      if [[ -z $AMBERHOME ]] ; then
        echo "Warning: AMBERHOME is not set."
        # Assume we are running in $AMBERHOME/AmberTools/src/test/Test_X
        DIRPREFIX=../../../../
      else
        DIRPREFIX=$AMBERHOME
      fi
      if [[ $USEDACDIF -eq 1 ]] ; then
        DACDIF=$DIRPREFIX/test/dacdif
      fi
      NCDUMP=$DIRPREFIX/bin/ncdump
      CPPTRAJ=$DIRPREFIX/bin/cpptraj$SFX
      AMBPDB=$DIRPREFIX/bin/ambpdb
      NPROC=$DIRPREFIX/AmberTools/test/numprocs
    else
      # Standalone: GitHub etc
      if [[ ! -z $CPPTRAJHOME ]] ; then
        CPPTRAJ=$CPPTRAJHOME/bin/cpptraj$SFX
        AMBPDB=$CPPTRAJHOME/bin/ambpdb
        NPROC=$CPPTRAJHOME/test/nproc
      else
        CPPTRAJ=../../bin/cpptraj$SFX
        AMBPDB=../../bin/ambpdb
        NPROC=../../test/nproc
      fi
    fi
  fi
  # Print DEBUG info
  if [[ ! -z $DEBUG ]] ; then
    if [[ $STANDALONE -eq 1 ]] ; then
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
  fi
  # Check binaries
  if [[ ! -f "$NCDUMP" ]] ; then
    echo "Warning: 'ncdump' not found; NetCDF file comparisons cannot be performed."
  fi
  if [[ ! -f "$DIFFCMD" ]] ; then
    echo "Error: diff command '$DIFFCMD' not found." > /dev/stderr
    exit 1
  fi
  if [[ $STANDALONE -eq 0 && $USEDACDIF -eq 1 && ! -f "$DACDIF" ]] ; then
    echo "Error: dacdiff command '$DACDIFF' not found." > /dev/stderr
    exit 1
  fi
  if [[ ! -f "$CPPTRAJ" ]] ; then
    echo "Error: cpptraj binary '$CPPTRAJ' not found." > /dev/stderr
    exit 1
  fi
  if [[ ! -z $DEBUG || -z $DACDIF ]] ; then
    ls -l $CPPTRAJ
  fi
  if [[ ! -f "$AMBPDB" ]] ; then
    # Try to locate it based on the location of CPPTRAJ
    DIRPREFIX=`dirname $CPPTRAJ`
    AMBPDB=$DIRPREFIX/ambpdb
    if [[ ! -f "$AMBPDB" ]] ; then
      echo "Warning: ambpdb binary '$AMBPDB' not found."
    fi
  fi
}

#-------------------------------------------------------------------------------
# CheckDefines(): Check how CPPTRAJ was compiled.
CheckDefines() {
  DEFINES=`$CPPTRAJ --defines | grep Compiled`
  ZLIB=`echo $DEFINES | grep DHASGZ`
  BZLIB=`echo $DEFINES | grep DHASBZ2`
  NETCDFLIB=`echo $DEFINES | grep DBINTRAJ`
  MPILIB=`echo $DEFINES | grep DMPI`
  NOMATHLIB=`echo $DEFINES | grep DNO_MATHLIB`
  OPENMP=`echo $DEFINES | grep D_OPENMP`
  PNETCDFLIB=`echo $DEFINES | grep DHAS_PNETCDF`
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
