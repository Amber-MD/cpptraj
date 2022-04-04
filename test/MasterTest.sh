# This should be sourced at the top of CPPTRAJ test run scripts.
#
# The purpose of this script is to provide all common functionality required
# by all tests, including environment setup and collecting test results.
# Requires awk, grep, sed, diff (if not AmberTools), rm, make (multiple tests)
#
# ----- Environment variables ------------------------------
# Binary locations
#   CPPTRAJ              : Set to cpptraj binary being tested.
#   VALGRIND             : Set to 'valgrind' command if memory check requested.
#   CPPTRAJ_DIFF         : Command used to check for test differences.
#   CPPTRAJ_DACDIF       : Set if testing inside AmberTools.
#   CPPTRAJ_NDIFF        : Set to Nelson H. F. Beebe's ndiff.awk for numerical diff calc.
#   CPPTRAJ_NCDUMP       : Set to ncdump command; needed for NcTest()
#   CPPTRAJ_RM           : Command used to remove files
#   CPPTRAJ_TIME         : Set to the 'time' command if timing requested.
#   CPPTRAJ_NPROC        : nproc binary for counting threads in parallel tests.
# Test output locations
#   CPPTRAJ_TEST_RESULTS : File to record individual test results to.
#   CPPTRAJ_TEST_ERROR   : File to record individual test errors/diffs to.
#   CPPTRAJ_OUTPUT       : File to direct cpptraj STDOUT to.
#   CPPTRAJ_ERROR        : File to direct cpptraj STDERR to.
# Other variables
#   CPPTRAJ_TEST_ROOT    : Test root directory. FIXME make local?
#   CPPTRAJ_TEST_SETUP   : 'yes' if setup is complete.
#   CPPTRAJ_TEST_CLEAN   : If 1, only cleaning tests; do not run them.
#   CPPTRAJ_TEST_OS      : Operating system on which tests are being run. If blank assume linux.
#   N_THREADS            : Number of MPI threads if parallel.
#   OMP_NUM_THREADS      : Max number of OpenMP threads.
#   DO_PARALLEL          : MPI run command (e.g. 'mpirun -n 11')
#   CPPTRAJ_DEBUG        : Can be set to pass global debug flag to cpptraj.
#   DIFFOPTS             : Additional options to pass to CPPTRAJ_DIFF
#   CPPTRAJ_PROFILE      : If 1, end of test profiling with gprof performed.
#   CPPTRAJ_LONG_TEST    : If 1, enable long tests
#   CPPTRAJ_RNG          : Hold the --rng arg for CPPTRAJ.
# Cpptraj binary characteristics
#   CPPTRAJ_ZLIB         : If set CPPTRAJ has zlib support.
#   CPPTRAJ_BZLIB        : If set CPPTRAJ has bzip support.
#   CPPTRAJ_NETCDFLIB    : If set CPPTRAJ has NetCDF support.
#   CPPTRAJ_C11_SUPPORT  : If set CPPTRAJ has C++11 support.
#   CPPTRAJ_LIBPME       : If set CPPTRAJ was compiled with libPME.
#   CPPTRAJ_MPILIB       : If set CPPTRAJ has MPI support.
#   CPPTRAJ_MATHLIB      : If set CPPTRAJ was compiled with math libraries.
#   CPPTRAJ_OPENMP       : If set CPPTRAJ has OpenMP support.
#   CPPTRAJ_PNETCDFLIB   : If set CPPTRAJ has parallel NetCDF support.
#   CPPTRAJ_SANDERLIB    : If set CPPTRAJ was compiled with the sander API from AT.
#   CPPTRAJ_FFTW_FFT     : If set CPPTRAJ was compiled with fftw
#   CPPTRAJ_CUDA         : If set CPPTRAJ has CUDA support.
#   CPPTRAJ_XDRFILE      : If set CPPTRAJ has XDR file support.
#   CPPTRAJ_TNGFILE      : If set CPPTRAJ has TNG file support.
#   CPPTRAJ_SINGLE_ENS   : If set CPPTRAJ has single ensemble support.
#   CPPTRAJ_OPENMM       : If set CPPTRAJ has OpenMM support.
# ----- Variables that can be set by individual scripts ----
#   TOP                  : Topology file for cpptraj
#   INPUT                : Input file for cpptraj
#   CPPTRAJ_TEST_MODE    : Set to 'master' if executed from CpptrajTest.sh.
# ----- Local setup variables ------------------------------
VGMODE=0                 # Valgrind mode: 0 none, 1 memcheck, 2 helgrind
USE_DACDIF=1             # If 0 do not use dacdif even if in AmberTools
GET_TIMING=0             # If 1 time cpptraj with the CPPTRAJ_TIME binary
SUMMARY=0                # If 1 print summary of CPPTRAJ_TEST_RESULTS only
SHOWERRORS=0             # If 1, print test errors to STDOUT after summary.
SFX=""                   # CPPTRAJ binary suffix
TARGET=""                # Make target if multiple tests being run
STANDALONE=1             # If 0, part of AmberTools. If 1, stand-alone (e.g. from GitHub).
TEST_DIRS=''             # Specific test directories to run.
TEST_SKIPPED=0           # If 1 this test has been skipped via SkipTest.
EXIT_ON_ERROR=0          # (Debug) If 1 exit when a specified test fails.
# ----- Variables local to single test ---------------------
TEST_WORKDIR=''          # Test working directory
NUMCOMPARISONS=0         # Total number of times DoTest has been called this test.
ERRCOUNT=0               # Total number of errors detected by DoTest this test.
SKIPCOUNT=0              # Total number of times SkipCheck called by this test.
WARNCOUNT=0              # Total number of warnings detected by DoTest this test.
PROGCOUNT=0              # Total number of times RunCpptraj has been called this test.
PROGERROR=0              # Total number of program errors this test
CHECKERR=0               # Total errors this test from CheckEnv routine.
TESTNAME=''              # Current test name for Requires routine.
UNITNAME=''              # Current unit name for CheckFor routine.
DESCRIP=''               # Current test/unit name for CheckEnv routine.
FNC1=''                  # First file to NcTest(); for output in DoTest()
FNC2=''                  # Second file to NcTest(); for output in DoTest()

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

# ErrMsg() <message>
# Write out error message to stderr prefaced with 'Error:'.
ErrMsg() {
  >&2 echo "Error: $*"
}

# OUT() <message>
#   Write message to CPPTRAJ_TEST_RESULTS
OUT() {
  echo "$1" >> $CPPTRAJ_TEST_RESULTS
}

# OutBoth() <message>
#   Send <message> to CPPTRAJ_TEST_RESULTS and CPPTRAJ_TEST_ERROR
OutBoth() {
  OUT "$1"
  # Write header to CPPTRAJ_TEST_ERROR if this is the first error.
  if [ $ERRCOUNT -eq 0 ] ; then
    TestHeader "$CPPTRAJ_TEST_ERROR"
  fi
  echo "$1" >> $CPPTRAJ_TEST_ERROR
}

# ------------------------------------------------------------------------------
# CheckTestFiles() <test save> <test output>
#   Check that <test save> and <test output> exist.
#   \return 0 if both present, 1 if either one absent.
CheckTestFiles() {
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    # Standalone output.
    if [ ! -f "$1" ] ; then
      OutBoth "  Save file '$1' not found."
    elif [ ! -f "$2" ] ; then
      OutBoth "  Test output '$2' not found."
    else
      return 0
    fi
  else
    # AmberTools output
    if [ ! -f "$1" ] ; then
      echo "possible FAILURE:  file $1 does not exist."
    elif [ ! -f "$2" ] ; then
      echo "possible FAILURE:  file $2 does not exist."
    else
      return 0
    fi
  fi
  return 1
}

# ------------------------------------------------------------------------------
# DoTest() <File1> <File2> [-r <relative err>] [-a <absolute err>] [<arg1>] ... [<argN>]
#   Compare File1 (the 'save' file) to File2 (test output), print an error if
#   they differ. If '-r' or '-a' are specified the test will pass as long as the
#   maximum relative or absolulte errors are below the given value (using Nelson
#   H. F. Beebe's ndiff.awk script. The remaining args can be used to pass
#   options to CPPTRAJ_DIFF.
DoTest() {
  # Use diff, or ndiff where '-r' or '-a' specified.
  ((NUMCOMPARISONS++))
  DIFFARGS='--strip-trailing-cr'
  NDIFFARGS=""
  # First two arguments are files to compare.
  F1=$1 ; shift
  F2=$1 ; shift
  # Names for output
  if [ "$F1" = 'nc0.save' -a "$F2" = 'nc0' ] ; then
    FNAME1=$FNC1
    FNAME2=$FNC2
  else
    FNAME1=$F1
    FNAME2=$F2
  fi
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
  CheckTestFiles $F1 $F2
  if [ $? -ne 0 ] ; then
    ((ERRCOUNT++))
  else
    if [ ! -z "$CPPTRAJ_DACDIF" ] ; then
      # Print AT test header.
      echo "diffing $FNAME1 with $FNAME2"
    fi
    if [ $USE_NDIFF -eq 0 ] ; then
      $CPPTRAJ_DIFF $DIFFARGS $DIFFOPTS $F1 $F2 > temp.diff 2>&1
    else
      awk -f $CPPTRAJ_NDIFF $NDIFFARGS $F1 $F2 > temp.diff 2>&1
    fi
    if [ -s 'temp.diff' ] ; then
      if [ -z "$CPPTRAJ_DACDIF" ] ; then
        OutBoth "  $FNAME1 $FNAME2 are different."
        cat temp.diff >> $CPPTRAJ_TEST_ERROR
      else
        mv temp.diff $FNAME2.dif
        echo "possible FAILURE:  check $FNAME2.dif"
        echo "possible FAILURE:  check $FNAME2.dif" >> $CPPTRAJ_TEST_ERROR
        echo "$TEST_WORKDIR" >> $CPPTRAJ_TEST_ERROR
        cat $FNAME2.dif >> $CPPTRAJ_TEST_ERROR
        echo "---------------------------------------" >> $CPPTRAJ_TEST_ERROR
      fi
      ((ERRCOUNT++))
    else
      if [ -z "$CPPTRAJ_DACDIF" ] ; then
        # Standalone pass.
        OUT  "  $FNAME2 OK."
      else
        # AmberTools pass.
        echo "PASSED"
      fi
    fi
    $CPPTRAJ_RM temp.diff
  fi
  if [ ! -z "$CPPTRAJ_DACDIF" ] ; then
    # Print AT test footer.
    echo "=============================================================="
  fi
}

# ------------------------------------------------------------------------------
# NcTest() <1> <2>
#   Compare NetCDF files <1> and <2>. Use CPPTRAJ_NCDUMP to convert to ASCII
#   first, removing ==> line and :programVersion attribute.
NcTest() {
  #echo "DEBUG: NcTest $1 $2"
  if [ -z "$1" -o -z "$2" ] ; then
    ErrMsg "NcTest(): One or both files not specified."
    exit 1
  fi
  # Save remaining args for DoTest
  FNC1=$1
  FNC2=$2
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
  CheckTestFiles $FNC1 $FNC2
  if [ $? -ne 0 ] ; then
    ((NUMCOMPARISONS++))
    ((ERRCOUNT++))
  else
    # Prepare files.
    if [ $CALC_NUM_ERR -eq 1 ] ; then
      # FIXME: Must remove commas here because I cannot figure out how to pass
      # the regular expression to ndiff.awk FS without the interpreter giving
      # this error for FS='[ \t,()]':
      # awk: fatal: Invalid regular expression: /'[/
      $CPPTRAJ_NCDUMP -n nctest $FNC1 | grep -v "==>\|:program" | sed 's/,/ /g' > nc0.save
      $CPPTRAJ_NCDUMP -n nctest $FNC2 | grep -v "==>\|:program" | sed 's/,/ /g' > nc0
    else
      $CPPTRAJ_NCDUMP -n nctest $FNC1 | grep -v "==>\|:program" > nc0.save
      $CPPTRAJ_NCDUMP -n nctest $FNC2 | grep -v "==>\|:program" > nc0
    fi
    DoTest $DIFFARGS
    $CPPTRAJ_RM nc0.save nc0
  fi
}

# ------------------------------------------------------------------------------
# GetResultsFiles() <mode> <base>
#   Consolidate files with name <base>. If <mode> is 'single' assume one file in
#   the current directory, otherwise assume files are in subdirectories.
GetResultsFiles() {
  if [ ! -z "$TEST_DIRS" ] ; then
    RESULTSFILES=''
    for DIR in $TEST_DIRS ; do
      if [ -f "$DIR/$2" ] ; then
        RESULTSFILES="$RESULTSFILES $DIR/$2"
      fi
    done
  else
    if [ "$1" = 'single' ] ; then
      RESULTSFILES=`ls $2 2> tmp.cpptrajtest.devnull`
    else
      RESULTSFILES=`ls */$2 2> tmp.cpptrajtest.devnull`
    fi
    rm tmp.cpptrajtest.devnull
  fi
}

# ------------------------------------------------------------------------------
ParseValgrindOut() {
  awk 'BEGIN{
    ntests = 0;
    #n_vg_err = 0;
    #n_vg_allocs = 0;
    #n_vg_frees = 0;
    #n_vg_allocated = 0;
    #n_vg_leaks = 0;
    currentTest = "";
    byte_to_mb = 1.0 / (1024.0 * 1024.0);
  }
  function WriteOut(outfile)
  {
    printf("Valgrind summary for %i program executions:\n", ntests) >> outfile;
    for (i = 1; i <= ntests; i++) {
      printf("  %s %i errors,", vg_id[i], n_vg_err[i]) >> outfile;
      printf(" %i allocs, %i frees,",  n_vg_allocs[i], n_vg_frees[i]) >> outfile;
      if (n_vg_allocated[i] > 1000) {
        size_in_mb = n_vg_allocated[i] * byte_to_mb;
        printf(" %.3f MB allocated,", size_in_mb) >> outfile;
      } else
        printf(" %i bytes allocated,", n_vg_allocated[i]) >> outfile;
      printf(" %i memory leaks.\n", n_vg_leaks[i]) >> outfile;
    }
  }
  {
    if (index($1,"==") != 0) {
      if ($1 != currentTest) {
        ntests++;
        currentTest = $1;
        vg_id[ntests] = currentTest;
      }
      if ($2 == "total" && $3 == "heap") {
        gsub(/,/,"");
        n_vg_allocs[ntests] += $5;
        n_vg_frees[ntests] += $7;
        n_vg_allocated[ntests] += $9;
      } else if ($2 == "ERROR" && $3 == "SUMMARY:") {
        gsub(/,/,"");
        n_vg_err[ntests] += $4;
      } else if ($2 == "LEAK" && $3 == "SUMMARY:") {
        gsub(/,/,"");
        n_vg_leaks[ntests]++;
      }
    }
  }END{
    WriteOut("vg.summary");
  }' $1
  cat vg.summary >> $CPPTRAJ_TEST_RESULTS
  cat vg.summary
  rm vg.summary
}

# ------------------------------------------------------------------------------
# Summary()
#  Print a summary of results in all CPPTRAJ_TEST_RESULTS files and exit.
#  Optionally print an error summary and/or valgrind summary.
Summary() {
  MODE=''
  if [ "$CPPTRAJ_TEST_MODE" = 'master' ] ; then
    MODE='multi'
  else
    MODE='single'
  fi
  ERR_STATUS=0
  PRINT_FOOTER=0
  if [ ! -z "$CPPTRAJ_TEST_RESULTS" ] ; then
    GetResultsFiles $MODE $CPPTRAJ_TEST_RESULTS
    #echo "DEBUG: Getting results from $RESULTSFILES"
    if [ ! -z "$RESULTSFILES" ] ; then
      if [ "$MODE" = 'multi' ] ; then
        cat $RESULTSFILES > $CPPTRAJ_TEST_RESULTS
      fi
      PRINT_FOOTER=1
      echo "===================== TEST SUMMARY ======================"
      awk 'BEGIN{
        comparisons_ok = 0;   # Number of OK comparisons
        comparisons_warn = 0; # Number of comparisons with a warning
        comparisons_diff = 0; # Number of comparisons with diff/error
        comparisons_skip = 0; # Number of comparisons skipped.
        program_err = 0;      # Number of program errors.
        program_exe = 0;      # Number of program executions.
        ntests = 0;           # Number of tests.
        skip_tests = 0;       # Number of tests completely skipped.
        ok_tests = 0;         # Number of completely OK tests.
      }{
        if ( $1 == "TEST:" )
          ntests++;
        else if ($1 == "Warning:")
          comparisons_warn++;
        else if ($1 == "CPPTRAJ:")
          program_exe++;
        else if ($1 == "SKIP:")
          skip_tests++;
        else if ($1 == "Error:")
          program_err++;
        else if ($NF == "OK.")
          comparisons_ok++;
        else if ($NF == "different.")
          comparisons_diff++;
        else if ($(NF-1) == "not" && $NF == "found.")
          comparisons_diff++;
        else if ($(NF-1) == "comparisons") {
          if ($NF == "passed.")
            ok_tests++;
          else if ($NF == "skipped.")
            comparisons_skip += $1;
        }
      }END{
        n_comps = comparisons_ok + comparisons_warn + comparisons_diff + comparisons_skip;
        if (n_comps > 0)
          printf("  %i out of %i comparisons OK (%i failed, %i warnings, %i skipped).\n",
                 comparisons_ok, n_comps, comparisons_diff, comparisons_warn, comparisons_skip);
        if (program_exe > 0)
          printf("  %i out of %i program executions completed.\n",
                 program_exe - program_err, program_exe);
        printf("  %i out of %i tests completed with no issues.\n",
               ok_tests, ntests);
        if (skip_tests > 0)
          printf("  %i tests skipped.\n", skip_tests);
        exit (comparisons_diff + program_err);
      }' $CPPTRAJ_TEST_RESULTS
      ERR_STATUS=$?
      if [ $SHOWERRORS -eq 1 ] ; then
        grep "SKIP:" $CPPTRAJ_TEST_RESULTS > tmpout
        grep "Skipped test:" $CPPTRAJ_TEST_RESULTS >> tmpout
        if [ -s tmpout ] ; then
          echo ""
          echo "Skipped tests/comparisons:"
          echo "---------------------------------------------------------"
          cat tmpout
          echo "---------------------------------------------------------"
        fi
        $CPPTRAJ_RM tmpout
      fi
    fi
  fi
  # Error summary
  if [ ! -z "$CPPTRAJ_TEST_ERROR" ] ; then
    GetResultsFiles $MODE $CPPTRAJ_TEST_ERROR
    #echo "DEBUG: Getting errors from $RESULTSFILES"
    if [ ! -z "$RESULTSFILES" ] ; then
      if [ "$MODE" = 'multi' ] ; then
        cat $RESULTSFILES > $CPPTRAJ_TEST_ERROR
      fi
      if [ $SHOWERRORS -eq 1 ] ; then
        echo ""
        echo "Obtained the following errors:"
        echo "---------------------------------------------------------"
        cat $CPPTRAJ_TEST_ERROR
        echo "---------------------------------------------------------"
      fi
    fi
  fi
  # Valgrind summary
  if [ ! -z "$VALGRIND" ] ; then
    GetResultsFiles $MODE $CPPTRAJ_ERROR
    if [ ! -z "$RESULTSFILES" ] ; then
      echo "---------------------------------------------------------"
      echo "Valgrind summary:"
      cat $RESULTSFILES | awk 'BEGIN{
        n_vg_err = 0;
        n_vg_leaks = 0;
      }{
        if ($2 == "ERROR" && $3 == "SUMMARY:") {
          gsub(/,/,"");
          n_vg_err += $4;
        } else if ($2 == "LEAK" && $3 == "SUMMARY:") {
          n_vg_leaks++;
        }
      }END{
        printf("    %i errors,", n_vg_err);
        printf(" %i memory leak reports.\n", n_vg_leaks);
      }'
    fi
  fi
  if [ $PRINT_FOOTER -eq 1 ] ; then
    echo "========================================================="
  fi
  exit $ERR_STATUS
}

# ------------------------------------------------------------------------------
# ProgramError() <message>
ProgramError() {
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    ErrMsg " $1"
    OutBoth "Error: $1"
  else
    if [ -z "$2" ] ; then
      PNAME=$CPPTRAJ
    else
      PNAME=$2
    fi
    echo "$PNAME: Program error"
  fi
  ((PROGERROR++))
}

# ------------------------------------------------------------------------------
# RunCpptraj() <title>
#   Run cpptraj test with given title.
RunCpptraj() {
  # If only cleaning requested no run needed, exit now
  if [ $CPPTRAJ_TEST_CLEAN -eq 1 ] ; then
    exit 0
  fi
  ((PROGCOUNT++))
  echo ""
  echo "  CPPTRAJ: $1"
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    OUT "  CPPTRAJ: $1"
  fi
  cpptraj_cmd="$CPPTRAJ_TIME $DO_PARALLEL $VALGRIND $CPPTRAJ $TOP $INPUT $CPPTRAJ_RNG $CPPTRAJ_DEBUG"
  if [ ! -z "$CPPTRAJ_DEBUG" ] ; then
    echo "$cpptraj_cmd >> $CPPTRAJ_OUTPUT 2>>$CPPTRAJ_ERROR"
  fi
  # There are 4 different ways to redirect output depending on if
  # CPPTRAJ_OUTPUT/CPPTRAJ_ERROR are defined or not.
  if [ ! -z "$CPPTRAJ_OUTPUT" ] ; then
    # CPPTRAJ_OUTPUT is defined.
    if [ -z "$CPPTRAJ_ERROR" ] ; then
      $cpptraj_cmd >> $CPPTRAJ_OUTPUT
    else
      $cpptraj_cmd >> $CPPTRAJ_OUTPUT 2>> $CPPTRAJ_ERROR
    fi
  else
    # CPPTRAJ_OUTPUT is not defined.
    if [ -z "$CPPTRAJ_ERROR" ] ; then
      $cpptraj_cmd
    else
      $cpptraj_cmd 2>> $CPPTRAJ_ERROR
    fi
  fi
  STATUS=$?
  #echo "DEBUG: Cpptraj exited with status $STATUS"
  if [ $STATUS -ne 0 ] ; then
    ProgramError "cpptraj exited with status $STATUS"
  fi
}

# ------------------------------------------------------------------------------
# EndTest()
#   Print a summary of the current tests results. Should be called at the end of
#   every test script.
EndTest() {
  #echo "DEBUG: EndTest"
  # Report only when not using dacdif and test not skipped.
  if [ -z "$CPPTRAJ_DACDIF" -a $TEST_SKIPPED -eq 0 ] ; then
    echo ""
    if [ $PROGERROR -gt 0 ] ; then
      echo "  $PROGERROR out of $PROGCOUNT executions exited with an error."
    fi
    if [ $ERRCOUNT -gt 0 ] ; then
      echo    "  $ERRCOUNT out of $NUMCOMPARISONS comparisons failed."
      OutBoth "  $ERRCOUNT out of $NUMCOMPARISONS comparisons failed."
    elif [ $WARNCOUNT -gt 0 ] ; then
      ((PASSCOUNT = $NUMCOMPARISONS - $WARNCOUNT))
      echo "  $PASSCOUNT out of $NUMCOMPARISONS passing comparisons. $WARNCOUNT warnings."
      OUT  "  $PASSCOUNT out of $NUMCOMPARISONS passing comparisons. $WARNCOUNT warnings."
    elif [ $NUMCOMPARISONS -gt 0 ] ; then
      echo "All $NUMCOMPARISONS comparisons passed."
      OUT  "All $NUMCOMPARISONS comparisons passed."
    fi
    if [ $SKIPCOUNT -gt 0 ] ; then
      echo "$SKIPCOUNT comparisons skipped."
      OUT  "$SKIPCOUNT comparisons skipped."
    fi
    echo ""
    if [ ! -z "$VALGRIND" ] ; then
      ParseValgrindOut $CPPTRAJ_ERROR
      echo ""
    fi
    if [ $CPPTRAJ_PROFILE -eq 1 ] ; then
      if [ -e 'gmon.out' ] ; then
        gprof $CPPTRAJ > profiledata.txt
      fi
    fi
  fi
  if [ $SUMMARY -ne 0 ] ; then
    # NOTE: SUMMARY should only be set here if single test.
    Summary
  fi
  if [ $EXIT_ON_ERROR -eq 1 ] ; then
    if [ $PROGERROR -ne 0 -o $ERRCOUNT -ne 0 ] ; then
      exit 1
    elif [ $NUMCOMPARISONS -eq 0 -a $TEST_SKIPPED -eq 0 ] ; then
      ErrMsg "Zero comparisons and test not skipped."
      exit 1
    fi
  fi
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

#-------------------------------------------------------------------------------
# SetNthreads()
#   Use CPPTRAJ_NPROC to set N_THREADS if not already set.
SetNthreads() {
  if [ -z "$N_THREADS" ] ; then
    if [ ! -f "$CPPTRAJ_NPROC" ] ; then
      export N_THREADS=1
      return 1
    fi
    export N_THREADS=`$CPPTRAJ_NPROC`
    echo "  $N_THREADS MPI threads."
  fi
  return 0
}

#-------------------------------------------------------------------------------
# Help(): Print help
Help() {
  echo "Command line flags"
  if [ "$CPPTRAJ_TEST_MODE" = 'master' ] ; then
    echo "  <test dir 1> [<test dir 2> ...]"
  fi
  echo "  mpi             : Use MPI version of CPPTRAJ (automatically used if DO_PARALLEL set)."
  echo "  openmp          : Use OpenMP version of CPPTRAJ."
  echo "  cuda            : Use CUDA version of CPPTRAJ."
  echo "  hip             : Use HIP version of CPPTRAJ."
  echo "  time            : Time the test with 'time'."
  echo "  stdout          : Print CPPTRAJ test output to STDOUT."
  echo "  vg              : Run test with 'valgrind' memcheck."
  echo "  vgh             : Run test with 'valgrind' helgrind."
  echo "  clean           : Do not run test; clean test output."
  echo "  summary         : Print summary of tests run."
  echo "  showerrors      : (summary only) Print all test errors to STDOUT after summary."
  echo "  -nodacdif       : Do not use dacdif for test comparisons (AmberTools only)."
  echo "  -d              : Run CPPTRAJ with global debug level 4."
  echo "  -debug <#>      : Run CPPTRAJ with global debug level #."
  echo "  -cpptraj <file> : Use CPPTRAJ binary <file>."
  echo "  -profile        : Profile results with 'gprof' (requires special compile)."
  echo "Important environment variables"
  echo "  DO_PARALLEL     : MPI run command."
  echo "  N_THREADS       : Number of MPI threads (only needed if 'DO_PARALLEL nproc' fails)."
  echo "  OMP_NUM_THREADS : Number of OpenMP threads to use."
  echo "  CPPTRAJ_DEBUG   : Can be used to pass extra flags to CPPTRAJ."
  echo ""
}

#-------------------------------------------------------------------------------
# CmdLineOpts()
#   Process test script command line options. Only executed if
#   CPPTRAJ_TEST_SETUP is not already set.
CmdLineOpts() {
  CPPTRAJ_TEST_CLEAN=0 # Will be exported
  CPPTRAJ_LONG_TEST=0  # Will be exported
  CPPTRAJ_PROFILE=0    # Will be exported
  SFX_OMP=0
  SFX_CUDA=0
  SFX_HIP=0
  SFX_MPI=0
  GET_TIMING=0
  while [ ! -z "$1" ] ; do
    case "$1" in
      "clean"     ) CPPTRAJ_TEST_CLEAN=1 ;;
      "long"      ) CPPTRAJ_LONG_TEST=1 ;;
      "summary"   ) SUMMARY=1 ;;
      "showerrors") SHOWERRORS=1 ;;
      "stdout"    ) CPPTRAJ_OUTPUT='' ;;
      "openmp"    ) SFX_OMP=1 ;;
      "cuda"      ) SFX_CUDA=1 ;;
      "hip"       ) SFX_HIP=1 ;;
      "mpi"       ) SFX_MPI=1 ;;
      "vg"        ) VGMODE=1 ;;
      "vgh"       ) VGMODE=2 ;;
      "time"      ) GET_TIMING=1 ;;
      "timing"    ) GET_TIMING=2 ;;
      "-d"        ) CPPTRAJ_DEBUG="$CPPTRAJ_DEBUG -debug 4" ;;
      "-exitonerr") EXIT_ON_ERROR=1 ;;
      "-debug"    ) shift ; CPPTRAJ_DEBUG="$CPPTRAJ_DEBUG -debug $1" ;;
      "-nodacdif" ) USE_DACDIF=0 ;;
      "-cpptraj"  ) shift ; export CPPTRAJ=$1 ; echo "Using cpptraj: $CPPTRAJ" ;;
      "--target"  ) shift ; TARGET=$1 ;;
      "-profile"  ) CPPTRAJ_PROFILE=1 ; echo "Performing gnu profiling during EndTest." ;;
      "-h" | "--help" ) Help ; exit 0 ;;
      "alltests"  )
        echo "Running all tests in Test_* directories."
        for DIR in `ls -d Test_*` ; do
          if [ -f "$DIR/RunTest.sh" ] ; then
            TEST_DIRS="$TEST_DIRS $DIR"
          fi
        done
        ;;
      *           )
        if [ -d "$1" -a -f "$1/RunTest.sh" ] ; then
          # Assume this is a test we want to run.
          TEST_DIRS="$TEST_DIRS $1"
        else
          ErrMsg "Unknown option: $1"
          exit 1
        fi
        ;;
    esac
    shift
  done
  if [ $CPPTRAJ_PROFILE -eq 1 ] ; then
    Required "gprof"
  fi
  export CPPTRAJ_TEST_CLEAN
  export CPPTRAJ_LONG_TEST
  export CPPTRAJ_DEBUG
  export CPPTRAJ_PROFILE
  # If DO_PARALLEL has been set force MPI
  if [ ! -z "$DO_PARALLEL" ] ; then
    SFX_MPI=1
  elif [ $SFX_MPI -eq 1 ] ; then
    ErrMsg "'mpi' specified but DO_PARALLEL not set."
    exit 1
  fi
  # Warn if using OpenMP but OMP_NUM_THREADS not set.
  if [ $SFX_OMP -eq 1 -a -z "$OMP_NUM_THREADS" ] ; then
    echo "Warning: Using OpenMP but OMP_NUM_THREADS is not set."
  fi
  # Set up SFX
  if [ "$SFX_MPI"  -eq 1 ] ; then SFX="$SFX.MPI"  ; fi
  if [ "$SFX_OMP"  -eq 1 ] ; then SFX="$SFX.OMP"  ; fi
  if [ "$SFX_CUDA" -eq 1 ] ; then SFX="$SFX.cuda" ; fi
  if [ "$SFX_HIP"  -eq 1 ] ; then SFX="$SFX.hip"  ; fi
}

#-------------------------------------------------------------------------------
# Required() <binary>
#   Insure that specified binary is in the PATH
Required() {
  if [ -z "`which $1`" ] ; then
    ErrMsg "Required binary '$1' not found."
    exit 1
  fi
}

#-------------------------------------------------------------------------------
# CheckDefines()
#   Check how CPPTRAJ was compiled.
CheckDefines() {
  CPPTRAJ_XDRFILE='yes'
  CPPTRAJ_TNGFILE=''
  CPPTRAJ_MATHLIB='yes'
  CPPTRAJDEFINES=`$CPPTRAJ --defines`
  if [ $? -ne 0 ] ; then
    ErrMsg "Could not execute '$CPPTRAJ --defines'"
    exit 1
  fi
  for DEFINE in `$CPPTRAJ --defines` ; do
    case "$DEFINE" in
      '-DHASGZ'         ) export CPPTRAJ_ZLIB=$DEFINE ;;
      '-DHASBZ2'        ) export CPPTRAJ_BZLIB=$DEFINE ;;
      '-DBINTRAJ'       ) export CPPTRAJ_NETCDFLIB=$DEFINE ;;
      '-DC11_SUPPORT'   ) export CPPTRAJ_C11_SUPPORT=$DEFINE ;;
      '-DLIBPME'        ) export CPPTRAJ_LIBPME=$DEFINE ;;
      '-DMPI'           ) export CPPTRAJ_MPILIB=$DEFINE ;;
      '-DNO_MATHLIB'    ) CPPTRAJ_MATHLIB='' ;;
      '-D_OPENMP'       ) export CPPTRAJ_OPENMP=$DEFINE ;;
      '-DHAS_PNETCDF'   ) export CPPTRAJ_PNETCDFLIB=$DEFINE ;;
      '-DUSE_SANDERLIB' ) export CPPTRAJ_SANDERLIB=$DEFINE ;;
      '-DFFTW_FFT'      ) export CPPTRAJ_FFTW_FFT=$DEFINE ;;
      '-DCUDA'          ) export CPPTRAJ_CUDA=$DEFINE ;;
      '-DNO_XDRFILE'    ) CPPTRAJ_XDRFILE='' ;;
      '-DHAS_TNGFILE'   ) CPPTRAJ_TNGFILE='$DEFINE' ;;
      '-DENABLE_SINGLE_ENSEMBLE' ) export CPPTRAJ_SINGLE_ENS=$DEFINE ;;
      '-DHAS_OPENMM'    ) export CPPTRAJ_OPENMM=$DEFINE ;;
    esac
  done
  export CPPTRAJ_XDRFILE
  export CPPTRAJ_TNGFILE
  export CPPTRAJ_MATHLIB
  #echo "DEBUG: $ZLIB $BZLIB $NETCDFLIB $MPILIB $NOMATHLIB $OPENMP $PNETCDFLIB $SANDERLIB $CUDA $HIP $XDRFILE"
}

#-------------------------------------------------------------------------------
# SetBinaries()
#   Set paths for all binaries if not already set.
SetBinaries() {
  # Guess where things might be depending on if we are in AT or not, etc.
  PATH_TYPE='absolute'
  if [ $STANDALONE -eq 0 ] ; then
    DIRPREFIX=$AMBERHOME
  elif [ -z "$CPPTRAJHOME" ] ; then
    # Need to do relative instead of absolute path.
    PATH_TYPE='relative'
    if [ "$CPPTRAJ_TEST_MODE" = 'master' ] ; then
      DIRPREFIX=..
    else
      DIRPREFIX=../..
    fi
  else
    DIRPREFIX=$CPPTRAJHOME
  fi
  # Determine location of CPPTRAJ if not specified.
  if [ -z "$CPPTRAJ" ] ; then
    CPPTRAJ=$DIRPREFIX/bin/cpptraj$SFX
    if [ ! -f "$CPPTRAJ" ] ; then
      ErrMsg "CPPTRAJ $CPPTRAJ not found."
      exit 1
    fi
    export CPPTRAJ
  fi
  # Check how cpptraj was compiled.
  CheckDefines
  # If compiled with NetCDF support ensure ncdump is available
  if [ ! -z "$CPPTRAJ_NETCDFLIB" ] ; then
    Required "ncdump"
    export CPPTRAJ_NCDUMP=`which ncdump`
  fi
  # Set up valgrind if necessary.
  if [ $VGMODE -ne 0 ] ; then
    echo "  Using valgrind."
    CPPTRAJ_ERROR='valgrind.out'
    if [ ! -z "$CPPTRAJ_CUDA" ] ; then
      if [ $VGMODE -eq 1 ] ; then
        Required "cuda-memcheck"
        VALGRIND='cuda-memcheck'
      else
        ErrMsg "'helgrind' not supported for CUDA."
        exit 1
      fi
    elif [ ! -z "$CPPTRAJ_HIP" ] ; then
      if [ $VGMODE -eq 1 ] ; then
        Required "hip-memcheck"
        VALGRIND='hip-memcheck'
      else
        ErrMsg "'helgrind' not supported for HIP."
        exit 1
      fi
    else
      Required "valgrind"
      if [ $VGMODE -eq 1 ] ; then
        VALGRIND="valgrind --tool=memcheck --leak-check=yes --show-reachable=yes"
      elif [ $VGMODE -eq 2 ] ; then
        VALGRIND="valgrind --tool=helgrind"
      else
        ErrMsg "Unsupported VGMODE $VGMODE"
        exit 1
      fi
    fi
    export VALGRIND
  fi
  # Determine location of ndiff.awk
  if [ -z "$CPPTRAJ_NDIFF" ] ; then
    CPPTRAJ_NDIFF=$CPPTRAJ_TEST_ROOT/utilities/ndiff/ndiff.awk
    if [ ! -f "$CPPTRAJ_NDIFF" ] ; then
      ErrMsg "'ndiff.awk' not present: $CPPTRAJ_NDIFF"
      exit 1
    fi
    export CPPTRAJ_NDIFF
  fi
  # Determine location of nproc/numprocs
  if [ -z "$CPPTRAJ_NPROC" ] ; then
    if [ ! -z "$DO_PARALLEL" ] ; then
      CPPTRAJ_NPROC=$CPPTRAJ_TEST_ROOT/utilities/nproc.sh
      if [ ! -f "$CPPTRAJ_NPROC" ] ; then
        ErrMsg "Error: nproc $CPPTRAJ_NPROC not found."
        exit 1
      fi
      export CPPTRAJ_NPROC
      SetNthreads
    fi
  fi
  # Determine timing if necessary
  if [ $GET_TIMING -gt 0 ] ; then
    Required "time"
    export CPPTRAJ_TIME=`which time`
  fi
  # Report binary details
  if [ $STANDALONE -eq 1 ] ; then
    ls -l $CPPTRAJ
  fi
  # Print DEBUG info
  if [ ! -z "$CPPTRAJ_DEBUG" ] ; then
    if [ $STANDALONE -eq 1 ] ; then
      echo "DEBUG: Standalone mode."
    else
      echo "DEBUG: AmberTools mode."
    fi
    echo "DEBUG: CPPTRAJ: $CPPTRAJ"
    echo "DEBUG: NPROC:   $CPPTRAJ_NPROC"
    echo "DEBUG: NCDUMP:  $CPPTRAJ_NCDUMP"
    echo "DEBUG: DIFFCMD: $CPPTRAJ_DIFF"
    echo "DEBUG: DACDIF:  $CPPTRAJ_DACDIF"
    echo "DEBUG: NDIFF:   $CPPTRAJ_NDIFF"
  fi
  # If path type is relative and we are not in an indivdual test
  # directories the path needs to be incremented one dir up.
  if [ "$PATH_TYPE" = 'relative' -a "$CPPTRAJ_TEST_MODE" = 'master' ] ; then
    CPPTRAJ="../$CPPTRAJ"
  #  CPPTRAJ_NDIFF="../$CPPTRAJ_NDIFF"
  #  CPPTRAJ_NPROC="../$CPPTRAJ_NPROC"
  fi
}

# ------------------------------------------------------------------------------
# Library/environment Check functions
# ------------------------------------------------------------------------------

# SkipTest() <description>
#  Skip an entire test directory.
SkipTest() {
  echo "  SKIP: $1"
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    OUT "  SKIP: $1"
  fi
  echo ""
  TEST_SKIPPED=1
  EndTest
  exit 0
}

# SkipCheck() <description>
#  Skip an individual test.
SkipCheck() {
  echo "  Skipped test: $1"
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    OUT "  Skipped test: $1"
  fi
  ((SKIPCOUNT++))
  # Reset check count FIXME needed?
  CHECKERR=0
}

# TestLibrary() <name> <env. var>
#   Test for presence of cpptraj capability by given environment var, set by
#   CheckDefines.
TestLibrary() {
  if [ -z "$2" ] ; then
    echo "  $DESCRIP requires $1."
    echo "  Cpptraj was compiled without $1 support."
    ((CHECKERR++))
  fi
}

# CheckEnv() <list>
# Check for the given list of requirements. For each requirement not meant
# increment CHECKERR by 1.
# netcdf         : NetCDF support
# zlib           : Zlib (gzip) support
# bzlib          : Bzlib support
# xdr            : XDR file support
# mathlib        : BLAS/LAPACK/ARPACK support
# sanderlib      : SANDER API support
# pnetcdf        : Parallel NetCDF support
# notparallel    : Test should not be run in parallel
# parallel       : Test should only be run in parallel
# maxthreads <#> : Test should not be run with more than <#> MPI threads
# nthreads <#>   : Test requires multiples of <#> MPI threads in parallel
# threads <#>    : Test requires exactly <#> threads in parallel.
# amberhome      : Test requires AMBERHOME set
# amberorcpptraj : Test requires AMBERHOME or CPPTRAJHOME set
# inpath <name>  : Test requires <name> to be in PATH
# testos <os>    : Test requires specific OS
# file <file>    : Test requires specified file
# openmp         : Test requires openmp
# singleensemble : Single NetCDF ensemble file support.
# long           : Test requires the 'long' option be specified.
# disabled       : Test is disabled, just count as skipped.
CheckEnv() {
  #echo "DEBUG: CheckEnv() $*"
  if [ -z "$DESCRIP" ] ; then
    ErrMsg "CheckEnv() called with TESTNAME/UNITNAME unset."
    exit 1
  fi
  CHECKERR=0
  while [ ! -z "$1" ] ; do
    #echo "DEBUG: $DESCRIP: Checking requirement: $1"
    case "$1" in
      'netcdf'    ) TestLibrary "NetCDF"             "$CPPTRAJ_NETCDFLIB" ;;
      'c++11'     ) TestLibrary "C++11"              "$CPPTRAJ_C11_SUPPORT" ;;
      'libpme'    ) TestLibrary "libPME"             "$CPPTRAJ_LIBPME" ;;
      'zlib'      ) TestLibrary "Zlib"               "$CPPTRAJ_ZLIB" ;;
      'bzlib'     ) TestLibrary "Bzlib"              "$CPPTRAJ_BZLIB" ;;
      'xdr'       ) TestLibrary "XDR file"           "$CPPTRAJ_XDRFILE" ;;
      'tng'       ) TestLibrary "TNG file"           "$CPPTRAJ_TNGFILE" ;;
      'mathlib'   ) TestLibrary "BLAS/LAPACK/ARPACK" "$CPPTRAJ_MATHLIB" ;;
      'sanderlib' ) TestLibrary "SANDER API from AmberTools" "$CPPTRAJ_SANDERLIB" ;;
      'fftw'      ) TestLibrary "FFTW"               "$CPPTRAJ_FFTW_FFT" ;;
      'openmp'    ) TestLibrary "OpenMP"             "$CPPTRAJ_OPENMP" ;;
      'singleensemble' ) TestLibrary "Single ensemble support" "$CPPTRAJ_SINGLE_ENS" ;;
      'cuda'      ) TestLibrary "CUDA"               "$CPPTRAJ_CUDA" ;;
      'hip'       ) TestLibrary "HIP"                "$CPPTRAJ_HIP" ;;
      'openmm'    ) TestLibrary "OpenMM"             "$CPPTRAJ_OPENMM" ;;
      'notcuda'   )
	if [ ! -z "$CPPTRAJ_CUDA" ]; then
          echo "  $DESCRIP cannot be run on CUDA."
	  ((CHECKERR++))
        fi
	;;
      'pnetcdf'   )
        if [ ! -z "$DO_PARALLEL" ] ; then
          TestLibrary "Parallel NetCDF" "$CPPTRAJ_PNETCDFLIB"
        fi
        ;;
      'notparallel' )
        if [ ! -z "$DO_PARALLEL" ] ; then
          echo "  $DESCRIP cannot be run in parallel."
          ((CHECKERR++))
        fi
        ;;
      'parallel' )
        if [ -z "$DO_PARALLEL" ] ; then
          echo "  $DESCRIP must be run in parallel."
          ((CHECKERR++))
        fi
        ;;
      'maxthreads' )
        shift
        if [ ! -z "$DO_PARALLEL" ] ; then
          if [ $N_THREADS -gt $1 ] ; then
            echo "  $DESCRIP can only run with $1 or fewer parallel threads."
            ((CHECKERR++))
          fi
        fi
        ;;
      'nthreads' )
        shift
        if [ ! -z "$DO_PARALLEL" ] ; then
          REMAINDER=`expr $N_THREADS % $1`
          if [ -z "$REMAINDER" -o $REMAINDER -ne 0 ] ; then
            echo "  $DESCRIP requires a multiple of $1 parallel threads."
            ((CHECKERR++))
          fi
        fi
        ;;
      'threads' )
        shift
        if [ ! -z "$DO_PARALLEL" ] ; then
          if [ $N_THREADS -ne $1 ] ; then
            echo "  $DESCRIP requires exactly $1 parallel threads."
            ((CHECKERR++))
          fi
        fi
        ;;
      'amberhome' )
        if [ -z "$AMBERHOME" ] ; then
          echo "  $DESCRIP requires AMBERHOME to be set."
          ((CHECKERR++))
        fi
        ;;
      'amberorcpptraj' )
        if [ -z "$CPPTRAJHOME" -a -z "$AMBERHOME" ] ; then
          echo "  $DESCRIP requires CPPTRAJHOME or AMBERHOME to be set."
          ((CHECKERR++))
        fi
        ;;
      'inpath' )
        shift
        if [ -z "`which $1`" ] ; then
          echo "  $DESCRIP requires $1 to be in PATH."
          ((CHECKERR++))
        fi
        ;;
      'testos' )
        shift
        if [ "$CPPTRAJ_TEST_OS" != "$1" ] ; then
          echo "  $DESCRIP requires $1 OS."
          ((CHECKERR++))
        fi
        ;;
      'notos' )
        shift
        if [ "$CPPTRAJ_TEST_OS" = "$1" ] ; then
          echo "  $DESCRIP cannot run on $1 OS."
          ((CHECKERR++))
        fi
        ;;
      'long' )
        if [ -z "$CPPTRAJ_LONG_TEST" -o $CPPTRAJ_LONG_TEST -eq 0 ] ; then
          echo "  $DESCRIP is a long test and long tests disabled. Use 'long' to run."
          ((CHECKERR++))
        fi
        ;;
      'file' )
        shift
        if [ ! -f "$1" ] ; then
          echo "  $DESCRIP requires file $1"
          ((CHECKERR++))
        fi
        ;;
      'disabled' )
        echo "  $DESCRIP is disabled."
        ((CHECKERR++))
        ;;
      * ) ErrMsg "Unknown CheckEnv() option: $1" ; exit 1 ;;
    esac
    shift
  done
}

# Requires() <list>
# If list of requirements is not met, skip entire test.
Requires() {
  DESCRIP=$TESTNAME
  CheckEnv $*
  if [ "$CHECKERR" -ne 0 ] ; then
    SkipTest "$TESTNAME"
  fi
}

# CheckFor() <list>
# \return 1 if unit should be skipped, 0 otherwise.
CheckFor() {
  DESCRIP=$UNITNAME
  CheckEnv $*
  if [ "$CHECKERR" -ne 0 ] ; then
    SkipCheck "$UNITNAME"
    return 1
  fi
  return 0
}

# ------------------------------------------------------------------------------
# TODO remove all deprecated functions below.
# FIXME This is a stub and should be removed
CheckTest() {
  ErrMsg "CHECKTEST DISABLED!"
  exit 1
}

Disabled() {
  ErrMsg "FUNCTION DISABLED. FIX IT!"
  exit 1
}

CheckZlib() {
  Disabled
}

CheckBzlib() {
  Disabled
}

CheckNetcdf() {
  Disabled
}

RequiresNetcdf() {
  Disabled
}

CheckXdr() {
  Disabled
}

RequiresXdr() {
  Disabled
}

CheckMathlib() {
  Disabled
}

CheckPtrajAnalyze() {
  Disabled
}

RequiresMathlib() {
  Disabled
}

CheckSanderlib() {
  Disabled
}

RequiresSanderlib() {
  Disabled
}

CheckPnetcdf() {
  Disabled
}

NotParallel() {
  Disabled
}

RequiresNotParallel() {
  Disabled
}

CheckNthreads() {
  Disabled
}

MaxThreads() {
  Disabled
}

RequiresMaxThreads() {
  Disabled
}

RequiresThreads() {
  Disabled
}

# ==============================================================================
# M A I N
# ==============================================================================
#echo "DEBUG: Begin MasterTest.sh. $*"
#echo "DEBUG: CPPTRAJ_TEST_MODE: $CPPTRAJ_TEST_MODE"
if [ -z "$CPPTRAJ_TEST_SETUP" ] ; then
  #echo "DEBUG: Initial test setup."
  # MasterTest.sh has not been called yet; set up test environment.
  export CPPTRAJ_TEST_ROOT=`pwd`
  # If invocation is "./CpptrajTest.sh", all tests. Otherwise assume individual test dir.
  if [ "$0" != './CpptrajTest.sh' ] ; then
    CPPTRAJ_TEST_ROOT=`dirname $CPPTRAJ_TEST_ROOT`
  fi
  #echo "DEBUG: CPPTRAJ_TEST_ROOT= $CPPTRAJ_TEST_ROOT $0"
  # If CPPTRAJ_TEST_OS is not set, try to determine.
  if [ -z "$CPPTRAJ_TEST_OS" ] ; then
    export CPPTRAJ_TEST_OS=`uname -s | awk '{print $1}'`
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
  Required "grep"
  Required "sed"
  Required "awk"
  Required "rmdir"
  # Set some defaults
  export CPPTRAJ_TEST_RESULTS='Test_Results.dat'
  export CPPTRAJ_TEST_ERROR='Test_Error.dat'
  CPPTRAJ_OUTPUT='test.out'
  CPPTRAJ_ERROR=''
  #export CPPTRAJ_RNG='--rng marsaglia'
  CPPTRAJ_RNG='' # Placeholder
  # Process command line options
  CmdLineOpts $*
  # Determine standalone or AmberTools
  if [ ! -z "`echo "$CPPTRAJ_TEST_ROOT" | grep AmberTools`" ] ; then
    # Assume AmberTools. Need AMBERHOME.
    STANDALONE=0
    if [ -z "$AMBERHOME" ] ; then
      ErrMsg "In AmberTools and AMBERHOME is not set. Required for tests."
      exit 1
    fi
  else
    # Standalone. Never use dacdif.
    STANDALONE=1
    USE_DACDIF=0
  fi
  # Determine if diff or dacdif style will be used.
  CPPTRAJ_DIFF=''
  CPPTRAJ_DACDIF=''
  if [ $USE_DACDIF -eq 1 ] ; then
    CPPTRAJ_DACDIF="$AMBERHOME/test/dacdif"
    if [ ! -f "$CPPTRAJ_DACDIF" ] ; then
      echo "$CPPTRAJ_DACDIF not found. Required for AmberTools tests."
      exit 1
    fi
    export CPPTRAJ_TEST_ERROR="$AMBERHOME/test/TEST_FAILURES.diff"
  fi
  #else
    Required "diff"
    CPPTRAJ_DIFF=`which diff`
  #fi
  export CPPTRAJ_DACDIF
  export CPPTRAJ_DIFF
  # If not cleaning see what else needs to be set up.
  if [ $CPPTRAJ_TEST_CLEAN -eq 0 -a -z "$IS_LIBCPPTRAJ" ] ; then
    # Determine binary locations
    SetBinaries
    if [ ! -z "$DIFFOPTS" ] ; then
      echo "Warning: DIFFOPTS is set to '$DIFFOPTS'"
    fi
  fi # END if not cleaning
  # Export test output and error file names
  export CPPTRAJ_OUTPUT
  export CPPTRAJ_ERROR
  # Initial setup complete
  #echo "DEBUG: Initial test setup complete."
  export CPPTRAJ_TEST_SETUP='yes'
fi # END initial setup

# Determine mode of execution: individual test or multiple tests.
if [ "$CPPTRAJ_TEST_MODE" = 'master' ] ; then
  # Executed from CpptrajTest.sh. Assume we are executing multiple
  # tests.
  #echo "DEBUG: Executing multiple tests."
  # Clean any existing test results files.
  if [ -f "$CPPTRAJ_TEST_RESULTS" ] ; then
    $CPPTRAJ_RM $CPPTRAJ_TEST_RESULTS
  fi
  if [ -f "$CPPTRAJ_TEST_ERROR" -a -z "$CPPTRAJ_DACDIF" ] ; then
    $CPPTRAJ_RM $CPPTRAJ_TEST_ERROR
  fi
  if [ ! -z "$TEST_DIRS" ] ; then
    #echo "DEBUG: Running tests in specified directories."
    if [ $GET_TIMING -eq 2 ] ; then
      TIMING_FILE="$CPPTRAJ_TEST_ROOT/Timing.dat"
      if [ -f "$TIMING_FILE" ] ; then
        $CPPTRAJ_RM $TIMING_FILE
      fi
      for DIR in $TEST_DIRS ; do
        cd $CPPTRAJ_TEST_ROOT/$DIR
        $CPPTRAJ_TIME -f "%e" -o testtime ./RunTest.sh
        printf "%.2f %s\n" `cat testtime` $DIR >> $TIMING_FILE
        $CPPTRAJ_RM testtime
      done
      echo "Test timing data written to $TIMING_FILE"
    elif [ $EXIT_ON_ERROR -ne 0 ] ; then
      TEST_OK_FILE=$CPPTRAJ_TEST_ROOT/OK
      if [ ! -f "$TEST_OK_FILE" ] ; then
        touch $TEST_OK_FILE
      fi
      for DIR in $TEST_DIRS ; do
        if [ -z "`grep $DIR $TEST_OK_FILE`" ] ; then
          cd $CPPTRAJ_TEST_ROOT/$DIR && ./RunTest.sh
          if [ $EXIT_ON_ERROR -ne 0 -a $? -ne 0 ] ; then
            break
          fi
          echo "$DIR" >> $TEST_OK_FILE
        fi
      done
      echo "DEBUG: OK tests listed in $TEST_OK_FILE"
    else
      for DIR in $TEST_DIRS ; do
        cd $CPPTRAJ_TEST_ROOT/$DIR && ./RunTest.sh
      done
    fi
    if [ $SUMMARY -ne 0 ] ; then
      cd $CPPTRAJ_TEST_ROOT
      Summary
    fi
  else
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
    Required "make"
    if [ -z "$TARGET" ] ; then
      # If no target specified, probably not executed via make.
      ErrMsg "No test directories specified."
      exit 1
    fi
    make $TARGET
    if [ $? -ne 0 ] ; then
      exit 1
    fi
  fi
else
  # Assume we are only executing a single test.
  #echo "DEBUG: Executing single test."
  # Single test.
  # Always clean up individual test output and error files
  if [ ! -z "$CPPTRAJ_OUTPUT" -a -f "$CPPTRAJ_OUTPUT" ] ; then
    $CPPTRAJ_RM $CPPTRAJ_OUTPUT
  fi
  if [ ! -z "$CPPTRAJ_ERROR" -a -f "$CPPTRAJ_ERROR" ] ; then
    $CPPTRAJ_RM $CPPTRAJ_ERROR
  fi
  if [ -f "$CPPTRAJ_TEST_RESULTS" ] ; then
    $CPPTRAJ_RM $CPPTRAJ_TEST_RESULTS
  fi
  if [ -z "$CPPTRAJ_DACDIF" ] ; then
    # Standalone. Only remove previous error file if it exists.
    if [ -f "$CPPTRAJ_TEST_ERROR" ] ; then
      $CPPTRAJ_RM $CPPTRAJ_TEST_ERROR
    fi
  else
    # AmberTools - remove previous .dif files
    $CPPTRAJ_RM *.dif 2> tmp.cpptrajtest.devnull
    rm tmp.cpptrajtest.devnull
  fi
  if [ -f 'valgrind.out' ] ; then
    $CPPTRAJ_RM valgrind.out
  fi
  if [ -f 'test.out' ] ; then
    $CPPTRAJ_RM test.out
  fi
  THREADFILES=`ls Thread.* 2> tmp.cpptrajtest.devnull`
  rm tmp.cpptrajtest.devnull
  if [ ! -z "$THREADFILES" ] ; then
    for FILE in $THREADFILES ; do
      $CPPTRAJ_RM $FILE
    done
  fi
  if [ "$CPPTRAJ_TEST_CLEAN" -eq 0 ] ; then
    TEST_WORKDIR=`pwd`
    TestHeader
    if [ -z "$CPPTRAJ_DACDIF" ] ; then
      TestHeader "$CPPTRAJ_TEST_RESULTS"
    fi
  fi
fi

#echo "DEBUG: End MasterTest.sh."
