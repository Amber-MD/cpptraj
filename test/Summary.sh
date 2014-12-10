#!/bin/bash

. MasterTest.sh

if [[ $CLEAN -eq 1 ]] ; then
  echo Clean.
  exit 0
fi
RESULTFILES=""
if [[ ! -z $TEST_RESULTS ]] ; then
  RESULTFILES=`ls */$TEST_RESULTS`
else
  exit 0
fi
echo "===================== TEST SUMMARY ======================"
if [[ ! -z $RESULTFILES ]] ; then
  cat $RESULTFILES > $TEST_RESULTS
  # DoTest - Number of comparisons OK
  OK=`cat $TEST_RESULTS | grep OK | wc -l`
  # DoTest - Number of comparisons different
  ERR=`cat $TEST_RESULTS | grep different | wc -l`
  NOTFOUND=`cat $TEST_RESULTS | grep "not found" | wc -l`
  ((ERR = $ERR + $NOTFOUND))
  # Number of tests run
  NTESTS=`cat $TEST_RESULTS | grep "TEST:" | wc -l`
  # Number of tests successfully finished
  PASSED=`cat $TEST_RESULTS | grep "comparisons passed" | wc -l`
  #CERR=`cat $RESULTFILES | grep failed | awk 'BEGIN{sum=0;}{sum+=$1;}END{print sum}'`
  #NOERR=`cat $RESULTFILES | grep failed | awk 'BEGIN{sum=0;}{sum+=$4;}END{print sum}'`
  #PASSED=`cat $RESULTFILES | grep passed | wc -l`
  ((NCOMPS = $OK + $ERR))
  echo "  $OK out of $NCOMPS comparisons OK ($ERR failed)."
  echo "  $PASSED out of $NTESTS tests completed with no issues."
  #echo "  $PASSED tests passed."
  #echo "  $ERR tests failed."
  RESULTFILES=`ls */$TEST_ERROR`
  if [[ ! -z $RESULTFILES ]] ; then
    cat $RESULTFILES > $TEST_ERROR
  fi 
else
  echo "No Test Results files (./*/$TEST_RESULTS) found."
fi

if [[ ! -z $VALGRIND ]] ; then
  RESULTFILES=`ls */$ERROR`
  if [[ ! -z $RESULTFILES ]] ; then
    echo "---------------------------------------------------------"
    echo "Valgrind summary:"
    NUMVGERR=`cat $RESULTFILES | grep ERROR | awk 'BEGIN{sum=0;}{sum+=$4;}END{print sum;}'`
    echo "    $NUMVGERR errors."
    NUMVGOK=`cat $RESULTFILES | grep "All heap" | wc -l`
    echo "    $NUMVGOK memory leak checks OK."
    NUMVGLEAK=`cat $RESULTFILES | grep LEAK | wc -l`
    echo "    $NUMVGLEAK memory leak reports."
#    echo "---------------------------------------------------------"
  else
    echo "No valgrind test results found."
    exit 0
  fi    
fi 
echo "========================================================="

exit 0
