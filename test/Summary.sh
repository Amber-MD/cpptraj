#!/bin/bash

. MasterTest.sh

if [[ $CLEAN -eq 1 ]] ; then
  echo Clean.
  exit 0
fi

RESULTFILES=`ls */$TEST_RESULTS`
if [[ ! -z $RESULTFILES ]] ; then
  ERR=`cat $RESULTFILES | grep failed | wc -l`
  CERR=`cat $RESULTFILES | grep failed | awk 'BEGIN{sum=0;}{sum+=$1;}END{print sum}'`
  NOERR=`cat $RESULTFILES | grep failed | awk 'BEGIN{sum=0;}{sum+=$4;}END{print sum}'`
  OK=`cat $RESULTFILES | grep OK | wc -l`
  PASSED=`cat $RESULTFILES | grep passed | wc -l`
  cat $RESULTFILES > $TEST_RESULTS
else
  echo "No Test Results found."
  exit 0
fi

echo "===================== TEST SUMMARY ======================"
echo "  $OK comparisons OK."
echo "  $CERR comparisons failed."
echo "  $PASSED tests passed."
echo "  $ERR tests failed."

if [[ ! -z $VALGRIND ]] ; then
  RESULTFILES=`ls */$ERROR`
  if [[ ! -z $RESULTFILES ]] ; then
    echo "---------------------------------------------------------"
    echo "Valgrind summary:"
    NUMVGERR=`cat $RESULTFILES | grep ERROR | awk 'BEGIN{sum=0;}{sum+=$4;}END{print sum;}'`
    echo "    $NUMVGERR errors."
    NUMVGOK=`cat $RESULTFILES | grep heap | wc -l`
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
