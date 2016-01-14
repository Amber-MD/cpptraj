#!/bin/bash

. ../MasterTest.sh

CleanFiles check.in report.dat

MaxThreads 1 "Structure check test"
if [[ $? -ne 0 ]] ; then
  EndTest
  exit 0
fi

INPUT="-i check.in"
cat > check.in <<EOF
parm ../tz2.parm7
trajin tz2.stretched.pdb
check reportfile report.dat offset 0.7
EOF
RunCpptraj "Structure Check"
if [[ ! -z $OPENMP ]] ; then
  # OpenMP check does not bother sorting for efficiency
  SORT=`which sort`
  if [[ -z $SORT ]] ; then
    echo "Warning: 'sort' required for Structure Check OpenMP test. Skipping."
    exit 0
  fi
  sort report.dat.save > sortedReport.save
  sort report.dat      > sortedReport
  DoTest sortedReport.save sortedReport
  CleanFiles sortedReport.save sortedReport
else
  DoTest report.dat.save report.dat
fi
CheckTest
EndTest

exit 0
