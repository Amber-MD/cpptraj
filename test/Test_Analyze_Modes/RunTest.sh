#!/bin/bash

. ../MasterTest.sh

CleanFiles modes.in fluct.dat displ.dat corr.dat modestest.2.crd eigenval.dat rmsip.dat

INPUT='modes.in'
CheckMathlib
if [ $? -ne 0 ] ; then
  SkipTest "Modes Analysis"
fi

# Test modes fluct and mwcovar matrix generation
TestFluct() {
  TESTNAME='Modes analysis, RMS fluctuations'
  CheckNetcdf "$TESTNAME"
  if [ $? -ne 0 ] ; then
    SkipCheck "$TESTNAME"
  else
    TOP=../tz2.parm7
    cat > modes.in <<EOF
trajin ../tz2.nc
matrix mwcovar name tz2 @CA
diagmatrix tz2 name tz2modes vecs 20
modes fluct name tz2modes out fluct.dat setname Fluct prec 10.3
EOF
    RunCpptraj "$TESTNAME"
    DoTest fluct.dat.save fluct.dat
  fi
}

# Test modes displ and modes file read
TestDispl() {
  TOP=../tz2.parm7
  # Since the displacement test can be fooled by eigenvector
  # sign flips, read in a previously generated set of modes.
  cat > modes.in <<EOF
readdata tz2.evecs.dat name tz2modes
analyze modes displ name tz2modes out displ.dat setname Displ prec 10.3
EOF
  RunCpptraj "Modes analysis, displacements"
  DoTest displ.dat.save displ.dat
}

# Test modes corr and mwcovar matrix generation
TestCorr() {
  TESTNAME='Modes analysis, dipole correlation'
  CheckNetcdf "$TESTNAME"
  if [ $? -ne 0 ] ; then
    SkipCheck "$TESTNAME"
  else
    TOP=../tz2.parm7
    cat > modes.in <<EOF
trajin ../tz2.nc
matrix mwcovar name tz2
analyze matrix tz2 name tz2modes vecs 20
analyze modes corr name tz2modes out corr.dat mask1 :2-13@N mask2 :2-13@H
EOF
    RunCpptraj "$TESTNAME"
    DoTest corr.dat.save corr.dat
  fi
}

# Test modes trajout
TestTrajout() {
  TOP="INPpYLYP.FF14SB.parm7"
  cat > modes.in <<EOF
readdata evecs.dat name evecs
modes name evecs trajout modestest.2.crd pcmin -33 pcmax 46 tmode 2 trajoutmask !@H=
EOF
  RunCpptraj "Modes analysis, pseudo trajectory creation test"
  DoTest modestest.2.crd.save modestest.2.crd
}

# Test modes RMSIP and eigenval
TestRMSIP() {
  cat > modes.in <<EOF
readdata evecs.dat name evecs
readdata evecs2.dat name evecs2
modes name evecs name2 evecs2 rmsip beg 1 end 5 out rmsip.dat setname RMSIP
modes name evecs eigenval out eigenval.dat setname EV prec 12.6
EOF
  RunCpptraj "Modes analysis, RMSIP and eigenvalue fraction test"
  DoTest rmsip.dat.save rmsip.dat
  DoTest eigenval.dat.save eigenval.dat
}

TestFluct
TestDispl
TestCorr
TestTrajout
TestRMSIP

EndTest

exit 0
