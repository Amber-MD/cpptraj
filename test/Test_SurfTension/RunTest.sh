#!/bin/bash
# Smoke tests for the surftension Action.
# tz2.ortho is a solvated protein box, not a liquid slab, so the 1-frame run
# only checks that Init/Setup do not crash. Numeric γ comparison needs a
# dedicated slab trajectory.

. ../MasterTest.sh

CleanFiles st.in st.summary.dat

TESTNAME='Surface tension (surftension) tests'

# Command is registered and Help() prints.
UNITNAME='surftension help'
cat > st.in <<EOF
help surftension
EOF
INPUT='-i st.in'
RunCpptraj "$UNITNAME"

# One ortho frame: parse keywords, bind :WAT@O, require a box.
UNITNAME='surftension smoke (Init/Setup)'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > st.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
surftension :WAT@O temp 300.0
run
EOF
  INPUT='-i st.in'
  RunCpptraj "$UNITNAME"
fi

UNITNAME='surftension smoke (interface itim)'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > st.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
surftension :WAT@O temp 300.0 interface itim
run
EOF
  INPUT='-i st.in'
  RunCpptraj "$UNITNAME"
fi

UNITNAME='surftension smoke (normal x)'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > st.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
surftension :WAT@O temp 300.0 normal x
run
EOF
  INPUT='-i st.in'
  RunCpptraj "$UNITNAME"
fi

UNITNAME='surftension smoke (normal y, dnormal alias)'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > st.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
surftension :WAT@O temp 300.0 normal y dnormal 1.0 sigmanormal 1.5
run
EOF
  INPUT='-i st.in'
  RunCpptraj "$UNITNAME"
fi

UNITNAME='surftension smoke (nsurf 1)'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > st.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
surftension :WAT@O temp 300.0 nsurf 1 side upper
run
EOF
  INPUT='-i st.in'
  RunCpptraj "$UNITNAME"
fi

UNITNAME='surftension smoke (mask2, summaryout, blocktime)'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > st.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
surftension :WAT@O mask2 :WAT@O temp 300.0 dt 100 blocktime 50000 summaryout st.summary.dat
run
EOF
  INPUT='-i st.in'
  RunCpptraj "$UNITNAME"
fi

EndTest
exit 0
