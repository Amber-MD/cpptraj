#!/bin/bash
<<<<<<< HEAD
# surftension tests.
# slab.pdb is a 20 A cubic lattice of O atoms with a small z-corrugation
# so ITIM min/max is not flat. Stretching one box length creates vacuum
# on both sides of a slab so Willard/ITIM can find interfaces.
# The PDB is loaded twice so serial and MPI (2 ranks) both analyze 2
# identical frames; SyncAction then reduces |h_q|^2 and the summaries match.

. ../MasterTest.sh

CleanFiles st.in st.willard.dat st.itim.dat st.normalx.dat st.normaly.dat \
           st.nsurf1.dat st2_summary.dat

TESTNAME='Surface tension (surftension) tests'
# Two frames: Jenkins mpiexec -n 2 assigns one frame per rank.
Requires maxthreads 2

INPUT='-i st.in'

UNITNAME='surftension Willard-Chandler (default)'
cat > st.in <<EOF
noprogress
parm slab.pdb
trajin slab.pdb
trajin slab.pdb
box z 50.0
surftension :WAT@O temp 300.0 qmax 0.80 summaryout st.willard.dat
run
EOF
RunCpptraj "$UNITNAME"
DoTest st.willard.dat.save st.willard.dat -a 0.001

UNITNAME='surftension ITIM'
cat > st.in <<EOF
noprogress
parm slab.pdb
trajin slab.pdb
trajin slab.pdb
box z 50.0
surftension :WAT@O temp 300.0 qmax 0.80 interface itim summaryout st.itim.dat
run
EOF
RunCpptraj "$UNITNAME"
DoTest st.itim.dat.save st.itim.dat -a 0.001

UNITNAME='surftension normal x'
cat > st.in <<EOF
noprogress
parm slab.pdb
trajin slab.pdb
trajin slab.pdb
box x 50.0
surftension :WAT@O temp 300.0 qmax 0.80 normal x summaryout st.normalx.dat
run
EOF
RunCpptraj "$UNITNAME"
DoTest st.normalx.dat.save st.normalx.dat -a 0.001

UNITNAME='surftension normal y, dnormal alias'
cat > st.in <<EOF
noprogress
parm slab.pdb
trajin slab.pdb
trajin slab.pdb
box y 50.0
surftension :WAT@O temp 300.0 qmax 0.80 normal y dnormal 1.0 sigmanormal 1.5 summaryout st.normaly.dat
run
EOF
RunCpptraj "$UNITNAME"
DoTest st.normaly.dat.save st.normaly.dat -a 0.001

UNITNAME='surftension nsurf 1'
cat > st.in <<EOF
noprogress
parm slab.pdb
trajin slab.pdb
trajin slab.pdb
box z 50.0
surftension :WAT@O temp 300.0 qmax 0.80 nsurf 1 side upper summaryout st.nsurf1.dat
run
EOF
RunCpptraj "$UNITNAME"
DoTest st.nsurf1.dat.save st.nsurf1.dat -a 0.001

UNITNAME='surftension mask2, fprefix, blocktime'
cat > st.in <<EOF
noprogress
parm slab.pdb
trajin slab.pdb
trajin slab.pdb
box z 50.0
surftension :WAT@O mask2 :WAT@O temp 300.0 qmax 0.80 dt 1.0 blocktime 10.0 \
  fprefix st2_ summaryout summary.dat
run
EOF
RunCpptraj "$UNITNAME"
DoTest st2_summary.dat.save st2_summary.dat -a 0.001
=======
# Smoke tests for the surftension Action.
# tz2.ortho is a solvated protein box, not a liquid slab, so the 1-frame run
# only checks that Init/Setup do not crash. Numeric γ comparison needs a
# dedicated slab trajectory.

. ../MasterTest.sh

CleanFiles st.in

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
>>>>>>> c51400c7 (Add 'surftension' command to calculate capillary-wave surface tension of a liquid slab)

EndTest
exit 0
