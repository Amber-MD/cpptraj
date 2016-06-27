#!/bin/bash

. ../MasterTest.sh

CleanFiles cluster.in random.out random.info.dat random.summary.dat \
           random.crd.c? random.cpop.agr random.rep.*.pdb CpptrajPairDist
CheckNetcdf
PREFIX="random"
INPUT="-i cluster.in"

Ctest() {
  cat > cluster.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
cluster crd1 @CA clusters 5 rms out $PREFIX.out includesieved_cdist \
        sieve 5 random sieveseed 1 includesieveincalc \
        summary $PREFIX.summary.dat info $PREFIX.info.dat bestrep cumulative \
        clusterout $PREFIX.crd cpopvtime $PREFIX.cpop.agr \
        repout $PREFIX.rep repfmt pdb savepairdist loadpairdist
EOF
  RunCpptraj "Cluster test with random sieve ($1)."
  DoTest random.info.dat.save random.info.dat
  DoTest random.summary.dat.save random.summary.dat
}

Ctest "Save Pairwise Distances"
Ctest "Load Pairwise Distances"

EndTest
exit 0
