#!/bin/bash

. ../MasterTest.sh

CleanFiles cluster.in CpptrajPairDist *.rmsd.dat *.summary.dat *.info.dat \
           *.half.dat *.nc.c? *.pdb *.agr *.out *.sil.*.dat

INPUT="-i cluster.in"
CheckNetcdf
Cluster() {
  PREFIX=$1
  SIEVEARG=$2
  SAVEARG=$3
  cat > cluster.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
#debug analysis 2
cluster crd1 @CA clusters 5 rms out $PREFIX.out summary $PREFIX.summary.dat info $PREFIX.info.dat \
        clusterout $PREFIX.nc clusterfmt netcdf summaryhalf $PREFIX.half.dat cpopvtime $PREFIX.cpop.agr splitframe 24 \
        repout $PREFIX.rep.pdb repfmt pdb $SIEVEARG $SAVEARG sil $PREFIX.sil
EOF
  RunCpptraj "Cluster, $PREFIX $SAVEARG"
  DoTest $PREFIX.info.dat.save $PREFIX.info.dat
  DoTest $PREFIX.half.dat.save $PREFIX.half.dat
  DoTest $PREFIX.sil.cluster.dat.save $PREFIX.sil.cluster.dat
}

SIEVEARGS="sieve 5 bestrep cumulative includesieveincalc"

# Test sieving
Cluster nosieve
Cluster sieve5 "$SIEVEARGS"

# Test loading/saving of pairdist file with/without sieve
Cluster nosieve " " savepairdist
Cluster nosieve " " loadpairdist
Cluster sieve5 "$SIEVEARGS" savepairdist
NcTest sieve5.nc.c0.save sieve5.nc.c0
Cluster sieve5 "$SIEVEARGS" loadpairdist

# Test pairwise no memory
Cluster nosieve " "       "pairwisecache none"
Cluster sieve5  "$SIEVEARGS" "pairwisecache none"

# Test random sieving
Cluster random "$SIEVEARGS random sieveseed 1"

EndTest

exit 0
