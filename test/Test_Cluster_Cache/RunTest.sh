#!/bin/bash

. ../MasterTest.sh

CleanFiles cluster.in *.cnumvtime.dat *.info.dat *.summary.dat PW0

INPUT="-i cluster.in"
TESTNAME='Cluster pairwise cache tests'

# <prefix> <sieve> <save>
Cluster() {
  PREFIX=$1
  if [ "$PREFIX" = 'random' ] ; then
    RNG='rng setdefault marsaglia'
  else
    RNG=''
  fi
  SIEVEARG=$2
  SAVEARG=$3
  cat > cluster.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
strip !@CA
$RNG

cluster C1 \
        hieragglo clusters 5 rms $SIEVEARG $SAVEARG \
        out     $PREFIX.cnumvtime.dat \
        summary $PREFIX.summary.dat \
        info    $PREFIX.info.dat
EOF
  RunCpptraj "Cluster, $PREFIX $SAVEARG"
  #DoTest $PREFIX.cnumvtime.dat.save $PREFIX.cnumvtime.dat
  #DoTest $PREFIX.info.dat.save      $PREFIX.info.dat
  #DoTest $PREFIX.summary.dat.save   $PREFIX.summary.dat
}

SIEVEARGS="sieve 5 bestrep cumulative includesieveincalc"

# Test in-memory cache save/load with no sieve
Cluster nosieve.mem.save " " "savepairdist pairdist PW0 pairwisecache mem"
Cluster nosieve.mem.load " " "loadpairdist pairdist PW0"
DoTest nosieve.mem.save.info.dat.save nosieve.mem.save.info.dat
DoTest nosieve.mem.save.cnumvtime.dat nosieve.mem.load.cnumvtime.dat
DoTest nosieve.mem.save.info.dat      nosieve.mem.load.info.dat
DoTest nosieve.mem.save.summary.dat   nosieve.mem.load.summary.dat

# Test sieving
#Cluster nosieve
#Cluster sieve5 "$SIEVEARGS"

# Test loading/saving of pairdist file with/without sieve
#Cluster nosieve " " savepairdist
#Cluster nosieve " " loadpairdist
#Cluster sieve5 "$SIEVEARGS" savepairdist
#NcTest sieve5.nc.c0.save sieve5.nc.c0
#Cluster sieve5 "$SIEVEARGS" loadpairdist

# Test pairwise no memory
#Cluster nosieve " "       "pairwisecache none"
#Cluster sieve5  "$SIEVEARGS" "pairwisecache none"

# Test random sieving
#Cluster random "$SIEVEARGS random sieveseed 1"

EndTest

exit 0
