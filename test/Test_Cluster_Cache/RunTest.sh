#!/bin/bash

. ../MasterTest.sh

CleanFiles cluster.in *.cnumvtime.dat *.info.dat *.summary.dat PW0 PW1 \
           CpptrajPairwiseCache

INPUT="-i cluster.in"
TESTNAME='Cluster pairwise cache tests'
Requires netcdf

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

# Test in-memory cache save/load with no sieve
Cluster nosieve.mem.save " " "savepairdist pairdist PW0 pairwisecache mem"
Cluster nosieve.mem.load " " "loadpairdist pairdist PW0"
DoTest nosieve.mem.save.info.dat.save nosieve.mem.save.info.dat
DoTest nosieve.mem.save.cnumvtime.dat nosieve.mem.load.cnumvtime.dat
DoTest nosieve.mem.save.info.dat      nosieve.mem.load.info.dat
DoTest nosieve.mem.save.summary.dat   nosieve.mem.load.summary.dat

# Test on-disk cache with no sieve
Cluster nosieve.disk.save " " "savepairdist pairdist PW1 pairwisecache disk"
DoTest nosieve.mem.save.info.dat.save nosieve.disk.save.info.dat
Cluster nosieve.disk.load " " "loadpairdist pairdist PW1"
DoTest nosieve.mem.save.info.dat.save nosieve.disk.load.info.dat

# Test no cache, no sieve
Cluster nosieve.nocache.save " " "pairwisecache none"
DoTest nosieve.mem.save.info.dat.save nosieve.nocache.save.info.dat

SIEVEARGS="sieve 5"
# Test in-memory cache save/load with sieve
Cluster sieve5.mem.save "$SIEVEARGS" "savepairdist pairdist PW0 pairwisecache mem"
Cluster sieve5.mem.load "$SIEVEARGS" "loadpairdist pairdist PW0"
DoTest sieve5.mem.save.info.dat.save sieve5.mem.save.info.dat
DoTest sieve5.mem.save.cnumvtime.dat sieve5.mem.load.cnumvtime.dat
DoTest sieve5.mem.save.info.dat      sieve5.mem.load.info.dat
DoTest sieve5.mem.save.summary.dat   sieve5.mem.load.summary.dat

# Test on-disk cache with sieve
Cluster sieve5.disk.save "$SIEVEARGS" "savepairdist pairdist PW1 pairwisecache disk"
DoTest sieve5.mem.save.info.dat.save sieve5.disk.save.info.dat

# Test no cache, with sieve
Cluster sieve5.nocache.save "$SIEVEARGS" "pairwisecache none"
DoTest sieve5.mem.save.info.dat.save sieve5.nocache.save.info.dat

# Test loading pairwise cache with different sieve
Cluster sieve3.mem.load "sieve 3" "loadpairdist pairdist PW1 pwrecalc"
DoTest sieve3.mem.load.info.dat.save sieve3.mem.load.info.dat

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
