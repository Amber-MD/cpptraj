#!/bin/bash

if [ -z "$CPPTRAJ" ] ; then
  #CPPTRAJ=`which cpptraj`
  export OMP_NUM_THREADS=4
  CPPTRAJ=`which cpptraj.OMP`
fi

TOP=~/Cpptraj/Cpptraj-ExtendedTests/ChainA_1-268_NAD_TCL-gaff.tip3p.parm7
TRJ=~/Cpptraj/Cpptraj-ExtendedTests/run9.nc # 2000 frames
COUNT=2000

# Each trajectory is 100000 frames
if [ ! -f '../framelist.sh' ] ; then
  echo "no ../framelist.sh"
  exit 1
fi
source ../framelist.sh

for TOTAL in $FRAMELIST ; do
  if [ -f 'CpptrajPairDist' ] ; then
    rm CpptrajPairDist
  fi
  cat > cpptraj.in <<EOF
set PREFIX = test1
parm $TOP
EOF
  # Set trajin
  echo " ===== $TOTAL ===== "
  if [ $TOTAL -le $COUNT ] ; then
    echo "trajin $TRJ 1 $TOTAL 1" >> cpptraj.in
  else
    ((NTRAJ = $TOTAL / $COUNT))
    ((REM   = $TOTAL % $COUNT))
    if [ $REM -gt 0 ] ; then
      ((NTRAJ++))
    fi
    echo "$NTRAJ $REM"
    for (( i=1; i < $NTRAJ; i++ )) ; do
      echo "trajin $TRJ 1 $COUNT 1" >> cpptraj.in
    done
    if [ $REM -gt 0 ] ; then
      echo "trajin $TRJ 1 $REM 1" >> cpptraj.in
    else
      echo "trajin $TRJ 1 $COUNT 1" >> cpptraj.in
    fi
  fi
  
  # Get number of frames we want from any given trajectory
  ((NFRAMES = $TOTAL / 4))
  # Calculate split frames
  ((SPLIT1 = $NFRAMES * 1))
  ((SPLIT2 = $NFRAMES * 2))
  ((SPLIT3 = $NFRAMES * 3))
  echo "Total $TOTAL, Each section $NFRAMES splits $SPLIT1 $SPLIT2 $SPLIT3"

#  kmeans clusters 5 randompoint kseed 42 \
  cat >> cpptraj.in <<EOF

cluster \
  hieragglo clusters 5 \
  rms :10-260@CA sieve 20 random sieveseed 23 savepairdist \
  summaryhalf \$PREFIX.split.dat splitframe $SPLIT1,$SPLIT2,$SPLIT3 \
  info \$PREFIX.info.dat \
  summary \$PREFIX.summary.dat \
  sil \$PREFIX.Sil \
  cpopvtime \$PREFIX.cpop.agr normframe \
  bestrep cumulative_nosieve
EOF
  $CPPTRAJ -i cpptraj.in -o cpptraj.$TOTAL.out
done
