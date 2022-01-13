#!/bin/bash

CPPTRAJ=`which cpptraj`
#export OMP_NUM_THREADS=12
#CPPTRAJ=`which cpptraj.OMP`

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
  continue
  
  # Get number of frames we want from any given trajectory
  ((NFRAMES = $TOTAL / 4))
  # Calculate split frames
  ((SPLIT1 = $NFRAMES * 1))
  ((SPLIT2 = $NFRAMES * 2))
  ((SPLIT3 = $NFRAMES * 3))
  echo "Total $TOTAL, Each traj $NFRAMES splits $SPLIT1 $SPLIT2 $SPLIT3"

  cat > cpptraj.in <<EOF
set PREFIX = test1

parm /home/bergonzoc/WORK/HBV_rna/3ANALYSIS/nowat.noions.parm7
trajin ~/WORK/HBV_rna/2PROD/analysis/nowat.nc.0 1 $NFRAMES 1
trajin ~/WORK/HBV_rna/2PROD/analysis/nowat.nc.1 1 $NFRAMES 1
trajin ~/WORK/HBV_rna/2PROD/analysis/nowat.nc.2 1 $NFRAMES 1
trajin ~/WORK/HBV_rna/2PROD/analysis/nowat.nc.3 1 $NFRAMES 1


cluster \
  kmeans clusters 5 randompoint kseed 42 \
  rms :2-61@P sieve 20 random sieveseed 23 savepairdist \
  summaryhalf \$PREFIX.split.dat splitframe $SPLIT1,$SPLIT2,$SPLIT3 \
  info \$PREFIX.info.dat \
  cpopvtime \$PREFIX.cpop.agr normframe
EOF
  #summary \$PREFIX.summary.dat \
  $CPPTRAJ -i cpptraj.in -o cpptraj.$TOTAL.out
done
