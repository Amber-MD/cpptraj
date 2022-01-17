#!/bin/bash


rm time.pw.dat time.cluster.dat time.restore.dat time.total.dat time.summary.dat \
   time.calcoutput.dat
source ../framelist.sh
for TOTAL in $FRAMELIST ; do
  OUTFILE=cpptraj.$TOTAL.out
  grep TIME $OUTFILE | awk -v nframes=$TOTAL '{
    if ($2 == "Pairwise") printf("%i %f\n", nframes, $5) >> "time.pw.dat";
    if ($2 == "Clustering") printf("%i %f\n", nframes, $4) >> "time.cluster.dat";
    if ($4 == "restore") printf("%i %f\n", nframes, $6) >> "time.restore.dat";
    if ($2 == "Analyses") printf("%i %f\n", nframes, $4) >> "time.total.dat";
    if ($2 == "Summary") printf("%i %f\n", nframes, $5) >> "time.summary.dat";
    if ($2 == "Output") printf("%i %f\n", nframes, $5) >> "time.calcoutput.dat";
  }'
done
