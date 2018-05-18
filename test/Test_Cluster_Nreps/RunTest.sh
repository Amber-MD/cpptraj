#!/bin/bash

. ../MasterTest.sh

CleanFiles cluster.in 2drms.gnu clusters.*.dat summary.*.dat singlerep.nc \
           Rep.c*.nc PD

INPUT='-i cluster.in'

cat > cluster.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
#2drms !@H= out 2drms.gnu
cluster kmeans clusters 2 info clusters.2.dat summary summary.2.dat \
        savepairdist pairdist PD savenreps 5 

#\
#        singlerepout singlerep.nc \
#        repout Rep repfmt netcdf repframe
EOF
RunCpptraj "Clustering with specified # reps"

DoTest clusters.2.dat.save clusters.2.dat
DoTest summary.2.dat.save summary.2.dat

EndTest
exit 0
