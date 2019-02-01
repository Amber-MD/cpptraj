#!/bin/bash

. ../MasterTest.sh

CleanFiles matrix.dat hd.in hd.dat rms2d.dat rms2d.gnu hausdorff.matrix.dat \
           hausdoff.matrix.gnu hausdorff.fullmatrix.gnu

INPUT='-i hd.in'

# Simple Hausdorff distance from matrix.
cat > matrix.dat <<EOF
1.0000       2.2361       3.0000       4.1231
2.2361       1.0000       2.2361       5.0000
3.0000       2.2361       1.0000       4.1231
2.2361       3.0000       2.2361       3.0000
EOF

cat > hd.in <<EOF
readdata matrix.dat read2d name Matrix
runanalysis hausdorff Matrix out hd.dat name HD outab hd.dat outba hd.dat
EOF
RunCpptraj "Simple Hausdorff distance test."
DoTest hd.dat.save hd.dat

# 2D RMS
#cat > hd.in <<EOF
#parm ../DPDP.parm7
#trajin ../DPDP.nc
#rms2d DPDP out rms2d.gnu
#EOF
#RunCpptraj "2D RMS"

# Create 10 traectory chunks, do Hausdorff between 2D rms sets
cat > hd.in <<EOF
parm ../DPDP.parm7
for beg=1;beg<100;beg+=10 end=10;end+=10 i=1;i++
  loadcrd ../DPDP.nc \$beg \$end name Chunk\$i
done
# Do the 2drms in chunks
EOF
for ((i=1; i < 11; i++)) ; do
  ((start = $i + 1))
  for ((j=$start; j < 11; j++)) ; do
#for i=1;i<11;i++
#  for j=1;j<11;j++
  echo "2drms crdset Chunk$i reftraj Chunk$j M$i.$j" >> hd.in
  done
done
cat >> hd.in <<EOF
hausdorff M* out hausdorff.matrix.gnu outtype trimatrix nrows 10
runanalysis
list dataset
quit
EOF
RunCpptraj "Hausdorff distance of 2D rms output test."
DoTest hausdorff.matrix.gnu.save hausdorff.matrix.gnu

# Create 10 traectory chunks, do Hausdorff between 2D rms sets, full matrix
cat > hd.in <<EOF
parm ../DPDP.parm7
for beg=1;beg<100;beg+=10 end=10;end+=10 i=1;i++
  loadcrd ../DPDP.nc \$beg \$end name Chunk\$i
done
# Do the 2drms in chunks
for i=1;i<11;i++
  for j=1;j<11;j++
    2drms crdset Chunk\$i reftraj Chunk\$j M\$i.\$j
  done
done
hausdorff M* out hausdorff.fullmatrix.gnu title hausdorff.matrix.gnu outtype fullmatrix nrows 10
runanalysis
list dataset
quit
EOF
RunCpptraj "Hausdorff distance of 2D rms output test."
DoTest hausdorff.matrix.gnu.save hausdorff.fullmatrix.gnu

EndTest
exit 0
