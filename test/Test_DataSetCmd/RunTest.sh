#!/bin/bash

. ../MasterTest.sh

CleanFiles data.in avg.dat matrix.dat matrix2.dat VXYZ.dat

INPUT='-i data.in'
cat > data.in <<EOF
readdata data.dat name MyData1
readdata data.dat name MyData2
readdata data1.dat name D1
readdata data2.dat name D2
readdata ../Test_Vector/vtest.dat.6.save name Vec vector

# Test mode/type setting
dataset MyData2 mode torsion type alpha

# Test legend setting
dataset MyData2 legend Alpha

# Test concatenation
dataset cat D1 D2 name D12

# Test MakeXY
dataset makexy D1 D2 name Dxy

# Test Make2D
dataset make2d D12 name Dmat ncols 2 nrows 2

# Test extraction of vector coords
dataset vectorcoord Vec name VX X
dataset vectorcoord Vec name VY Y
dataset vectorcoord Vec name VZ Z
writedata VXYZ.dat VX VY VZ

list dataset
runanalysis avg MyData* D12 Dxy out avg.dat name MyAvg
writedata matrix.dat Dmat square2d

# Test remove
dataset remove ifaverage outside 0 and 4
dataset remove ifsize == 2
dataset remove ifmode notequal matrix
writedata matrix2.dat * square2d
list dataset
EOF
RunCpptraj "Data Set commands test"
DoTest avg.dat.save avg.dat
DoTest matrix.dat.save matrix.dat
DoTest matrix.dat.save matrix2.dat
DoTest VXYZ.dat.save VXYZ.dat

EndTest
