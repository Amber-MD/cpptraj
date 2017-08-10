#!/bin/bash

. ../MasterTest.sh

CleanFiles data.in avg.dat matrix.dat matrix2.dat VXYZ.dat Keep?.dat Drop?.dat

INPUT='-i data.in'
cat > data.in <<EOF
readdata data.dat name MyData1
readdata data.dat name MyData2
readdata data1.dat name D1
readdata data2.dat name D2
readdata data3.dat name D3
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

# Test removal of points
dataset keeppoints D3 range 3,5,8-10 name Keep1
writedata Keep1.dat Keep1
dataset keeppoints D3 start 1 stop 6 offset 2 name Keep2
writedata Keep2.dat Keep2
dataset droppoints D3 range 1-2,4,6-7 name Drop1
writedata Drop1.dat Drop1
dataset droppoints D3 start 2 stop 10 offset 2
dataset droppoints D3 range 4,5
writedata Drop2.dat D3

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
DoTest Keep1.dat.save Keep1.dat
DoTest Keep2.dat.save Keep2.dat
DoTest Keep1.dat.save Drop1.dat
DoTest Keep2.dat.save Drop2.dat
DoTest VXYZ.dat.save VXYZ.dat

EndTest
