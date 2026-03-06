#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o Frame.o a.out CpptrajStdio.o

UNITSOURCES='Frame.cpp CpptrajStdio.cpp Box.cpp CoordinateInfo.cpp Vec3.cpp Matrix_3x3.cpp'

CreateMakefile

RunMake "Frame class unit test."

EndTest
