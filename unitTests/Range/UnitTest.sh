#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o Range.o StringRoutines.o a.out CpptrajStdio.o ArgList.o

UNITSOURCES='Range.cpp StringRoutines.cpp CpptrajStdio.cpp ArgList.cpp'

CreateMakefile

RunMake "Range unit test."

EndTest
