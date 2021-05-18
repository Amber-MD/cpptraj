#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o ArgList.o a.out CpptrajStdio.o StringRoutines.o

UNITSOURCES='ArgList.cpp CpptrajStdio.cpp StringRoutines.cpp'

CreateMakefile

RunMake "ArgList class unit test."

EndTest
