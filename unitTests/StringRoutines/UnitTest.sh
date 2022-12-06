#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o StringRoutines.o a.out CpptrajStdio.o

UNITSOURCES='StringRoutines.cpp CpptrajStdio.cpp'

CreateMakefile

RunMake "StringRoutines unit test."

EndTest
