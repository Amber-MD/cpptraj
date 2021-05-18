#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o NameType.o a.out CpptrajStdio.o

UNITSOURCES='NameType.cpp CpptrajStdio.cpp'

CreateMakefile

RunMake "NameType class unit test."

EndTest
