#!/bin/bash

. ../UnitMaster.sh

CleanFiles Makefile main.o a.out CpptrajStdio.o GistEntropyUtils.o

UNITSOURCES='GistEntropyUtils.cpp'

CreateMakefile

RunMake "GIST Entropy Utilities unit tests."

EndTest
