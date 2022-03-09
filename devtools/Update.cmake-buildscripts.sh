#!/bin/bash

# Adapt pieces of the Amber cmake build system needed by CPPTRAJ. 
A2C_COUNT=0

if [ -z "$AMBERHOME" ] ; then
  echo "AMBERHOME is not set."
  exit 1
fi
AMBERCMAKE="$AMBERHOME/cmake"

WORKDIR=`pwd`

# AmberToCpptraj <amber filename> <cpptraj filename>
# This is a 1 to 1 copy of the file, same name, maybe different location.
# If file is already copied, report any differences.
AmberToCpptraj() {
  ((A2C_COUNT++))
  amberfile=$AMBERCMAKE/$1
  if [ -z "$2" ] ; then
    cpptrajfile=$1
  else
    cpptrajfile=$2
  fi
  #echo "AmberToCpptraj: $amberfile $cpptrajfile"
  if [ ! -f "$amberfile" ] ; then
    echo "  ERROR: Amber file $amberfile not present"
    exit 1
  elif [ ! -f "$cpptrajfile" ] ; then
    echo "  COPY: $amberfile to $cpptrajfile"
    cp -i $amberfile $cpptrajfile
    if [ $? -ne 0 ] ; then
      echo "Copy failed."
      exit 1
    fi
  else
    diff $amberfile $cpptrajfile > temp.diff
    if [ -s 'temp.diff' ] ; then
      ndiff=`cat temp.diff | wc -l`
      echo "  CHECK: $amberfile $cpptrajfile diff: $ndiff"
    fi
    rm temp.diff
  fi
}

# CompareToAmber <amber filename> <cpptraj filename>
# Assuming the cpptraj module is altered in necessary ways, just report the
# differences from the Amber module.
CompareToAmber() {
  ((A2C_COUNT++))
  amberfile=$AMBERCMAKE/$1
  if [ -z "$2" ] ; then
    cpptrajfile=$1
  else
    cpptrajfile=$2
  fi
  echo "CompareToAmber: $amberfile $cpptrajfile"
  if [ ! -f "$amberfile" ] ; then
    echo "  ERROR: Amber file $amberfile not present"
    exit 1
  elif [ ! -f "$cpptrajfile" ] ; then
    echo "  ERROR: Cpptraj file $cpptrajfile not present"
  else
    diff $amberfile $cpptrajfile > temp.diff
    if [ -s 'temp.diff' ] ; then
      ndiff=`cat temp.diff | wc -l`
      echo "  DIFF: $amberfile $cpptrajfile diff: $ndiff"
    else
      echo "  WARNING: no differences detected."
    fi
    rm temp.diff
  fi
}

# CpptrajOnly <cpptraj filename>
# Modules that should only exist here
CpptrajOnly() {
  ((A2C_COUNT++))
  amberfile=$AMBERCMAKE/$1
  cpptrajfile=$1
  echo "CpptrajOnly: $cpptrajfile"
  if [ -f "$amberfile" ] ; then
    echo "  ERROR: $amberfile exists in Amber cmake."
    exit 1
  elif [ ! -f "$cpptrajfile" ] ; then
    echo "  ERROR: Cpptraj file $cpptrajfile not present"
  fi
}


# ==============================================================================
# 1 to 1 files
AmberToCpptraj AmberCompilerConfig.cmake
AmberToCpptraj BuildReport.cmake
AmberToCpptraj CheckLinkerFlag.cmake
AmberToCpptraj ColorMessage.cmake
AmberToCpptraj CompilationOptions.cmake
AmberToCpptraj CompilerFlags.cmake
AmberToCpptraj CopyTarget.cmake
AmberToCpptraj CudaConfig.cmake
AmberToCpptraj patched-cmake-modules/FindOpenMPFixed.cmake FindOpenMPFixed.cmake
AmberToCpptraj InstallDirs.cmake
AmberToCpptraj LibraryBundling.cmake
AmberToCpptraj LibraryTracking.cmake
AmberToCpptraj LibraryUtils.cmake
AmberToCpptraj MPIConfig.cmake
AmberToCpptraj NetlibConfig.cmake
AmberToCpptraj OpenMPConfig.cmake
AmberToCpptraj PackageTypes.cmake
AmberToCpptraj Packaging.cmake
AmberToCpptraj ParallelizationConfig.cmake
AmberToCpptraj Policies.cmake
AmberToCpptraj Shorthand.cmake
AmberToCpptraj TargetArch.cmake
AmberToCpptraj TryLinkLibrary.cmake
AmberToCpptraj Utils.cmake
# ThirdPartyTools
AmberToCpptraj FindARPACK.cmake ThirdPartyTools/FindARPACK.cmake
AmberToCpptraj patched-cmake-modules/FindBLASFixed.cmake ThirdPartyTools/FindBLASFixed.cmake
AmberToCpptraj FindCMath.cmake ThirdPartyTools/FindCMath.cmake
AmberToCpptraj jedbrown/FindFFTW.cmake ThirdPartyTools/FindFFTW.cmake
AmberToCpptraj patched-cmake-modules/FindLAPACKFixed.cmake ThirdPartyTools/FindLAPACKFixed.cmake
AmberToCpptraj hanjianwei/FindMKL.cmake ThirdPartyTools/FindMKL.cmake
AmberToCpptraj jedbrown/FindNetCDF.cmake ThirdPartyTools/FindNetCDF.cmake
AmberToCpptraj FindPnetCDF.cmake ThirdPartyTools/FindPnetCDF.cmake
AmberToCpptraj FindReadline.cmake ThirdPartyTools/FindReadline.cmake

echo "$A2C_COUNT cmake modules direct from Amber"

# These are based on module in Amber but have necessary modifications
CompareToAmber 3rdPartyTools.cmake SetupThirdParty.cmake
CompareToAmber AmberBuildSystemInit.cmake BuildSystemInit.cmake
CompareToAmber AmberBuildSystem2ndInit.cmake BuildSystemSetup.cmake

# These should only exist in the cpptraj cmake system
CpptrajOnly DebugCpptrajCmake.cmake

echo "$A2C_COUNT total cmake modules."
echo "`ls *.cmake */*.cmake | wc -l` cmake modules present."
exit 0
