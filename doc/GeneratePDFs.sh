#!/bin/bash

LYX=`which lyx`
if [ -z "$LYX" ] ; then
  echo "No lyx present."
  exit 1
fi
MD5SUM=`which md5sum`
if [ -z "$MD5SUM" ] ; then
  echo "No md5sum present."
  exit 1
fi

GEN_MANUAL=0
C0=`grep CpptrajManual.lyx DocumentChecksums.txt`
M0=`md5sum CpptrajManual.lyx`
echo $C0
echo $M0
if [ "$C0" != "$M0" ] ; then
  GEN_MANUAL=1
fi
C1=`grep cpptraj.lyx DocumentChecksums.txt`
M1=`md5sum cpptraj.lyx`
echo $C1
echo $M1
if [ "$C1" != "$M1" ] ; then
  GEN_MANUAL=1
fi
if [ $GEN_MANUAL -eq 1 ] ; then
  lyx -batch --export pdf2 CpptrajManual.lyx
  if [ $? -ne 0 ] ; then
    echo "Generation of manual failed."
    exit 1
  fi
  ls -l CpptrajManual.pdf
fi

GEN_DEV=0
C2=`grep CpptrajDevelopmentGuide.lyx DocumentChecksums.txt`
M2=`md5sum CpptrajDevelopmentGuide.lyx`
echo $C2
echo $M2
if [ "$C2" != "$M2" ] ; then
  GEN_DEV=1
fi
if [ $GEN_DEV -eq 1 ] ; then
  lyx -batch --export pdf2 CpptrajDevelopmentGuide.lyx
  if [ $? -ne 0 ] ; then
    echo "Generation of dev guide failed."
    exit 1
  fi
  ls -l CpptrajDevelopmentGuide.pdf
fi

# Update checksums if needed
if [ $GEN_MANUAL -ne 0 -o $GEN_DEV -ne 0 ] ; then
  echo "Updating document checksums."
  md5sum *.lyx > DocumentChecksums.txt
  ls -l DocumentChecksums.txt
fi
exit 0
