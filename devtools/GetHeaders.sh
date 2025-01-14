#!/bin/bash

# Get all header files

WORKDIR=`pwd`
BASEDIR=`basename $WORKDIR`
if [ "$BASEDIR" != 'src' ] ; then
  echo "Should be executed from the CPPTRAJ src directory."
  exit 1
fi

OUTFILE='cpptrajheaders'
echo "# All cpptraj headers that should be installed for the cpptraj library." > $OUTFILE
echo "CPPTRAJ_HEADERS = \\" >> $OUTFILE

for DIR in . Cluster Structure Energy ; do
  for FILE in `ls $DIR/*.h` ; do
    echo "        $FILE \\" >> $OUTFILE
  done
done
echo "" >> $OUTFILE

exit 0
