#!/bin/bash

. ../MasterTest.sh

CleanFiles for.in TRP.vec.dat TRP.rms.dat TRP.CA.dist.dat TRP.tocenter.dat \
           nh.dat rms.nofit.dat last10.dat distance.dat nested.agr \
           EndToEnd0.dat EndToEnd1.dat EndToEnd2.agr temp.*.dat \
           DataOut.dat testset.dist.dat

TESTNAME='Loop tests'
Requires netcdf maxthreads 10

INPUT='-i for.in'
cat > for.in <<EOF
set FNAME= ../tz2.nc
parm ../tz2.parm7
trajin \$FNAME 1 10
for residues T inmask :TRP nvec=1;nvec++
  vector v\$nvec \$T center out TRP.vec.dat
  rms RMS\$T first \$T out TRP.rms.dat
  for atoms A0 inmask @CA
    distance d\$T\$A0 \$T \$A0 out TRP.CA.dist.dat
  done
  distance dcenter\$T \$T :1-12 out TRP.tocenter.dat
done

for atoms A0 inmask :2-4@N atoms A1 inmask :2-4@H
  distance d\$A0\$A1 \$A0 \$A1 out nh.dat
done

rms :1-12&!@H=
for residues MyRes inmask :1-12 r=1;r++
  rms R\$r \$MyRes&!@H= nofit out rms.nofit.dat
done
show

# Print info for the last 10 atoms. This tests using data set values
# as script variables and replacement of multiple script variables
# in a single argument.
set Natom = atoms inmask *
last10 = \$Natom - 10
show
atoms "@\$last10 - \$Natom" out last10.dat
run

EOF
RunCpptraj "$TESTNAME"
DoTest TRP.vec.dat.save TRP.vec.dat
DoTest TRP.rms.dat.save TRP.rms.dat
DoTest TRP.CA.dist.dat.save TRP.CA.dist.dat
DoTest TRP.tocenter.dat.save TRP.tocenter.dat
DoTest nh.dat.save nh.dat
DoTest rms.nofit.dat.save rms.nofit.dat
DoTest last10.dat.save last10.dat

UNITNAME='Test set command'
CheckFor maxthreads 10
if [ $? -eq 0 ] ; then
  cat > for.in <<EOF
parm ../FtuFabI.NAD.TCL.parm7
trajin ../FtuFabI.NAD.TCL.nc

set NATOM = atoms inmask ^1
set NRES = residues inmask ^1
set NMOL = molecules inmask *
show

distance @1 @\$NATOM out testset.dist.dat
distance :1 :\$NRES  out testset.dist.dat
distance ^1 ^\$NMOL  out testset.dist.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest testset.dist.dat.save testset.dist.dat
fi

UNITNAME='Test reading comma-separated list of strings'
CheckFor maxthreads 4
if [ $? -eq 0 ] ; then
  cat > for.in <<EOF
parm ../tz2.parm7
# Test reading in a comma-separated list of strings
clear trajin
for TRAJ in ../tz2.nc,../tz2.crd,../tz2.pdb,../tz2.rst7
  trajin \$TRAJ lastframe
done
distance D1-12 :1 :12 out distance.dat
run
EOF
  RunCpptraj "$UNITNAME"
  DoTest distance.dat.save distance.dat
fi

UNITNAME='Test nested loops and variables in loops'
cat > for.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc

set maxval = 4
set maxval1 = 3
set val = 0
for i=1;i<\$maxval;i++ myval=0;myval+=10
  startj = \$i + 1
  for j=\$startj;j<=\$maxval1;j++
    distance d\$i\$j @\$i @\$j out nested.agr
  done
done
run
show
list
EOF
RunCpptraj "$UNITNAME"
DoTest nested.agr nested.agr.save

UNITNAME='Test loop over data set blocks'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > for.in <<EOF
parm ../tz2.parm7

starttraj = 1
endtraj = 10
offset = 10
for i=1;i<=10;i++
  trajin ../tz2.nc \$starttraj \$endtraj
  distance EndToEnd\$endtraj :1 :12
  run
  starttraj = \$starttraj + \$offset
  endtraj = \$endtraj + \$offset
  clear trajin
done
writedata EndToEnd0.dat EndToEnd*

trajin ../tz2.nc 1 100
distance EndToEndAll :1 :12
run

for DS datasetblocks EndToEndAll blocksize 10 i=1;i++
  set WriteVars1 += " \$DS "
  show DS WriteVars1
done
writedata EndToEnd1.dat \$WriteVars1 

for ES datasetblocks EndToEndAll blocksize 10 firstblock 50 cumulative j=\$i;j++
  set WriteVars2 += " \$ES "
  show ES WriteVars2
done
writedata EndToEnd2.agr \$WriteVars2
list
EOF
  RunCpptraj "$UNITNAME"
  DoTest EndToEnd1.dat.save EndToEnd1.dat
  DoTest EndToEnd2.agr.save EndToEnd2.agr
fi

UNITNAME='Test for filename wildcards and oversets'
CheckFor notos windows
if [ $? -eq 0 ] ; then
  cat > for.in <<EOF
for FILE in distance.dat.save,last10.dat.save,doesnot*,TRP.*.dat.save
  readdata \$FILE
done

for SET1 oversets TRP.rms.dat.save* SET2 oversets TRP.tocenter.dat.save* i=1;i++
  show SET1 SET2
  NewSet\$i = \$SET1 + \$SET2
  set DataOut += " NewSet\$i "
done
writedata DataOut.dat \$DataOut
list datasets
EOF
  RunCpptraj "$UNITNAME"
  DoTest DataOut.dat.save DataOut.dat
fi

EndTest
exit 0
