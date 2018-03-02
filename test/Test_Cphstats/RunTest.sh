#!/bin/bash

. ../MasterTest.sh

CleanFiles cphstats.in sorted.pH_*.00 stats.dat frac.agr implicit.sorted.dat \
           explicit.004.cpout implicit.001.cpout smallimplicit.sorted.?.cpout \
           smallimplicit.stats.dat

INPUT='-i cphstats.in'
TESTNAME='Constant pH stats / data set sort test'

UNITNAME='Explicit pH REMD ensemble data read / sort'
SKIP='no'
if [ ! -z $N_THREADS ] ; then
  if [ $N_THREADS -lt 7 ] ; then
    if [ $N_THREADS -eq 4 -o $N_THREADS -eq 5 ] ; then
      echo "  $UNITNAME cannot be run with 4 or 5 threads."
      ((CHECKERR++))
      SkipCheck "$UNITNAME"
      SKIP='yes'
    fi
  else
    CheckFor nthreads 6
    if [ $? -eq 1 ] ; then
      SKIP='yes'
    fi
  fi
fi
if [ "$SKIP" = 'no' ] ; then
  cat > cphstats.in <<EOF
readensembledata ExplicitRemd/cpout.001 cpin ExplicitRemd/cpin name PH
#readensembledata cpout.001 filenames cpout.002,cpout.003,cpout.004,cpout.005,cpout.006 name PH
list dataset
sortensembledata PH
list dataset
for i=0;i<6;i++ j=1;j++
  writedata sorted.pH_\$j.00 cpout PH[*]%\$i noensextension
done
EOF
  RunCpptraj "$UNITNAME"
  DoTest sorted.pH_1.00.save sorted.pH_1.00
  DoTest sorted.pH_2.00.save sorted.pH_2.00
  DoTest sorted.pH_3.00.save sorted.pH_3.00
  DoTest sorted.pH_4.00.save sorted.pH_4.00
  DoTest sorted.pH_5.00.save sorted.pH_5.00
  DoTest sorted.pH_6.00.save sorted.pH_6.00
fi

UNITNAME='Constant pH stats test'
cat > cphstats.in <<EOF
#readensembledata sorted.pH_1.00.save filenames sorted.pH_2.00.save,sorted.pH_3.00.save,sorted.pH_4.00.save,sorted.pH_5.00.save,sorted.pH_6.00.save cpin cpin name PH
set CPIN = ExplicitRemd/cpin
readdata sorted.pH_1.00.save separate cpin \$CPIN name PH1
readdata sorted.pH_2.00.save separate cpin \$CPIN name PH2
readdata sorted.pH_3.00.save separate cpin \$CPIN name PH3
readdata sorted.pH_4.00.save separate cpin \$CPIN name PH4
readdata sorted.pH_5.00.save separate cpin \$CPIN name PH5
readdata sorted.pH_6.00.save separate cpin \$CPIN name PH6
list datasets

ensextension off
runanalysis cphstats PH*[*] statsout stats.dat fracplot fracplotout frac.agr deprot
list dataset
EOF
RunCpptraj "$UNITNAME"
DoTest stats.dat.save stats.dat
DoTest frac.agr.save frac.agr

UNITNAME='Sorted implicit constant pH stats test'
cat > cphstats.in <<EOF
readdata ImplicitRemd/md2_cpout.pH_2.00.save cpin ImplicitRemd/1AKI.dry.equil.cpin name PH
runanalysis cphstats PH[*] statsout implicit.sorted.dat
EOF
RunCpptraj "$UNITNAME"
DoTest implicit.sorted.dat.save implicit.sorted.dat

UNITNAME='Implicit pH REMD ensemble data read / sort / stats'
SKIP='no'
if [ ! -z $N_THREADS ] ; then
  if [ $N_THREADS -lt 5 ] ; then
    if [ $N_THREADS -eq 3 ] ; then
      echo "  $UNITNAME cannot be run with 3 threads."
      ((CHECKERR++))
      SkipCheck "$UNITNAME"
      SKIP='yes'
    fi
  else
    CheckFor nthreads 4
    if [ $? -eq 1 ] ; then
      SKIP='yes'
    fi
  fi
fi
if [ "$SKIP" = 'no' ] ; then
  cat > cphstats.in <<EOF
#readensembledata SmallImplicitRemd/run0.cpout.00* cpin SmallImplicitRemd/cpin name PH
# FIXME wildcard matching does not work on windows
readensembledata SmallImplicitRemd/run0.cpout.000 filenames SmallImplicitRemd/run0.cpout.001,SmallImplicitRemd/run0.cpout.002,SmallImplicitRemd/run0.cpout.003 cpin SmallImplicitRemd/cpin name PH

sortensembledata PH

writedata smallimplicit.sorted.0.cpout noensextension PH[*]%0
#writedata temp.dat noensextension PH[*]%0 xmin 0 xstep 1
writedata smallimplicit.sorted.1.cpout noensextension PH[*]%1
writedata smallimplicit.sorted.2.cpout noensextension PH[*]%2
writedata smallimplicit.sorted.3.cpout noensextension PH[*]%3

ensextension off
cphstats PH[*] statsout smallimplicit.stats.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest smallimplicit.sorted.0.cpout.save smallimplicit.sorted.0.cpout
  DoTest smallimplicit.sorted.1.cpout.save smallimplicit.sorted.1.cpout
  DoTest smallimplicit.sorted.2.cpout.save smallimplicit.sorted.2.cpout
  DoTest smallimplicit.sorted.3.cpout.save smallimplicit.sorted.3.cpout
  DoTest smallimplicit.stats.dat.save smallimplicit.stats.dat
fi

UNITNAME='Unsorted pH read/write test'
cat > cphstats.in <<EOF
readdata ExplicitRemd/cpout.004 cpin ExplicitRemd/cpin name PH4
writedata explicit.004.cpout PH4
readdata SmallImplicitRemd/run0.cpout.001 cpin SmallImplicitRemd/cpin name PH
writedata implicit.001.cpout PH
EOF
RunCpptraj "$UNITNAME"
DoTest ExplicitRemd/cpout.004 explicit.004.cpout
DoTest SmallImplicitRemd/run0.cpout.001 implicit.001.cpout

EndTest
exit 0
