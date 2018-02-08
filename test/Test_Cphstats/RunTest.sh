#!/bin/bash

. ../MasterTest.sh

CleanFiles cphstats.in sorted.pH_*.00 stats.dat frac.agr

INPUT='-i cphstats.in'
TESTNAME='Constant pH stats / data set sort test'

UNITNAME='Ensemble data read / sort'
cat > cphstats.in <<EOF
readensembledata cpout.001 cpin cpin name PH
#readensembledata cpout.001 filenames cpout.002,cpout.003,cpout.004,cpout.005,cpout.006 name PH
list dataset
sortensembledata PH
list dataset
for i=0;i<6;i++ j=1;j++
  writedata sorted.pH_\$j.00 as cpout PH[*]%\$i nwriteheader 50
done
EOF
RunCpptraj "$UNITNAME"
DoTest sorted.pH_1.00.save sorted.pH_1.00
DoTest sorted.pH_2.00.save sorted.pH_2.00
DoTest sorted.pH_3.00.save sorted.pH_3.00
DoTest sorted.pH_4.00.save sorted.pH_4.00
DoTest sorted.pH_5.00.save sorted.pH_5.00
DoTest sorted.pH_6.00.save sorted.pH_6.00

UNITNAME='Constant pH stats test'
cat > cphstats.in <<EOF
#readensembledata sorted.pH_1.00.save filenames sorted.pH_2.00.save,sorted.pH_3.00.save,sorted.pH_4.00.save,sorted.pH_5.00.save,sorted.pH_6.00.save cpin cpin name PH
readdata sorted.pH_*.00 separate cpin cpin name PH
list datasets
#runanalysis cphstats PH[*] statsout stats.dat fracplot fracplotout frac.agr deprot
runanalysis cphstats PH*[*] statsout stats.dat fracplot fracplotout frac.agr deprot
list dataset
EOF
RunCpptraj "$UNITNAME"
DoTest stats.dat.save stats.dat
DoTest frac.agr.save frac.agr

EndTest
exit 0
