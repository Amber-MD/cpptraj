#!/bin/bash

# Test of IRED 
. ../MasterTest.sh

CleanFiles orderparam ired.vec noe v0.cjt v0.cmt ptraj_ired_step1.in

CheckPtrajAnalyze

TOP="1IEE_A_prot.prmtop"
INPUT="ptraj_ired_step1.in"
echo "trajin 1IEE_A_test.mdcrd" > $INPUT
N=0
# Define N-H vectors. H atom is always N atom plus one.
for NATOM in 25   41   61   68   92   102  117  136  146  156  166  183  205  \
             229  246  253  272  284  298  319  343  350  371  382  401  408  \
             422  446  462  472  482  492  514  534  549  560  574  594  608  \
             622  639  649  663  677  701  715  729  741  748  759  773  785  \
             806  813  832  851  868  887  901  912  936  960  984  994  1008 \
             1020 1027 1051 1079 1086 1097 1121 1135 1154 1164 1178 1211 1221 \
             1232 1242 1261 1280 1291 1302 1314 1333 1347 1357 1368 1384 1398 \
             1408 1418 1440 1462 1481 1497 1508 1520 1527 1541 1548 1565 1579 \
             1589 1613 1629 1639 1663 1687 1701 1725 1735 1757 1764 1778 1790 \
             1806 1823 1833 1857 1876 1900 1907 1917 1941 ; do
  ((HATOM = NATOM + 1))
  echo "vector v$N @$NATOM ired @$HATOM" >> $INPUT
  ((N++))
done
cat >> $INPUT <<EOF
matrix ired name matired order 2
###Attention: adjust vecs to number of residues in your protein
analyze matrix matired out ired.vec name ired.vec vecs 126 
analyze ired relax freq 500.0 NHdist 1.02 tstep 1.0 tcorr 10000.0 norm out v0 noefile noe \
        order 2 modes ired.vec orderparamfile orderparam
EOF

RunCpptraj "IRED vector/matrix test"
DoTest ired.vec.save ired.vec
DoTest orderparam.save orderparam
DoTest noe.save noe
#DoTest v0.cjt.save v0.cjt
#DoTest v0.cmt.save v0.cmt
CheckTest
EndTest

exit 0
