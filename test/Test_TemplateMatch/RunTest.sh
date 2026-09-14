#!/bin/bash

. ../MasterTest.sh

CleanFiles templatematch.in \
           fle_to_ern.map ern_naorder.map \
           phenol_to_benzene.map sec_to_cys.map \
           daa_identity.map

TESTNAME='Template match tests'
Requires maxthreads 1

INPUT='-i templatematch.in'

# Nucleic acids: unmodified rA (FLE) onto unmodified dA (ERN).
# Shared atoms keep ERN indices; 2'-OH is an insertion; ERN H2'' is unmatched.
UNITNAME='Nucleic acid: FLE (rA) onto ERN (dA)'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > templatematch.in <<EOF
readdata FLE.lib name FLE
readdata ERN.lib name ERN
templatematch FLE[FLE] template ERN[ERN] mapout fle_to_ern.map maponly
EOF
  RunCpptraj "$UNITNAME"
  DoTest fle_to_ern.map.save fle_to_ern.map
fi

# Freeze a nucleotide template in the canonical sugar-phosphate walk.
UNITNAME='Nucleic acid: ERN naorder template walk'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > templatematch.in <<EOF
readdata ERN.lib name ERN
templatematch ERN[ERN] naorder mapout ern_naorder.map maponly
EOF
  RunCpptraj "$UNITNAME"
  DoTest ern_naorder.map.save ern_naorder.map
fi

# Small molecule: phenol onto benzene. OH replaces H1; ring is shared.
UNITNAME='Small molecule: phenol onto benzene'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > templatematch.in <<EOF
parm benzene.mol2 name benzene
parm phenol.mol2 name phenol
templatematch phenol template benzene mapout phenol_to_benzene.map maponly
EOF
  RunCpptraj "$UNITNAME"
  DoTest phenol_to_benzene.map.save phenol_to_benzene.map
fi

# Amino acid: selenocysteine onto cysteine.
# Must be Amber OFF so SE is selenium (elmnt 34). Mol2 atom name SE is sulfur.
UNITNAME='Amino acid: selenocysteine onto cysteine'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > templatematch.in <<EOF
readdata cys.lib name CYS
readdata sec.lib name SEC
templatematch SEC[SEC] template CYS[CYS] mapout sec_to_cys.map maponly
EOF
  RunCpptraj "$UNITNAME"
  DoTest sec_to_cys.map.save sec_to_cys.map
fi

# Smoke test: official ModXNA parent shipped in dat/templatematch/.
UNITNAME='Shipped ModXNA parent: DAA identity map'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > templatematch.in <<EOF
parm ../../dat/templatematch/DAA.mol2 name DAA
templatematch DAA template DAA mapout daa_identity.map maponly
EOF
  RunCpptraj "$UNITNAME"
  DoTest daa_identity.map.save daa_identity.map
fi

EndTest
exit 0
