#!/bin/bash

. ../MasterTest.sh

CleanFiles templatematch.in \
           fle_to_ern.map ern_naorder.map \
           phenol_to_benzene.map sec_to_cys.map \
           daa_identity.map \
           fle_ern.0.mol2 fle_ern.1.mol2 fle_ern.0.lib fle_ern.1.lib \
           fle_ern.scmask fle_ern.atoms fle_ern.map \
           phenol_bnz.0.mol2 phenol_bnz.1.mol2 phenol_bnz.0.lib phenol_bnz.1.lib \
           phenol_bnz.scmask phenol_bnz.atoms phenol_bnz.map \
           sec_cys.0.mol2 sec_cys.1.mol2 sec_cys.0.lib sec_cys.1.lib \
           sec_cys.scmask sec_cys.atoms sec_cys.map

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

# Dual-topology TI export: dummies, matching atom counts, mol2+lib+scmask.
# FLE (rA) onto ERN (dA): 31 shared + dummy O2'/HO2' on lambda 0 + dummy H2'' on lambda 1 = 34.
UNITNAME='TI export: FLE onto ERN (dummies, mol2, lib, scmask)'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > templatematch.in <<EOF
readdata FLE.lib name FLE
readdata ERN.lib name ERN
templatematch FLE[FLE] template ERN[ERN] tiout fle_ern maponly
EOF
  RunCpptraj "$UNITNAME"
  DoTest fle_ern.atoms.save fle_ern.atoms
  DoTest fle_ern.scmask.save fle_ern.scmask
fi

UNITNAME='TI export: phenol onto benzene'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > templatematch.in <<EOF
parm benzene.mol2 name benzene
parm phenol.mol2 name phenol
templatematch phenol template benzene tiout phenol_bnz maponly
EOF
  RunCpptraj "$UNITNAME"
  DoTest phenol_bnz.atoms.save phenol_bnz.atoms
  DoTest phenol_bnz.scmask.save phenol_bnz.scmask
fi

UNITNAME='TI export: selenocysteine onto cysteine'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > templatematch.in <<EOF
readdata cys.lib name CYS
readdata sec.lib name SEC
templatematch SEC[SEC] template CYS[CYS] tiout sec_cys maponly
EOF
  RunCpptraj "$UNITNAME"
  DoTest sec_cys.atoms.save sec_cys.atoms
  DoTest sec_cys.scmask.save sec_cys.scmask
fi

EndTest
exit 0
