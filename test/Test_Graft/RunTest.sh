#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in Final.graft.mol2 Nucleotide.pdb Nucleotide.charge.mol2 \
           IC.Final.graft.mol2 IC.Nucleotide.pdb
TESTNAME='Graft test'
Requires notparallel

INPUT="-i cpptraj.in"

# Combine Tyr FF14SB backbone + CB with PRY fragment
TyrPry() {
  cat > cpptraj.in <<EOF
set TYRFILE = ../Test_CombineCrd/Tyr.mol2
# Load the target
parm \$TYRFILE 
loadcrd \$TYRFILE parmindex 0 name TYR

# Load the source
set PRYFILE= ../Test_CombineCrd/PRY-gauss-fragment.mol2
parm \$PRYFILE
loadcrd \$PRYFILE parmindex 1 name PRY

graft \
  tgt TYR \
  src PRY \
  name Final \
  srcfitmask @O1,C5,C6,C4,H6,H5,C7,C3,H7,H4,C2 \
  tgtfitmask @OH,CZ,CE1,CE2,HE1,HE2,CD1,CD2,HD1,HD2,CG \
  srcmask !(@1-4) \
  tgtmask @C,O,CA,HA,N,H,CB,HB2,HB3 \
  bond :1@CB,:1@C2
crdout Final Final.graft.mol2
zmatrix Final out zFinal.dat
EOF
  RunCpptraj "$TESTNAME, Prenylated Tyrosine"
  DoTest Final.graft.mol2.save Final.graft.mol2
}

# Combine Tyr FF14SB backbone + CB with PRY fragment, use internal coords
TyrPryIC() {
  cat > cpptraj.in <<EOF
set TYRFILE = ../Test_CombineCrd/Tyr.mol2
# Load the target
parm \$TYRFILE 
loadcrd \$TYRFILE parmindex 0 name TYR

# Load the source
set PRYFILE= ../Test_CombineCrd/PRY-gauss-fragment.mol2
parm \$PRYFILE
loadcrd \$PRYFILE parmindex 1 name PRY

graft ic \
  tgt TYR \
  src PRY \
  name Final \
  srcmask !(@1-4) \
  tgtmask @C,O,CA,HA,N,H,CB,HB2,HB3 \
  bond :1@CB,:1@C2
crdout Final IC.Final.graft.mol2
zmatrix Final out zicFinal.dat
EOF
  RunCpptraj "$TESTNAME, Prenylated Tyrosine (internal coords)"
  DoTest IC.Final.graft.mol2.save IC.Final.graft.mol2
}

# Link nucleic acid base + sugar + phosphate
DNA() {
  cat > cpptraj.in <<EOF
for FILE in template-sugar.pdb,template-dimethylphosphate.pdb,template-base-adenine.pdb NAME in Sugar,Phos,Base
  parm \$FILE
  loadcrd \$FILE name \$NAME parm \$FILE
done
list
# Base + Sugar
graft \\
  tgt Base \\
  src Sugar \\
  name BaseSugar \\
  tgtfitmask @N9,C1',1H1' \\
  srcfitmask @C10,C1',O4' \\
  tgtmask !(@C1',1H1',2H1',3H1') \\
  srcmask !(@C10,1H10,2H10,3H10,HO3',O2',HO2') \\
  bond @N9,@C1'
#crdout BaseSugar BaseSugar.mol2
# BaseSugar + Phosphate
graft \\
  tgt BaseSugar \\
  src Phos \\
  name Nucleotide \\
  tgtfitmask @C5',O5',HO5' \\
  srcfitmask @C2,O5',P \\
  tgtmask !(@O5',HO5') \\
  srcmask !(@C2,H21,H22,H23,O3',C1,H11,H12,H13) \\
  bond @C5',@O5'
#crdout Nucleotide Nucleotide.mol2
# Format the PDB
change crdset Nucleotide resname from * to DA
change crdset Nucleotide oresnums of :1 min 1 max 1
change crdset Nucleotide oresnums of :2 min 1 max 1
change crdset Nucleotide oresnums of :3 min 1 max 1
crdout Nucleotide Nucleotide.pdb
EOF
  RunCpptraj "$TESTNAME, Construct Nucleic Acid"
  DoTest Nucleotide.pdb.save Nucleotide.pdb
}

# Link nucleic acid base + sugar + phosphate, IC
DNAic() {
  cat > cpptraj.in <<EOF
for FILE in template-sugar.pdb,template-dimethylphosphate.pdb,template-base-adenine.pdb NAME in Sugar,Phos,Base
  parm \$FILE
  loadcrd \$FILE name \$NAME parm \$FILE
done
list
# Base + Sugar
graft ic \\
  tgt Base \\
  src Sugar \\
  name BaseSugar \\
  tgtmask !(@C1',1H1',2H1',3H1') \\
  srcmask !(@C10,1H10,2H10,3H10,HO3',O2',HO2') \\
  bond @N9,@C1'
zmatrix BaseSugar out zicBaseSugar.dat
#crdout BaseSugar BaseSugar.mol2
#quit
# BaseSugar + Phosphate
graft ic \\
  tgt BaseSugar \\
  src Phos \\
  name Nucleotide \\
  tgtmask !(@O5',HO5') \\
  srcmask !(@C2,H21,H22,H23,O3',C1,H11,H12,H13) \\
  bond @C5',@O5'
#crdout Nucleotide Nucleotide.ic.mol2
# Format the PDB
change crdset Nucleotide resname from * to DA
change crdset Nucleotide oresnums of :1 min 1 max 1
change crdset Nucleotide oresnums of :2 min 1 max 1
change crdset Nucleotide oresnums of :3 min 1 max 1
crdout Nucleotide IC.Nucleotide.pdb
EOF
  RunCpptraj "$TESTNAME, Construct Nucleic Acid, IC"
  DoTest IC.Nucleotide.pdb.save IC.Nucleotide.pdb
}

# Link nucleic acid base + sugar + phosphate, fix charges.
DNAcharge() {
  cat > cpptraj.in <<EOF
parm DDD.names.mol2
loadcrd DDD.names.mol2 name Sugar parm DDD.names.mol2
parm MP1.names.mol2
loadcrd MP1.names.mol2 name Phos parm MP1.names.mol2
parm ADD.names.mol2
loadcrd ADD.names.mol2 name Base parm ADD.names.mol2
graft \
  tgt Base \
  tgtcharge 0.116 \
  src Sugar \
  srccharge 0.2933 \
  name BaseSugar \
  tgtfitmask @N1,C1,H1 \
  srcfitmask @C4,C3,O2 \
  tgtmask !(@C1,H1,H6,H7) \
  srcmask !(@C4,H6,H7,H8,H12,H11) \
  bond @N1,@C3
charge crdset BaseSugar *
graft \
  tgt BaseSugar \
  tgtcharge 0.4093 \
  src Phos \
  srccharge -1.4093 \
  name Nucleotide \
  tgtfitmask @C1,O1,H1 \
  srcfitmask @C3,O3,P1 \
  tgtmask !(@O1,H1) \
  srcmask !(@C3,H4,H5,H6,O1,C1,H1,H2,H3) \
  bond @C1,@O3
charge crdset Nucleotide *
crdout Nucleotide Nucleotide.charge.mol2
quit
EOF
  RunCpptraj "$TESTNAME, Construct Nucleic Acid and adjust charge"
  DoTest Nucleotide.charge.mol2.save Nucleotide.charge.mol2
}

TyrPry
TyrPryIC
DNA
DNAic
DNAcharge

EndTest
exit 0
