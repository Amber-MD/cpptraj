#!/bin/bash

. ../MasterTest.sh

CleanFiles rst.in output

INPUT='-i rst.in'

cat > rst.in <<EOF
parm ../tz2.parm7
reference ../tz2.rst7
# Single atom distance restraint
rst :1@CA :12@CA r1 4.0 r2 5.0 r3 5.0 r4 6.0 rk2 2.0 rk3 2.0 out output
# Multiple atom distance restraint using reference
rst :1&!@H= :12&!@H= reference offset 1.0 rk2 1.0 rk3 1.0 out output
# Angle restraint
rst :1@CA :5@CA :12@CA r1 80 r2 85 r3 90 r4 95 rk2 10 rk3 10 out output
# Torsion restraint
rst :4@CA :4@C :5@N :5@CA r1 175 r2 180 r3 180 r4 185 rk2 100 rk3 100 out output
EOF
RunCpptraj "Generate restraints (rst) test"
DoTest output.save output

EndTest
exit 0
