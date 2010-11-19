#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles dummy.rst7.0 dummy.pdb strip.ChainA-tip3p.parm7

# Test 1
INPUT="-i strip.1frame.in"
RunCpptraj "One frame strip command test."
DoTest dummy.pdb.save dummy.pdb
DoTest dummy.rst7.0.save dummy.rst7.0
DoTest strip.ChainA-tip3p.parm7.save strip.ChainA-tip3p.parm7 
CheckTest

EndTest

exit 0
