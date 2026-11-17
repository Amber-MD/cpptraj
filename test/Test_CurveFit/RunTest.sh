#!/bin/bash

. ../MasterTest.sh

CleanFiles cf.in curve.dat curve1.dat Kcurve.dat PKcurve.dat curve2.dat \
           Results.dat results1.dat Cubic.dat results.cubic.dat

INPUT="-i cf.in"
# General test
cat > cf.in <<EOF
readdata Data.dat index 1
runanalysis curvefit Data.dat "FitY = (A0 * exp(X * A1)) + (A2 * exp(X * A3))" \
  A0=1 A1=-1 A2=1 A3=-1 \
  out curve.dat tol 0.0001 maxit 5000

createset "Y = 3*X^3 + 2*X^2 + 0.5*X" xmin -5 xstep 0.1 nx 100
#writedata Cubic.dat.save Y
runanalysis curvefit Y "Cubic = (A0*X^3) + (A1*X^2) + (A2*X)" \
  A0=1 A1=1 A2=1 \
  out Cubic.dat resultsout results.cubic.dat tol 0.0001 maxit 5000
EOF
RunCpptraj "Curve fitting test."
DoTest curve.dat.save curve.dat
# Add 0.03 cushion for Chi squared, uncertainty, and RMS P.E.
DoTest results.cubic.dat.save results.cubic.dat -a 0.03

# Mexp test
cat > cf.in <<EOF
readdata Data.dat index 1
runanalysis curvefit Data.dat nexp 2 name FitY \
  A0=1 A1=-1 A2=1 A3=-1 \
  out curve1.dat resultsout results1.dat tol 0.0001 maxit 5000

runanalysis curvefit Data.dat nexp 2 name KFitY form mexpk \
  A0=0 A1=1 A2=-1 A3=1 A4=-1 \
  out Kcurve.dat tol 0.0001 maxit 5000

runanalysis curvefit Data.dat nexp 2 name PKFitY form mexpk_penalty \
  A0=0 A1=.5 A2=-1 A3=.5 A4=-1 \
  out PKcurve.dat tol 0.0001 maxit 5000 resultsout Results.dat
EOF
RunCpptraj "Curve fitting multi-exponential tests."
DoTest curve.dat.save curve1.dat
# Order of exponentials in results can change, so just check curve.
#DoTest results1.dat.save results1.dat
DoTest Kcurve.dat.save Kcurve.dat -a 0.0002
# Differences in windows seem like round-off
DoTest PKcurve.dat.save PKcurve.dat -a 0.0002
# Order of exponentials in results can change, so just check curve.
#DoTest Results.dat.save Results.dat -a 0.2

# Custom X output range.
cat > cf.in <<EOF
readdata Data.dat index 1
runanalysis curvefit Data.dat nexp 2 name FitY \
  A0=1 A1=-1 A2=1 A3=-1 \
  out curve2.dat tol 0.0001 maxit 5000 \
  outxbins 1000 outxmin 0.0 outxmax 500
EOF
RunCpptraj "Curve fitting, custom output range."
DoTest curve2.dat.save curve2.dat

EndTest
exit 0
