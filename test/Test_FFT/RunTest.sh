#!/bin/bash

. ../MasterTest.sh

CleanFiles fft.in fft.dat fft.agr sin1.agr fft?.agr

INPUT="-i fft.in"
cat > fft.in <<EOF
readdata sin_4pt.agr
write fft.dat sin_4pt.agr xmin 0 xstep 0.0001
runanalysis fft sin_4pt.agr out fft.agr dt 0.025

readdata sin_100pt.agr
runanalysis fft sin_100pt.agr out fft1.agr dt 0.001

readdata sin_128pt.agr
runanalysis fft sin_128pt.agr out fft2.agr dt 0.00078125

readdata cos2piNover10.agr
runanalysis fft cos2piNover10.agr out fft3.agr dt 1.0

readdata cos10t_cos5t.agr
runanalysis fft cos10t_cos5t.agr out fft4.agr dt .0033333333

quit
EOF
RunCpptraj "FFT test."
DoTest fft.agr.save fft.agr
DoTest fft1.agr.save fft1.agr
DoTest fft2.agr.save fft2.agr
DoTest fft3.agr.save fft3.agr
DoTest fft4.agr.save fft4.agr

EndTest

exit 0
