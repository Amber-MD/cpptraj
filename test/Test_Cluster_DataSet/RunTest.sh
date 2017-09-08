#!/bin/bash

. ../MasterTest.sh

CleanFiles ds.in *.summary.dat *.info.dat *.gnu *.d1.c1.dat \
           eps_v_n.dat twodihds.kmeans.info.dat onedihds.kmeans.info.dat \
           cvt.dat

INPUT="-i ds.in"
TESTNAME='Clustering via datasets tests'
Requires netcdf
OneDS() {
  TOP=../tz2.parm7
  cat > ds.in <<EOF
trajin ../tz2.nc
distance d1 :1 :13
createcrd crd1
#debug analysis 2
cluster crdset crd1 c1 data d1 clusters 5 epsilon 4.0 out oneds.gnu summary oneds.summary.dat info oneds.info.dat gracecolor epsilonplot eps_v_n.dat \
        clustersvtime cvt.dat cvtwindow 10
create oneds.d1.c1.dat d1 c1
EOF
  RunCpptraj "Clustering, One DataSet"
  DoTest oneds.info.dat.save oneds.info.dat
  DoTest eps_v_n.dat.save eps_v_n.dat
  DoTest cvt.dat.save cvt.dat
}

TwoDS() {
  TOP=../tz2.parm7
  cat > ds.in <<EOF
trajin ../tz2.nc
distance d1 :1 :13
hbond hb1
#debug analysis 2
cluster c1 data d1,hb1[UU] clusters 5 epsilon 4.0 out twods.gnu summary twods.summary.dat info twods.info.dat gracecolor
create twods.d1.c1.dat d1 hb1[UU] c1
EOF
  RunCpptraj "Clustering, Two DataSets"
  DoTest twods.info.dat.save twods.info.dat
}

OneDihDS() {
  TOP=../tz2.parm7
  cat > ds.in <<EOF
trajin ../tz2.nc
dihedral gly7phi :6@C :7@N :7@CA :7@C
#debug analysis 2
cluster c1 data gly7phi clusters 5 out onedihds.gnu summary onedihds.summary.dat info onedihds.info.dat gracecolor
create onedihds.d1.c1.dat gly7phi c1
EOF
  RunCpptraj "Clustering, dihedral DataSet"
  DoTest onedihds.info.dat.save onedihds.info.dat 
}

TwoDihDS() {
  TOP=../tz2.parm7
  cat > ds.in <<EOF
trajin ../tz2.nc
dihedral gly7phi :6@C :7@N :7@CA :7@C
dihedral gly7psi :7@N :7@CA :7@C :8@N
createcrd crd1
#debug analysis 2
cluster crdset crd1 c1 data gly7phi,gly7psi clusters 5 out twodihds.gnu summary twodihds.summary.dat info twodihds.info.dat gracecolor
create twodihds.d1.c1.dat gly7phi gly7psi c1 noxcol
EOF
  RunCpptraj "Clustering, two dihedral DataSets"
  DoTest twodihds.info.dat.save twodihds.info.dat
}

OneDihDSKmeans() {
  TOP=../tz2.parm7
  cat > ds.in <<EOF
trajin ../tz2.nc
dihedral gly7phi :6@C :7@N :7@CA :7@C
#debug analysis 1
cluster c1 kmeans randompoint kseed 10 data gly7phi clusters 5 out onedihds.gnu summary onedihds.summary.dat info onedihds.kmeans.info.dat gracecolor
create onedihds.d1.c1.dat gly7phi c1
EOF
  RunCpptraj "Clustering, dihedral DataSet, K-means"
  DoTest onedihds.kmeans.info.dat.save onedihds.kmeans.info.dat
}

TwoDihDSKmeans() {
  TOP=../tz2.parm7
  cat > ds.in <<EOF
trajin ../tz2.nc
dihedral gly7phi :6@C :7@N :7@CA :7@C
dihedral gly7psi :7@N :7@CA :7@C :8@N
createcrd crd1
#debug analysis 2
cluster crdset crd1 c1 kmeans data gly7phi,gly7psi clusters 5 out twodihds.gnu summary twodihds.summary.dat info twodihds.kmeans.info.dat gracecolor
EOF
  RunCpptraj "Clustering, two dihedral DataSets, K-means"
  DoTest twodihds.kmeans.info.dat.save twodihds.kmeans.info.dat
}

OneDS
TwoDS
OneDihDS
TwoDihDS
OneDihDSKmeans
TwoDihDSKmeans

EndTest

exit 0
