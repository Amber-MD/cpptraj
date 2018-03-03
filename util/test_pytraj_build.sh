#!/usr/bin/env bash

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../"; pwd;)
DOCKER_IMAGE=hainm/pytraj-build-box:17.0

# not clone pytraj inside the docker image
# https://githubengineering.com/crypto-removal-notice/
git clone http://github.com/amber-md/pytraj

docker info
cat << EOF | docker run -i \
                        -v ${FEEDSTOCK_ROOT}:/cpptraj \
                        -a stdin -a stdout -a stderr \
                        ${DOCKER_IMAGE}\
                        bash || exit $?


set -x
cd /cpptraj/

bash configure -shared -openmp gnu
make libcpptraj -j4
export CPPTRAJHOME=/cpptraj

cd pytraj
/opt/python/cp35-cp35m/bin/python setup.py install
# make sure we can run very simple test
/opt/python/cp35-cp35m/bin/python run_tests.py -s
EOF
