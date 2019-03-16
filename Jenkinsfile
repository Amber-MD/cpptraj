// Builds on a variety of platforms with a variety of options

pipeline {
    // Each of the build steps will define their own agents. We will use docker wherever
    // We can, but there is no Docker container available with either the Intel or PGI
    // compilers, so we'll have to make do with running on bare metal in those cases
    agent none

    stages {
        stage("Build and test") {
            when {
                // Testing needs to be done as a merge gate. Since this is a
                // somewhat pricey operation in the number of executors used,
                // assume sufficient testing was already done for the master
                // branch
                not {
                    branch "master"
                }
            }
            parallel {
                stage("Linux GNU serial build") {
                    agent {
                        dockerfile {
                            dir "devtools/ci/jenkins"
                        }
                    }

                    steps {
                        sh "./configure --with-netcdf --with-fftw3 gnu"
                        sh "make -j4 install"
                        sh "make check"
                    }
                }
                stage("Linux Intel Serial Build") {
                    agent {
                        dockerfile {
                            dir "devtools/ci/jenkins"
                            label "intel"
                            // There's no way to have a docker container installed with a licensed
                            // copy of the Intel compilers, so we need to mount it from the host
                            // machine. This introduces the constraint that *all* of the Jenkins
                            // slaves need to have the Intel compilers available at /opt/intel
                            args "-v /opt/intel:/opt/intel"
                        }
                    }

                    environment {
                        MKL_HOME = "/opt/intel/compilers_and_libraries/linux/mkl"
                    }

                    steps {
                        sh "./configure --with-netcdf -mkl intel"
                        sh "make -j4 install"
                        sh "make check"
                    }
                }
                stage("Linux PGI serial build") {
                    agent {
                        dockerfile {
                            dir "devtools/ci/jenkins"
                            label "pgi"
                            // Pull the licensed PGI compilers from the host machine (must have
                            // the compilers installed in /opt/pgi)
                            args "-v /opt/pgi:/opt/pgi"
                        }
                    }

                    steps {
                        sh "./configure --with-netcdf --with-fftw3 pgi"
                        sh "make -j6 install"
                        sh "make check"
                    }
                }
                stage("Linux GNU parallel build") {
                    agent {
                        dockerfile {
                            dir "devtools/ci/jenkins"
                        }
                    }

                    steps {
                        sh "./configure --with-netcdf --with-fftw3 -mpi gnu"
                        sh "make -j4 install"
                        sh "make -e DO_PARALLEL='mpiexec -n 2' check"
                        sh "make -e DO_PARALLEL='mpiexec -n 4' check"
                    }
                }
                stage("Linux CUDA build") {
                    agent {
                        dockerfile {
                            dir "devtools/ci/jenkins"
                            filename "Dockerfile.cuda"
                        }
                    }

                    steps {
                        sh "./configure --with-netcdf --with-fftw3 -cuda gnu"
                        sh "make -j4 install"
                        sh "make -e OPT=cuda check"
                    }
                }
                stage("Linux GNU OpenMP build") {
                    agent {
                        dockerfile {
                            dir "devtools/ci/jenkins"
                        }
                    }

                    steps {
                        sh "./configure --with-netcdf --with-fftw3 -openmp gnu"
                        sh "make -j4 install"
                        sh "make -e OPT=openmp check"
                    }
                }
            }
        }

        stage("Build libcpptraj docker container") {
            agent { label "linux && docker" }

            environment {
                DOCKER_IMAGE_TAG = "${env.BRANCH_NAME == "master" ? "master" : env.ghprbActualCommit}"
            }

            steps {
                echo "Building and pushing ambermd/libcpptraj:${env.DOCKER_IMAGE_TAG}"
                script {
                    def image = docker.build("ambermd/libcpptraj:${env.DOCKER_IMAGE_TAG}",
                                             "-f ./devtools/ci/jenkins/Dockerfile.libcpptraj .")
                    docker.withRegistry("", "amber-docker-credentials") {
                        image.push()
                    }
                }
            }
        }

        stage("Check pytraj with this version of libcpptraj") {
            environment {
                DOCKER_IMAGE_TAG = "${env.BRANCH_NAME == "master" ? "master" : env.ghprbActualCommit}"
            }

            steps {
                build job: '/amber-github/pytraj',
                           parameters: [string(name: 'LIBCPPTRAJ_IMAGE_TAG', value: "${DOCKER_IMAGE_TAG}"),
                                        string(name: 'BRANCH_TO_BUILD', value: 'master')]
            }
        }
    }
}
