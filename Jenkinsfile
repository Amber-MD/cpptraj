#!/usr/bin/env groovy
// Builds on a variety of platforms with a variety of options

pipeline {
    // Each of the build steps will define their own agents. We will use docker wherever
    // We can, but there is no Docker container available with either the Intel or PGI
    // compilers, so we'll have to make do with running on bare metal in those cases
    agent none

    options {
        skipDefaultCheckout()
        buildDiscarder(logRotator(numToKeepStr: "10"))
        timestamps()
    }

    stages {
        stage("Checkout and stash source code") {
            agent { label "linux" }

            steps {
                script {
                    env.GIT_COMMIT = checkout(scm).GIT_COMMIT
                }

                // We need to *not* exclude default excludes so we keep the .git
                // directory and its relevant information
                stash includes: "**", name: "source", useDefaultExcludes: false
            }
            post { cleanup { deleteDir() } }
        }

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
                        docker {
                            image 'ambermd/cpu-build:latest'
                            alwaysPull true
                        }
                    }

                    steps {
                        unstash "source"
                        sh "./configure --with-netcdf --with-fftw3 gnu"
                        sh "make -j4 install"
                        sh "cd test && make test.showerrors"
                    }

                    post { cleanup { deleteDir() } }
                }
                stage("Linux Intel Serial Build") {
                    agent {
                        docker {
                            image 'ambermd/cpu-build:latest'
                            alwaysPull true
                            label "docker && intel"
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
                        unstash "source"
                        sh "./configure --with-netcdf -mkl intel"
                        sh "make -j4 install"
                        sh "cd test && make test.showerrors"
                    }
                    post { cleanup { deleteDir() } }
                }
                stage("Linux PGI serial build") {
                    agent {
                        docker {
                            label "pgi && Batwoman"
                            image 'ambermd/cpu-build:latest'
                            alwaysPull true
                            // Pull the licensed PGI compilers from the host machine (must have
                            // the compilers installed in /opt/pgi)
                            args "-v /opt/pgi:/opt/pgi"
                        }
                    }

                    steps {
                        unstash "source"
                        script {
                            try {
                                sh "./configure --with-netcdf --with-fftw3 pgi"
                                sh "make -j6 install"
                                sh "cd test && make test.showerrors"
                            } catch (error) {
                                echo "PGI BUILD AND/OR TEST FAILED"
                            }
                        }
                    }
                    post { cleanup { deleteDir() } }
                }
                stage("Linux GNU parallel build") {
                    agent {
                        docker {
                            image 'ambermd/cpu-build:latest'
                            alwaysPull true
                        }
                    }

                    steps {
                        unstash "source"
                        sh "./configure --with-netcdf --with-fftw3 -mpi gnu"
                        sh "make -j4 install"
                        sh "make -e DO_PARALLEL='mpiexec -n 2' check"
                        sh "make -e DO_PARALLEL='mpiexec -n 4' check"
                    }
                    post { cleanup { deleteDir() } }
                }
                stage("Linux CUDA build") {
                    agent {
                        docker {
                            image 'ambermd/gpu-build:latest'
                            alwaysPull true
                            label "docker && cuda"
                        }
                    }

                    steps {
                        unstash "source"
                        sh "./configure --with-netcdf --with-fftw3 -cuda gnu"
                        sh "make -j4 install"
                        sh "make -e OPT=cuda check"
                    }
                    post { cleanup { deleteDir() } }
                }
                stage("Linux GNU OpenMP build") {
                    agent {
                        docker {
                            image 'ambermd/cpu-build:latest'
                            alwaysPull true
                        }
                    }

                    steps {
                        unstash "source"
                        sh "./configure --with-netcdf --with-fftw3 -openmp gnu"
                        sh "make -j4 install"
                        sh "make -e OPT=openmp check"
                    }
                    post { cleanup { deleteDir() } }
                }
            }
        }

        stage("Build libcpptraj docker container") {
            agent { label "linux && docker" }

            environment {
                DOCKER_IMAGE_TAG = env.BRANCH_NAME == "master" ? "master" : env.GIT_COMMIT
            }

            steps {
                unstash "source"
                echo "Building and pushing ambermd/libcpptraj:${env.DOCKER_IMAGE_TAG}"
                script {
                    def image = docker.build("ambermd/libcpptraj:${env.DOCKER_IMAGE_TAG}",
                                             "-f ./devtools/ci/jenkins/Dockerfile.libcpptraj .")
                    docker.withRegistry("", "amber-docker-credentials") {
                        image.push()
                    }
                }
            }
            post { cleanup { deleteDir() } }
        }

        stage("Check pytraj with this version of libcpptraj") {
            environment {
                DOCKER_IMAGE_TAG = env.BRANCH_NAME == "master" ? "master" : env.GIT_COMMIT
            }

            steps {
                build job: '/amber-github/pytraj',
                           parameters: [string(name: 'LIBCPPTRAJ_IMAGE_TAG', value: "${DOCKER_IMAGE_TAG}"),
                                        string(name: 'BRANCH_TO_BUILD', value: 'master')]
            }
        }
    }
}
