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
        timeout time: 2, unit: 'HOURS'
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
              stage("Linux GNU Unit Tests") {
                    agent {
                        docker {
                            image 'ambermd/cpu-build:latest'
                            alwaysPull true
                        }
                    }

                    steps {
                        unstash "source"
                        sh "./configure --with-netcdf --with-fftw3 gnu"
                        sh "cd unitTests && make test.all"
                    }

                    post { cleanup { deleteDir() } }
                }
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
                    when { expression { return false } }
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
                            // Until we figure out why PGI compilers fail on the master Jenkins node,
                            // force them to run on Intel CPU slaves
                            label "pgi && intel-cpu"
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
                                sh """#!/bin/sh -ex
                                unset CUDA_HOME
                                ./configure --with-netcdf --with-fftw3 pgi
                                make -j6 install
                                cd test && make test.showerrors
                                """
                            } catch (error) {
                                echo "PGI BUILD AND/OR TEST FAILED"
                                try {
                                    pullRequest.comment("The PGI build in Jenkins failed.")
                                } catch (err2) {
                                    echo "Could not post a PR comment: ${err2}"
                                }
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
                        sh "./configure --with-netcdf --with-fftw3 -mpi --buildlibs gnu"
                        sh "make -j4 install"
                        sh "make -e DO_PARALLEL='mpiexec -n 2' verbosecheck"
                        sh "make -e DO_PARALLEL='mpiexec -n 4' verbosecheck"
                    }
                    post { cleanup { deleteDir() } }
                }
                stage("Linux CUDA build") {
                    agent {
                        docker {
                            image 'ambermd/gpu-build:latest'
                            alwaysPull true
                            label "docker && cuda"
                            args '--gpus all'
                        }
                    }

                    steps {
                        unstash "source"
                        sh "./configure --with-netcdf --with-fftw3 -cuda gnu"
                        sh "make -j4 install"
                        sh "make -e OPT=cuda verbosecheck"
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
                        sh "make -e OPT=openmp verbosecheck"
                    }
                    post { cleanup { deleteDir() } }
                }
            }
        }

        stage("Post-test steps") {
            parallel {
                stage("Check pytraj") {
                    agent none
                    environment {
                        DOCKER_IMAGE_TAG = "${env.BRANCH_NAME == "master" ? "master" : env.GIT_COMMIT}"
                    }

                    stages {
                        stage("Build libcpptraj docker container") {
                            agent { label "linux && docker" }

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
                            steps {
                                build job: '/amber-github/pytraj',
                                           parameters: [string(name: 'LIBCPPTRAJ_IMAGE_TAG', value: env.DOCKER_IMAGE_TAG),
                                                        string(name: 'BRANCH_TO_BUILD', value: 'master')]
                            }
                        }
                    }
                }
                stage("Publish the manual") {
                    stages {
                        stage("Build the manual") {
                            agent {
                                docker {
                                    image 'ambermd/lyx:latest'
                                    alwaysPull true
                                }
                            }

                            steps {
                                unstash "source"
                                sh """#!/bin/sh -ex
                                    make docs
                                    cd doc
                                    lyx -batch --export pdf2 CpptrajManual.lyx
                                    lyx -batch --export pdf2 CpptrajDevelopmentGuide.lyx
                                """
                                stash includes: "doc/**", name: "documentation"
                            }

                            post { cleanup { deleteDir() } }
                        }
                        stage("Publish the manual") {
                            agent { label 'linux' }

                            steps {
                                unstash 'documentation'
                                // Eventually it would be nice to do something better than simply archive
                                // and make the artifacts available from Jenkins.
                                archiveArtifacts 'doc/CpptrajManual.pdf,doc/CpptrajDevelopmentGuide.pdf'
                                // It would be nice to get this in a better place, but for now this
                                // will suffice.
                                publishHTML([
                                    allowMissing: false,
                                    alwaysLinkToLastBuild: false,
                                    keepAll: false,
                                    reportDir: 'doc/html',
                                    reportFiles: 'index.html',
                                    reportName: 'Doxygen Documentation',
                                    reportTitles: ''])
                            }

                            post { cleanup { deleteDir() } }
                        }
                    }
                }
            }
        }
    }
}
