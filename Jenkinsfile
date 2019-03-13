// Builds on a variety of platforms with a variety of options

pipeline {
    // Each of the build steps will define their own agents. We will use docker wherever
    // We can, but there is no Docker container available with either the Intel or PGI
    // compilers, so we'll have to make do with running on bare metal in those cases
    agent none

    stages {
        stage("Build and test") {
            parallel {
                stage("Linux GNU serial build") {
                    agent {
                        dockerfile {
                            dir "devtools/ci/jenkins"
                        }
                    }

                    environment {
                        COMPILER = "gnu"
                        OPERATING_SYSTEM = "linux"
                    }

                    steps {
                        sh "bash -ex devtools/ci/jenkins/install.sh"
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
                        COMPILER = "intel"
                        OPERATING_SYSTEM = "linux"
                        COMPILER_FLAGS = "-mkl"
                    }

                    steps {
                        sh "bash -ex devtools/ci/jenkins/install.sh"
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

                    environment {
                        COMPILER = "pgi"
                        OPERATING_SYSTEM = "linux"
                    }

                    steps {
                        sh "bash -ex devtools/ci/jenkins/install.sh"
                    }
                }
                stage("Linux GNU parallel build") {
                    // Insert parallel build here
                    agent {
                        dockerfile {
                            dir "devtools/ci/jenkins"
                        }
                    }

                    environment {
                        COMPILER = "gnu"
                        OPERATING_SYSTEM = "linux"
                        COMPILER_FLAGS = "-mpi"
                    }

                    steps {
                        sh "bash -ex devtools/ci/jenkins/install.sh"
                    }
                }
                stage("Linux CUDA build") {
                    agent {
                        dockerfile {
                            dir "devtools/ci/jenkins"
                            additionalBuildArgs "--build-arg BASEIMAGE=nvidia/cuda:10.0-devel-ubuntu18.04"
                        }
                    }

                    environment {
                        COMPILER = "gnu"
                        COMPILER_FLAGS = "-cuda"
                        OPERATING_SYSTEM = "linux"
                    }

                    steps {
                        sh "bash -ex devtools/ci/jenkins/install.sh"
                    }
                }
                stage("macOS build") {
                    // No docker on Macs :'(
                    agent { label "mac" }

                    environment {
                        OPERATING_SYSTEM = "macOS"
                    }

                    steps {
                        sh "bash -ex devtools/ci/jenkins/install.sh"
                    }
                }
                stage("Linux GNU OpenMP build") {
                    agent {
                        dockerfile {
                            dir "devtools/ci/jenkins"
                        }
                    }

                    environment {
                        COMPILER = "gnu"
                        OPERATING_SYSTEM = "linux"
                        COMPILER_FLAGS = "-openmp"
                    }

                    steps {
                        sh "bash -ex devtools/ci/jenkins/install.sh"
                    }
                }
            }
        }
    }
}