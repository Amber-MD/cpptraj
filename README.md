CPPTRAJ
=======

Fast, parallelized molecular dynamics trajectory data analysis.

Build Status
=============
* Travis-CI: [![Travis Build Status](https://travis-ci.org/Amber-MD/cpptraj.svg?branch=master)](https://travis-ci.org/Amber-MD/cpptraj)
* AppVeyor: [![Windows Build Status](https://ci.appveyor.com/api/projects/status/github/Amber-MD/cpptraj?branch=master&svg=true&retina=true)](https://ci.appveyor.com/project/drroe/cpptraj-aof9y/branch/master)
* Jenkins: [![Jenkins Build Status](https://jenkins.jasonswails.com/buildStatus/icon?job=amber-github%2Fcpptraj-master?style=plastic)](https://jenkins.jasonswails.com/job/amber-github/job/cpptraj-master/)

Description
===========

CPPTRAJ is a program designed to process and analyze molecular dynamics
trajectories and relevant data sets derived from their analysis. CPPTRAJ
supports many popular MD software packages including [Amber](http://ambermd.org/),
[CHARMM](https://www.charmm.org/charmm/), [Gromacs](http://www.gromacs.org/),
and [NAMD](http://www.ks.uiuc.edu/Research/namd/).

CPPTRAJ is also distributed as part of the freely available AmberTools
software package. The official AmberTools release version of CPPTRAJ can
be found at the [Amber website](http://ambermd.org).

For those wanting to use CPPTRAJ in their Python scripts, see
[Pytraj](https://github.com/Amber-MD/pytraj).

See what's [new in CPPTRAJ](https://github.com/Amber-MD/cpptraj/wiki). For those
just starting out you may want to check out some [CPPTRAJ tutorials](http://ambermd.org/tutorials/TrajectoryAnalysis.php#cpptraj).

For more information (or to cite CPPTRAJ) see the following publication:

[Daniel R. Roe and Thomas E. Cheatham, III, "PTRAJ and CPPTRAJ: Software for
  Processing and Analysis of Molecular Dynamics Trajectory Data". J. Chem.
  Theory Comput., 2013, 9 (7), pp 3084-3095](http://pubs.acs.org/doi/abs/10.1021/ct400341p).

For more information regarding trajectory/ensemble parallelism via MPI in
CPPTRAJ see the following publication:

[Daniel R. Roe and Thomas E. Cheatham, III, \"Parallelization of CPPTRAJ enables
  large scale analysis of molecular dynamics trajectory data\". J. Comp.
  Chem., 2018, DOI: 10.1002/jcc25382](https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.25382).


Disclaimer and Copyright
========================

CPPTRAJ is Copyright (c) 2010-2019 Daniel R. Roe.
The terms for using, copying, modifying, and distributing CPPTRAJ are
specified in the file LICENSE.

Documentation
=============
  The `/doc` subdirectory contains PDF and LyX versions of the CPPTRAJ manual.
An HTML version can be found [here](https://amber-md.github.io/cpptraj/). There
is also limited help for commands in interactive mode via `help [<command>]`;
`help` with no arguments lists all known commands.

  Code documentation can be generated via Doxygen by typing `make docs`. This
will install HTML and Latex documentation at `doc/html/index.html` and in
the `doc/latex` respectively. A limited developers guide is available in
Lyx/PDF formats in the `doc/` subdirectory and in HTML format
[here](https://amber-md.github.io/cpptraj/).

Installation & Testing
======================
Run `./configure --help` for a short list of configure options. `./configure --full-help`
will list all available configure options. For full functionality, CPPTRAJ makes use of
the following libraries:

* NetCDF
* BLAS
* LAPACK
* ARPACK (now bundled with CPPTRAJ)
* Bzip2
* Gzip
* Parallel NetCDF (-mpi build only, for NetCDF trajectory output in parallel)
* CUDA (-cuda build only)
* FFTW (mostly optional; required for PME functionality)

`./configure gnu` should be adequate to set up compilation for most systems.
For systems without BLAS/LAPACK/ARPACK and/or NetCDF libraries installed,
the `-amberlib` flag can be specified to use the ones already compiled in
an AmberTools installation (`$AMBERHOME` must be set), e.g.
`./configure -amberlib gnu`. For multicore systems, the `-openmp` flag can
be specified to enable OpenMP parallelization, e.g. `./configure -openmp gnu`.
An MPI-parallelized version of CPPTRAJ can also be built using the `-mpi` flag.
CPPTRAJ can be built with both MPI and OpenMP; when running this build users
should take care to properly set OMP_NUM_THREADS if using more than 1 MPI
thread per node. A CUDA build is now also available via the `-cuda` flag.
By default CPPTRAJ will be configured for multiple shader models; to restrict
the CUDA build to a single shader model use the SHADER_MODEL environment variable.
Any combination of `-cuda`, `-mpi`, and `-openmp` may be used.

The configure script by default sets everything up to link dynamically. The
`-static` flag can be used to force static linking. If linking errors are
encountered you may need to specify library locations using the `--with-LIB=`
options. For example, to use NetCDF compiled in `/opt/netcdf` use the option
`--with-netcdf=/opt/netcdf`. Alternatively, individual libraries can be
disabled with the `-no<LIB>` options. The `-libstatic` flag
can be used to static link only libraries that have been specified.

After `configure` has been successfully run, `make install` will
compile and place the cpptraj binary in the `bin/` subdirectory. Note that
on multithreaded systems `make -j X install` (where X is an integer > 1
and less than the max # cores on your system) will run much faster.
After installation, It is highly recommended that `make check` be run as
well to test the basic functionality of CPPTRAJ.

CPPTRAJ Authors
===============
**Lead Author:** Daniel R. Roe (<daniel.r.roe@gmail.com>)
Laboratory of Computational Biology
National Heart Lung and Blood Institute
National Institutes of Health, Bethesda, MD.

  CPPTRAJ began as a C++ rewrite of PTRAJ by Thomas E. Cheatham, III
(Department of Medicinal Chemistry, University of Utah, Salt Lake City,
 UT, USA) and many routines from PTRAJ were adapted for
use in CPPTRAJ, including code used in the following classes:
Analysis\_CrankShaft, Analysis\_Statistics, Action\_DNAionTracker,
Action\_RandomizeIons, Action\_Principal, Action\_Grid, GridAction,
Action\_Image, and ImageRoutines.

## Contributors to CPPTRAJ

* James Maier (Stony Brook University, Stony Brook, NY, USA)
Code for calculating J-couplings (used in Action\_Jcoupling).

* Jason M. Swails (University of Florida, Gainesville, FL, USA)
Action\_LIE, Analysis\_RunningAvg, Action\_Volmap, Grid OpenDX output.

* Jason M. Swails (University of Florida, Gainesville, FL, USA)
Guanglei Cui (GlaxoSmithKline, Upper Providence, PA, USA)
Action\_SPAM.

* Mark J. Williamson (Unilever Centre for Molecular Informatics, Department of Chemistry, Cambridge, UK)
Action\_GridFreeEnergy.

* Hannes H. Loeffler (STFC Daresbury, Scientific Computing Department, Warrington, WA4 4AD, UK)
Action\_Density, Action\_OrderParameter, Action\_PairDist.

* Crystal N. Nguyen (University of California, San Diego)
Romelia F. Salomon (University of California, San Diego)
Original Action\_Gist.

* Pawel Janowski (Rutgers University, NJ, USA)
Normal mode wizard (nmwiz) output, original code for ADP calculation in Action\_AtomicFluct.

* Zahra Heidari (Faculty of Chemistry, K. N. Toosi University of Technology, Tehran, Iran)
Original code for Analysis\_Wavelet.

* Chris Lee (University of California, San Diego)
Support for processing force information in NetCDF trajectories.

* Steven Ramsey (CUNY Lehman College, Bronx, NY)
Enhancements to entropy calculation in original Action\_Gist.

* Amit Roy (University of Utah, UT)
Code for the CUDA version of the 'closest' Action.

* Andrew Simmonett (National Institutes of Health)
Code for the reciprocal part of the particle mesh Ewald calculation (electrostatic and Lennard-Jones).

* Christina Bergonzo (National Institute of Standards and Technology, Gaithersburg, MD)
Fixes and improvements to nucleic acid dihedral angle definitions (DihedralSearch).

* David S. Cerutti (Rutgers University, Piscataway, NJ, USA)
Original code for the 'xtalsymm' Action.

#### Various Contributions
* David A. Case (Rutgers University, Piscataway, NJ, USA)
* Hai Nguyen (Rutgers University, Piscataway, NJ, USA)
* Robert T. McGibbon (Stanford University, Stanford, CA, USA)

## Code in CPPTRAJ that originated in PTRAJ

* Holger Gohlke (Heinrich-Heine-University, Düsseldorf, Germany)
Alrun N. Koller (Heinrich-Heine-University, Düsseldorf, Germany)
Original implementation of matrix/vector functionality in PTRAJ, including matrix diagonalization, IRED analysis, eigenmode analysis, and vector time correlations.

* Holger Gohlke (Heinrich-Heine-University, Düsseldorf, Germany)
Original code for DSSP (secstruct).

* Michael Crowley (University of Southern California, Los Angeles, CA, USA)
Original code for dealing with truncated octahedral unit cells.

* Viktor Hornak (Merck, NJ, USA)
Original code for mask expression parser.

* John Mongan (UCSD, San Diego, CA, USA)
Original implementation of the Amber NetCDF trajectory format.

* Hannes H. Loeffler (STFC Daresbury, Scientific Computing Department, Warrington, WA4 4AD, UK)
Diffusion calculation code adapted for use in Action\_STFC\_Diffusion.

External libraries bundled with CPPTRAJ
=======================================

* CPPTRAJ makes use of the [GNU readline](https://tiswww.case.edu/php/chet/readline/rltop.html) library for the interactive command line.

* CPPTRAJ uses the [ARPACK](https://www.caam.rice.edu//software/ARPACK/) library to calculate eigenvalues/eigenvectors from large sparse matrices.

* CPPTRAJ uses the [xdrfile](http://www.gromacs.org/Developer\_Zone/Programming\_Guide/XTC\_Library) library for reading XTC file; specifically a somewhat updated version from [MDTRAJ](https://github.com/mdtraj/mdtraj) that includes some bugfixes and enhancements. See `src/xdrfile/README` for details.

* The reciprocal part of the PME calculation is handled by the [helPME](https://github.com/andysim/helpme) library by Andy Simmonett.
