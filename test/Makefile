# Test makefile for cpptraj
test: test.all 

test.general:
	@-cd Test_General && ./RunTest.sh $(OPT)
	@-cd Test_Parallel && ./RunTest.sh $(OPT)

test.strip:
	@-cd Test_Strip && ./RunTest.sh $(OPT)

test.brokentraj:
	@-cd Test_BrokenTraj && ./RunTest.sh $(OPT)

test.trajinoffset:
	@-cd Test_TrajinOffset && ./RunTest.sh $(OPT)

test.remdnetcdf:
	@-cd Test_RemdTraj && ./RunTest.sh $(OPT)

test.framerange:
	@-cd Test_trajout_framerange && ./RunTest.sh $(OPT)

test.perres:
	@-cd Test_PerResRMSD && ./RunTest.sh $(OPT)

test.dssp:
	@-cd Test_DSSP && ./RunTest.sh $(OPT)

test.center:
	@-cd Test_Center && ./RunTest.sh $(OPT)

test.hbond:
	@-cd Test_Hbond && ./RunTest.sh $(OPT)

test.image:
	@-cd Test_Image && ./RunTest.sh $(OPT)

test.netcdftraj:
	@-cd Test_NetcdfTraj && ./RunTest.sh $(OPT)

test.surf:
	@-cd Test_Surf && ./RunTest.sh $(OPT)

test.radgyr:
	@-cd Test_Radgyr && ./RunTest.sh $(OPT)

test.mask:
	@-cd Test_Mask && ./RunTest.sh $(OPT)

test.closest:
	@-cd Test_Closest && ./RunTest.sh $(OPT)

test.rmsd:
	@-cd Test_RMSD && ./RunTest.sh $(OPT)

test.outtraj:
	@-cd Test_Outtraj && ./RunTest.sh $(OPT)

test.datafile:
	@-cd Test_Datafile && ./RunTest.sh $(OPT)

test.avgref:
	@-cd Test_AvgRef && ./RunTest.sh $(OPT)

test.mol2:
	@-cd Test_Mol2 && ./RunTest.sh $(OPT)

test.nastruct:
	@-cd Test_NAstruct && ./RunTest.sh $(OPT)

test.pucker:
	@-cd Test_Pucker && ./RunTest.sh $(OPT)

test.atommap:
	@-cd Test_AtomMap && ./RunTest.sh $(OPT)

test.2drms:
	@-cd Test_2DRMS && ./RunTest.sh $(OPT)

test.average:
	@-cd Test_Average && ./RunTest.sh $(OPT)

test.radial:
	@-cd Test_Radial && ./RunTest.sh $(OPT)

test.hist:
	@-cd Test_Hist && ./RunTest.sh $(OPT)

test.ncrestart:
	@-cd Test_NCrestart && ./RunTest.sh $(OPT)

test.drmsd:
	@-cd Test_DRMSD && ./RunTest.sh $(OPT)

test.jcoupling:
	@-cd Test_Jcoupling && ./RunTest.sh $(OPT)

test.cluster:
	@-cd Test_Cluster && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_Sieve && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_RandomSieve && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_DataSet && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_DBSCAN && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_Kmeans && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_SymmRMSD && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_TrajWrites && ./RunTest.sh $(OPT)
	@-cd Test_Cluster_Dpeaks && ./RunTest.sh $(OPT)

test.ired:
	@-cd Test_IRED && ./RunTest.sh $(OPT)

test.molsurf:
	@-cd Test_Molsurf && ./RunTest.sh $(OPT)

test.charmm:
	@-cd Test_Charmm && ./RunTest.sh $(OPT)

test.rotdif:
	@-cd Test_Rotdif && ./RunTest.sh $(OPT)

test.corr:
	@-cd Test_Corr && ./RunTest.sh $(OPT)

test.checkstructure:
	@-cd Test_CheckStructure && ./RunTest.sh $(OPT)

test.runavg:
	@-cd Test_Runavg && ./RunTest.sh $(OPT)

test.analrunavg:
	@-cd Test_AnalysisRunAvg && ./RunTest.sh $(OPT)

test.rmsavgcorr:
	@-cd Test_RmsAvgCorr && ./RunTest.sh $(OPT)

test.contacts:
	@-cd Test_Contacts && ./RunTest.sh $(OPT)

test.principal:
	@-cd Test_Principal && ./RunTest.sh $(OPT)

test.clusterdihedral:
	@-cd Test_ClusterDihedral && ./RunTest.sh $(OPT)

test.vector:
	@-cd Test_Vector && ./RunTest.sh $(OPT)

test.matrix:
	@-cd Test_Matrix && ./RunTest.sh $(OPT)

test.atomicfluct:
	@-cd Test_AtomicFluct && ./RunTest.sh $(OPT)

test.grid:
	@-cd Test_Grid && ./RunTest.sh $(OPT)

test.unwrap:
	@-cd Test_Unwrap && ./RunTest.sh $(OPT)

test.project:
	@-cd Test_Project && ./RunTest.sh $(OPT)

test.analyzemodes:
	@-cd Test_Analyze_Modes && ./RunTest.sh $(OPT)

test.watershell:
	@-cd Test_Watershell && ./RunTest.sh $(OPT)

test.diffusion:
	@-cd Test_Diffusion && ./RunTest.sh $(OPT)

test.volmap:
	@-cd Test_Volmap && ./RunTest.sh $(OPT)

test.lie:
	@-cd Test_LIE && ./RunTest.sh $(OPT)

test.spam:
	@-cd Test_SPAM && ./RunTest.sh $(OPT)

test.randomizeions:
	@-cd Test_RandomizeIons && ./RunTest.sh $(OPT)

test.density:
	@-cd Test_Density && ./RunTest.sh $(OPT)

test.lipidorder:
	@-cd Test_LipidOrder && ./RunTest.sh $(OPT)
	@-cd Test_LipidSCD && ./RunTest.sh $(OPT)

test.pairdist:
	@-cd Test_PairDist && ./RunTest.sh $(OPT)

test.gist:
	@-cd Test_GIST && ./RunTest.sh $(OPT)
	@-cd Test_GIST_Order && ./RunTest.sh $(OPT)

test.readdata:
	@-cd Test_ReadData && ./RunTest.sh $(OPT)

test.ambpdb:
	@-cd Test_AmbPDB && ./RunTest.sh $(OPT)

test.systemVF:
	@-cd Test_systemVF && ./RunTest.sh $(OPT)

test.atomiccorr:
	@-cd Test_AtomicCorr && ./RunTest.sh $(OPT)

test.autocorr:
	@-cd Test_AutoCorr && ./RunTest.sh $(OPT)

test.box:
	@-cd Test_Box && ./RunTest.sh $(OPT)

test.atomtypeselect:
	@-cd Test_AtomTypeSelect && ./RunTest.sh $(OPT)

test.crankstat:
	@-cd Test_Crank_Stat && ./RunTest.sh $(OPT)

test.createcrd:
	@-cd Test_CreateCrd && ./RunTest.sh $(OPT)

test.autoimage:
	@-cd Test_AutoImage && ./RunTest.sh $(OPT)

test.loadcrd:
	@-cd Test_LoadCrd && ./RunTest.sh $(OPT)

test.multivector:
	@-cd Test_MultiVector && ./RunTest.sh $(OPT)

test.lifetime:
	@-cd Test_Lifetime && ./RunTest.sh $(OPT)

test.fft:
	@-cd Test_FFT && ./RunTest.sh $(OPT)

test.timecorr:
	@-cd Test_Timecorr && ./RunTest.sh $(OPT)

test.crosscorr:
	@-cd Test_CrossCorr && ./RunTest.sh $(OPT)

test.symmrmsd:
	@-cd Test_SymmRmsd && ./RunTest.sh $(OPT)

test.multidihedral:
	@-cd Test_Multidihedral && ./RunTest.sh $(OPT)

test.mremd:
	@-cd Test_Ensemble_MREMD && ./RunTest.sh $(OPT)

test.makestructure:
	@-cd Test_MakeStructure && ./RunTest.sh $(OPT)

test.rotate:
	@-cd Test_Rotate && ./RunTest.sh $(OPT)

test.integrate:
	@-cd Test_Integrate && ./RunTest.sh $(OPT)

test.nativecontacts:
	@-cd Test_NativeContacts && ./RunTest.sh $(OPT)

test.checkchirality:
	@-cd Test_CheckChiral && ./RunTest.sh $(OPT)

test.energy:
	@-cd Test_Energy && ./RunTest.sh $(OPT)

test.lessplit:
	@-cd Test_LESsplit && ./RunTest.sh $(OPT)

test.permutedih:
	@-cd Test_PermuteDihedrals && ./RunTest.sh $(OPT)

test.ccp4:
	@-cd Test_CCP4 && ./RunTest.sh $(OPT)

test.dihcovar:
	@-cd Test_DihCovar && ./RunTest.sh $(OPT)

test.curvefit:
	@-cd Test_CurveFit && ./RunTest.sh $(OPT)
	@-cd Test_Multicurve && ./RunTest.sh $(OPT)

test.datasetavg:
	@-cd Test_DataSetAverage && ./RunTest.sh $(OPT)

test.vectormath:
	@-cd Test_VectorMath && ./RunTest.sh $(OPT)

test.esander:
	@-cd Test_Esander && ./RunTest.sh $(OPT)
	@-cd Test_Esander_Short && ./RunTest.sh $(OPT)

test.temperature:
	@-cd Test_Temperature && ./RunTest.sh $(OPT)

test.spline:
	@-cd Test_Spline && ./RunTest.sh $(OPT)

test.kde:
	@-cd Test_KDE && ./RunTest.sh $(OPT)

test.wavelet:
	@-cd Test_Wavelet && ./RunTest.sh $(OPT)

test.cif:
	@-cd Test_CIF && ./RunTest.sh $(OPT)

test.pdb:
	@-cd Test_PDB && ./RunTest.sh $(OPT)

test.align:
	@-cd Test_Align && ./RunTest.sh $(OPT)

test.velocityautocorr:
	@-cd Test_VelocityAutoCorr && ./RunTest.sh $(OPT)

test.velfrc:
	@-cd Test_VelFrc && ./RunTest.sh $(OPT)

test.rotatedihedral:
	@-cd Test_RotateDihedral && ./RunTest.sh $(OPT)

test.dataset:
	@-cd Test_DataSetCmd && ./RunTest.sh $(OPT)

test.remap:
	@-cd Test_Remap && ./RunTest.sh $(OPT)

test.comparetop:
	@-cd Test_CompareTop && ./RunTest.sh $(OPT)

test.filter:
	@-cd Test_Filter && ./RunTest.sh $(OPT)

test.combinecrd:
	@-cd Test_CombineCrd && ./RunTest.sh $(OPT)

test.topinfo:
	@-cd Test_TopInfo && ./RunTest.sh $(OPT)

test.gromacs:
	@-cd Test_GromacsTrr && ./RunTest.sh $(OPT)
	@-cd Test_GromacsXtc && ./RunTest.sh $(OPT)

test.ewald:
	@-cd Test_Ewald && ./RunTest.sh $(OPT)

test.ti:
	@-cd Test_TI && ./RunTest.sh $(OPT)

test.fiximagedbonds:
	@-cd Test_FixImagedBonds && ./RunTest.sh $(OPT)

test.distance:
	@-cd Test_Distance && ./RunTest.sh $(OPT)

test.volume:
	@-cd Test_Volume && ./RunTest.sh $(OPT)

test.setvelocity:
	@-cd Test_SetVelocity && ./RunTest.sh $(OPT)

test.remlog:
	@-cd Test_Remlog && ./RunTest.sh $(OPT)

test.areapermol:
	@-cd Test_AreaPerMol && ./RunTest.sh $(OPT)

test.distbasedmask:
	@-cd Test_DistBasedMask && ./RunTest.sh $(OPT)

# Every test target should go here.
COMPLETETESTS=test.general \
              test.strip \
              test.brokentraj \
              test.trajinoffset \
              test.remdnetcdf \
              test.framerange \
              test.perres \
              test.dssp \
              test.center \
              test.hbond \
              test.ncrestart \
              test.image \
              test.surf \
              test.radgyr \
              test.mask \
              test.jcoupling \
              test.diffusion \
              test.outtraj \
              test.charmm \
              test.rotdif \
              test.netcdftraj \
              test.runavg \
              test.rmsd \
              test.2drms \
              test.mol2 \
              test.nastruct \
              test.average \
              test.pucker \
              test.atommap \
              test.datafile \
              test.radial \
              test.hist \
              test.closest \
              test.drmsd \
              test.cluster \
              test.ired \
              test.randomizeions \
              test.molsurf \
              test.corr \
              test.checkstructure \
              test.rmsavgcorr \
              test.contacts \
              test.principal \
              test.clusterdihedral \
              test.vector \
              test.matrix \
              test.atomicfluct \
              test.grid \
              test.unwrap \
              test.project \
              test.analyzemodes \
              test.watershell \
              test.analrunavg \
              test.volmap \
              test.spam \
              test.lie \
              test.density \
              test.lipidorder \
              test.pairdist \
              test.readdata \
              test.gist \
              test.ambpdb \
              test.systemVF \
              test.atomiccorr \
              test.autocorr \
              test.box \
              test.atomtypeselect \
              test.avgref \
              test.crankstat \
              test.createcrd \
              test.autoimage \
              test.loadcrd \
              test.multivector \
              test.lifetime \
              test.fft \
              test.timecorr \
              test.crosscorr \
              test.symmrmsd \
              test.multidihedral \
              test.mremd \
              test.makestructure \
              test.rotate \
              test.integrate \
              test.nativecontacts \
              test.checkchirality \
              test.energy \
              test.lessplit \
              test.permutedih \
              test.ccp4 \
              test.dihcovar \
              test.curvefit \
              test.datasetavg \
              test.vectormath \
              test.esander \
              test.temperature \
              test.spline \
              test.kde \
              test.wavelet \
              test.cif \
              test.pdb \
              test.align \
              test.velocityautocorr \
              test.velfrc \
              test.rotatedihedral \
              test.dataset \
              test.remap \
              test.comparetop \
              test.filter \
              test.combinecrd \
              test.topinfo \
              test.gromacs \
              test.ewald \
              test.ti \
              test.fiximagedbonds \
              test.distance \
              test.volume \
              test.setvelocity \
              test.remlog \
              test.areapermol \
              test.distbasedmask

test.all:
	$(MAKE) test.complete summary

test.mpi:
	$(MAKE) test.complete summary OPT="$(OPT) mpi"

test.openmp:
	$(MAKE) test.complete summary OPT="$(OPT) openmp"

test.cuda:
	$(MAKE) test.complete summary OPT="$(OPT) cuda"

test.libcpptraj:
	@-cd Test_Libcpptraj && ./RunTest.sh $(OPT)

test.complete: CpptrajTest.sh MasterTest.sh
	@./CpptrajTest.sh --target test.cpptraj $(OPT)

test.cpptraj: $(COMPLETETESTS)

test.showerrors:
	$(MAKE) test.complete summary OPT="$(OPT) showerrors long"

summary: CpptrajTest.sh MasterTest.sh
	@./CpptrajTest.sh summary $(OPT)

test.vg:
	$(MAKE) test.complete summary OPT="$(OPT) vg"

clean:
	@echo "Cleaning all test directories:"
	$(MAKE) test.complete test.libcpptraj OPT=clean

cleanprofile:
	@echo "Removing GNU profile files:"
	(find . -name 'libstdcxx-profile.*' | \
        while read dif ;\
        do \
            rm -f $$dif ;\
        done ;\
        )

