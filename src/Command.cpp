#include <cctype> // isalnum
#include <cstdarg>
#include <algorithm> // std::sort()
#include "Command.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h" // ProcessInput()
#include "CmdInput.h"     // ProcessInput()
#include "RPNcalc.h"
#include "Deprecated.h"
#include "Control.h"
// ----- GENERAL ---------------------------------------------------------------
#include "Exec_Analyze.h"
#include "Exec_Calc.h"
#include "Exec_ClusterMap.h"
#include "Exec_Commands.h"
#include "Exec_CreateSet.h"
#include "Exec_DataFile.h"
#include "Exec_DataFilter.h"
#include "Exec_DataSetCmd.h"
#include "Exec_GenerateAmberRst.h"
#include "Exec_Help.h"
#include "Exec_Precision.h"
#include "Exec_PrintData.h"
#include "Exec_ReadData.h"
#include "Exec_ReadEnsembleData.h"
#include "Exec_ReadInput.h"
#include "Exec_RunAnalysis.h"
#include "Exec_SortEnsembleData.h"
#include "Exec_SequenceAlign.h"
#include "Exec_ViewRst.h"
// ----- SYSTEM ----------------------------------------------------------------
#include "Exec_System.h"
// ----- COORDS ----------------------------------------------------------------
#include "Exec_CombineCoords.h"
#include "Exec_CrdAction.h"
#include "Exec_CrdOut.h"
#include "Exec_LoadCrd.h"
#include "Exec_LoadTraj.h"
#include "Exec_PermuteDihedrals.h"
#include "Exec_RotateDihedral.h"
// ----- TRAJECTORY ------------------------------------------------------------
#include "Exec_Traj.h"
// ----- TOPOLOGY --------------------------------------------------------------
#include "Exec_Change.h"
#include "Exec_CompareTop.h"
#include "Exec_ParmBox.h"
#include "Exec_ParmSolvent.h"
#include "Exec_ParmStrip.h"
#include "Exec_ParmWrite.h"
#include "Exec_ScaleDihedralK.h"
#include "Exec_Top.h"
#include "Exec_UpdateParameters.h"
// ----- ACTION ----------------------------------------------------------------
#include "Action_Angle.h"
#include "Action_Distance.h"
#include "Action_Rmsd.h"
#include "Action_Dihedral.h"
#include "Action_AtomMap.h"
#include "Action_Strip.h"
#include "Action_Unstrip.h"
#include "Action_DSSP.h"
#include "Action_Center.h"
#include "Action_Image.h"
#include "Action_Surf.h"
#include "Action_Radgyr.h"
#include "Action_Mask.h"
#include "Action_Closest.h"
#include "Action_NAstruct.h"
#include "Action_Pucker.h"
#include "Action_Outtraj.h"
#include "Action_Average.h"
#include "Action_Radial.h"
#include "Action_DistRmsd.h"
#include "Action_Jcoupling.h"
#include "Action_Pairwise.h"
#include "Action_Molsurf.h"
#include "Action_CheckStructure.h"
#include "Action_RunningAvg.h"
#include "Action_AtomicFluct.h"
#include "Action_Watershell.h"
#include "Action_Contacts.h"
#include "Action_Vector.h"
#include "Action_Principal.h"
#include "Action_Matrix.h"
#include "Action_LIE.h"
#include "Action_Grid.h"
#include "Action_GridFreeEnergy.h"
#include "Action_Dipole.h"
#include "Action_Projection.h"
#include "Action_ClusterDihedral.h"
#include "Action_Unwrap.h"
#include "Action_Diffusion.h"
#include "Action_DNAionTracker.h"
#include "Action_Scale.h"
#include "Action_RandomizeIons.h"
#include "Action_AutoImage.h"
#include "Action_STFC_Diffusion.h"
#include "Action_AtomicCorr.h"
#include "Action_Bounds.h"
#include "Action_Rotate.h"
#include "Action_Translate.h"
#include "Action_Box.h"
#include "Action_CreateCrd.h"
#include "Action_MultiDihedral.h"
#include "Action_MakeStructure.h"
#include "Action_SymmetricRmsd.h"
#include "Action_Volmap.h"
#include "Action_Spam.h"
#include "Action_Temperature.h"
#include "Action_GIST.h"
#include "Action_CreateReservoir.h"
#include "Action_Density.h"
#include "Action_PairDist.h"
#include "Action_OrderParameter.h"
#include "Action_FixAtomOrder.h"
#include "Action_NMRrst.h"
#include "Action_FilterByData.h"
#include "Action_LESsplit.h"
#include "Action_NativeContacts.h"
#include "Action_VelocityAutoCorr.h"
#include "Action_SetVelocity.h"
#include "Action_MultiVector.h"
#include "Action_MinImage.h"
#include "Action_ReplicateCell.h"
#include "Action_AreaPerMol.h"
#include "Action_Energy.h"
#include "Action_Esander.h"
#include "Action_CheckChirality.h"
#include "Action_Channel.h" // EXPERIMENTAL
#include "Action_Volume.h"
#include "Action_Align.h"
#include "Action_Remap.h"
#include "Action_HydrogenBond.h"
#include "Action_FixImagedBonds.h"
#include "Action_LipidOrder.h"
#include "Action_InfraredSpectrum.h"
// ----- ANALYSIS --------------------------------------------------------------
#include "Analysis_Hist.h"
#include "Analysis_Corr.h"
#include "Analysis_Matrix.h"
#include "Analysis_Timecorr.h"
#include "Analysis_IRED.h"
#include "Analysis_Modes.h"
#include "Analysis_CrankShaft.h"
#include "Analysis_Statistics.h"
#include "Analysis_CrossCorr.h"
#include "Analysis_AutoCorr.h"
#include "Analysis_Lifetime.h"
#include "Analysis_FFT.h"
#include "Analysis_CrdFluct.h"
#include "Analysis_RmsAvgCorr.h"
#include "Analysis_Rms2d.h"
#include "Analysis_Clustering.h"
#include "Analysis_RunningAvg.h"
#include "Analysis_MeltCurve.h"
#include "Analysis_Overlap.h"
#include "Analysis_AmdBias.h"
#include "Analysis_RemLog.h"
#include "Analysis_Integrate.h"
#include "Analysis_Spline.h"
#include "Analysis_Average.h"
#include "Analysis_KDE.h"
#include "Analysis_MultiHist.h"
#include "Analysis_Divergence.h"
#include "Analysis_VectorMath.h"
#include "Analysis_Regression.h"
#include "Analysis_LowestCurve.h"
#include "Analysis_CurveFit.h"
#include "Analysis_PhiPsi.h"
#include "Analysis_Rotdif.h"
#include "Analysis_Wavelet.h"
#include "Analysis_State.h"
#include "Analysis_Multicurve.h"
#include "Analysis_TI.h"
#include "Analysis_ConstantPHStats.h"

CmdList Command::commands_ = CmdList();

const Cmd Command::EMPTY_ = Cmd();

Command::Carray Command::names_ = Command::Carray();

Command::CtlArray Command::control_ = Command::CtlArray();

int Command::ctlidx_ = -1;

VariableArray Command::CurrentVars_ = VariableArray();

/** Initialize all commands. Should only be called once as program starts. */
void Command::Init() {
  // GENERAL
  Command::AddCmd( new Exec_ActiveRef(),       Cmd::EXE, 1, "activeref" );
  Command::AddCmd( new Exec_Analyze(),         Cmd::EXE, 1, "analyze" ); // HIDDEN
  Command::AddCmd( new Exec_Calc(),            Cmd::EXE, 1, "calc" );
  Command::AddCmd( new Exec_Clear(),           Cmd::EXE, 1, "clear" );
  Command::AddCmd( new Exec_ClusterMap(),      Cmd::EXE, 1, "clustermap" ); // HIDDEN
  Command::AddCmd( new Exec_CreateDataFile(),  Cmd::EXE, 1, "create" );
  Command::AddCmd( new Exec_CreateSet(),       Cmd::EXE, 1, "createset" );
  Command::AddCmd( new Exec_DataFileCmd(),     Cmd::EXE, 1, "datafile" );
  Command::AddCmd( new Exec_DataFilter(),      Cmd::EXE, 1, "datafilter" );
  Command::AddCmd( new Exec_DataSetCmd(),      Cmd::EXE, 1, "dataset" );
  Command::AddCmd( new Exec_EnsFileExt(),      Cmd::EXE, 1, "ensextension" );
  Command::AddCmd( new Exec_GenerateAmberRst(),Cmd::EXE, 1, "rst" );
  Command::AddCmd( new Exec_Help(),            Cmd::EXE, 1, "help" );
  Command::AddCmd( new Exec_ListAll(),         Cmd::EXE, 1, "list" );
  Command::AddCmd( new Exec_NoExitOnError(),   Cmd::EXE, 1, "noexitonerror" );
  Command::AddCmd( new Exec_NoProgress(),      Cmd::EXE, 1, "noprogress" );
  Command::AddCmd( new Exec_Precision(),       Cmd::EXE, 1, "precision" );
  Command::AddCmd( new Exec_PrintData(),       Cmd::EXE, 1, "printdata" );
  Command::AddCmd( new Exec_QuietBlocks(),     Cmd::EXE, 1, "quietblocks" );
  Command::AddCmd( new Exec_Quit(),            Cmd::EXE, 2, "exit", "quit" );
  Command::AddCmd( new Exec_ReadData(),        Cmd::EXE, 1, "readdata" );
  Command::AddCmd( new Exec_ReadEnsembleData(),Cmd::EXE, 1, "readensembledata" );
  Command::AddCmd( new Exec_ReadInput(),       Cmd::EXE, 1, "readinput" );
  Command::AddCmd( new Exec_RemoveData(),      Cmd::EXE, 1, "removedata" );
  Command::AddCmd( new Exec_Run(),             Cmd::EXE, 2, "go", "run" );
  Command::AddCmd( new Exec_RunAnalysis(),     Cmd::EXE, 1, "runanalysis" );
  Command::AddCmd( new Exec_SelectAtoms(),     Cmd::EXE, 1, "select" );
  Command::AddCmd( new Exec_SelectDS(),        Cmd::EXE, 1, "selectds" );
  Command::AddCmd( new Exec_SetListDebug(),    Cmd::EXE, 2, "debug", "prnlev" );
  Command::AddCmd( new Exec_SilenceActions(),  Cmd::EXE, 1, "silenceactions" );
  Command::AddCmd( new Exec_SequenceAlign(),   Cmd::EXE, 1, "sequencealign" );
  Command::AddCmd( new Exec_SortEnsembleData(),Cmd::EXE, 1, "sortensembledata" );
  Command::AddCmd( new Exec_WriteDataFile(),   Cmd::EXE, 2, "write", "writedata" );
  Command::AddCmd( new Exec_ViewRst(),         Cmd::EXE, 1, "viewrst" ); // HIDDEN
# ifdef MPI
  Command::AddCmd( new Exec_ForceParaEnsemble(), Cmd::EXE, 1, "forceparaensemble" );
# endif
  // SYSTEM
  Command::AddCmd( new Exec_System(), Cmd::EXE, 6, "gnuplot", "head", "less", "ls", "pwd", "xmgrace" );
  // COORDS
  Command::AddCmd( new Exec_CombineCoords(),    Cmd::EXE, 1, "combinecrd" ); 
  Command::AddCmd( new Exec_CrdAction(),        Cmd::EXE, 1, "crdaction" );
  Command::AddCmd( new Exec_CrdOut(),           Cmd::EXE, 1, "crdout" );
  Command::AddCmd( new Exec_LoadCrd(),          Cmd::EXE, 1, "loadcrd" );
  Command::AddCmd( new Exec_LoadTraj(),         Cmd::EXE, 1, "loadtraj" );
  Command::AddCmd( new Exec_PermuteDihedrals(), Cmd::EXE, 1, "permutedihedrals" );
  Command::AddCmd( new Exec_RotateDihedral(),   Cmd::EXE, 1, "rotatedihedral" );
  // TRAJECTORY
  Command::AddCmd( new Exec_Ensemble(),     Cmd::EXE, 1, "ensemble" );
  Command::AddCmd( new Exec_EnsembleSize(), Cmd::EXE, 1, "ensemblesize" );
  Command::AddCmd( new Exec_Reference(),    Cmd::EXE, 1, "reference" );
  Command::AddCmd( new Exec_Trajin(),       Cmd::EXE, 1, "trajin" );
  Command::AddCmd( new Exec_Trajout(),      Cmd::EXE, 1, "trajout" );
  // TOPOLOGY COMMANDS
  Command::AddCmd( new Exec_AngleInfo(),     Cmd::EXE, 3, "angles", "angleinfo", "printangles" );
  Command::AddCmd( new Exec_AtomInfo(),      Cmd::EXE, 3, "atoms", "atominfo", "printatoms" );
  Command::AddCmd( new Exec_BondInfo(),      Cmd::EXE, 3, "bonds", "bondinfo", "printbonds" );
  Command::AddCmd( new Exec_Change(),        Cmd::EXE, 1, "change" );
  Command::AddCmd( new Exec_ChargeInfo(),    Cmd::EXE, 1, "charge" );
  Command::AddCmd( new Exec_CompareTop(),    Cmd::EXE, 1, "comparetop" );
  Command::AddCmd( new Exec_DihedralInfo(),Cmd::EXE, 3,"dihedrals","dihedralinfo","printdihedrals");
  Command::AddCmd( new Exec_ImproperInfo(),Cmd::EXE, 3,"impropers","improperinfo","printimpropers");
  Command::AddCmd( new Exec_MassInfo(),      Cmd::EXE, 1, "mass" );
  Command::AddCmd( new Exec_MolInfo(),       Cmd::EXE, 1, "molinfo" );
  Command::AddCmd( new Exec_LoadParm(),      Cmd::EXE, 1, "parm" );
  Command::AddCmd( new Exec_ParmBox(),       Cmd::EXE, 1, "parmbox" );
  Command::AddCmd( new Exec_ParmInfo(),      Cmd::EXE, 1, "parminfo" );
  Command::AddCmd( new Exec_ParmSolvent(),   Cmd::EXE, 1, "solvent" );
  Command::AddCmd( new Exec_ParmStrip(),     Cmd::EXE, 1, "parmstrip" );
  Command::AddCmd( new Exec_ParmWrite(),     Cmd::EXE, 1, "parmwrite" );
  Command::AddCmd( new Exec_ResInfo(),       Cmd::EXE, 1, "resinfo" );
  Command::AddCmd( new Exec_ScaleDihedralK(),Cmd::EXE, 1, "scaledihedralk" );
  Command::AddCmd( new Exec_UBInfo(),        Cmd::EXE, 2, "ubinfo", "printub" );
  Command::AddCmd( new Exec_UpdateParameters(), Cmd::EXE, 1, "updateparameters"); // HIDDEN
  // ACTION
  Command::AddCmd( new Action_Align(),         Cmd::ACT, 1, "align" );
  Command::AddCmd( new Action_Angle(),         Cmd::ACT, 1, "angle" );
  Command::AddCmd( new Action_AreaPerMol(),    Cmd::ACT, 1, "areapermol" );
  Command::AddCmd( new Action_AtomicCorr(),    Cmd::ACT, 1, "atomiccorr" );
  Command::AddCmd( new Action_AtomicFluct(),   Cmd::ACT, 2, "atomicfluct", "rmsf" );
  Command::AddCmd( new Action_AtomMap(),       Cmd::ACT, 1, "atommap" );
  Command::AddCmd( new Action_AutoImage(),     Cmd::ACT, 1, "autoimage" );
  Command::AddCmd( new Action_Average(),       Cmd::ACT, 1, "average" );
  Command::AddCmd( new Action_Bounds(),        Cmd::ACT, 1, "bounds" );
  Command::AddCmd( new Action_Box(),           Cmd::ACT, 1, "box" );
  Command::AddCmd( new Action_Center(),        Cmd::ACT, 1, "center" );
  Command::AddCmd( new Action_Channel(),       Cmd::ACT, 1, "channel" ); // HIDDEN
  Command::AddCmd( new Action_CheckStructure(),Cmd::ACT, 3,"check","checkoverlap","checkstructure");
  Command::AddCmd( new Action_CheckChirality(),Cmd::ACT, 1, "checkchirality" );
  Command::AddCmd( new Action_Closest(),       Cmd::ACT, 2, "closest", "closestwaters" );
  Command::AddCmd( new Action_ClusterDihedral(),Cmd::ACT,1, "clusterdihedral" );
  Command::AddCmd( new Action_Contacts(),      Cmd::ACT, 1, "contacts" );
  Command::AddCmd( new Action_CreateCrd(),     Cmd::ACT, 1, "createcrd" );
  Command::AddCmd( new Action_CreateReservoir(),Cmd::ACT,1, "createreservoir" );
  Command::AddCmd( new Action_Density(),       Cmd::ACT, 1, "density" );
  Command::AddCmd( new Action_Diffusion(),     Cmd::ACT, 1, "diffusion" );
  Command::AddCmd( new Action_Dihedral(),      Cmd::ACT, 1, "dihedral" );
  Command::AddCmd( new Action_Dipole(),        Cmd::ACT, 1, "dipole" );
  Command::AddCmd( new Action_Distance(),      Cmd::ACT, 1, "distance" );
  Command::AddCmd( new Action_DNAionTracker(), Cmd::ACT, 1, "dnaiontracker" ); // HIDDEN
  Command::AddCmd( new Action_DistRmsd(),      Cmd::ACT, 2, "drms", "drmsd" );
  Command::AddCmd( new Action_DSSP(),          Cmd::ACT, 2, "dssp", "secstruct" );
  Command::AddCmd( new Action_Energy(),        Cmd::ACT, 1, "energy" );
  Command::AddCmd( new Action_Esander(),       Cmd::ACT, 1, "esander" );
  Command::AddCmd( new Action_FilterByData(),  Cmd::ACT, 1, "filter" );
  Command::AddCmd( new Action_FixAtomOrder(),  Cmd::ACT, 1, "fixatomorder" );
  Command::AddCmd( new Action_FixImagedBonds(),Cmd::ACT, 1, "fiximagedbonds" );
  Command::AddCmd( new Action_GIST(),          Cmd::ACT, 1, "gist" );
  Command::AddCmd( new Action_GridFreeEnergy(),Cmd::ACT, 1, "gfe" ); // HIDDEN
  Command::AddCmd( new Action_Grid(),          Cmd::ACT, 1, "grid" );
  Command::AddCmd( new Action_HydrogenBond(),  Cmd::ACT, 1, "hbond" );
  Command::AddCmd( new Action_Image(),         Cmd::ACT, 1, "image" );
  Command::AddCmd( new Action_InfraredSpectrum(),Cmd::ACT,2,"irspec","infraredspec");
  Command::AddCmd( new Action_Jcoupling(),     Cmd::ACT, 1, "jcoupling" );
  Command::AddCmd( new Action_LESsplit(),      Cmd::ACT, 1, "lessplit" );
  Command::AddCmd( new Action_LIE(),           Cmd::ACT, 1, "lie" );
  Command::AddCmd( new Action_OrderParameter(),Cmd::ACT, 1, "lipidorder" );
  Command::AddCmd( new Action_LipidOrder(),    Cmd::ACT, 1, "lipidscd" );
  Command::AddCmd( new Action_MakeStructure(), Cmd::ACT, 1, "makestructure" );
  Command::AddCmd( new Action_Mask(),          Cmd::ACT, 1, "mask" );
  Command::AddCmd( new Action_Matrix(),        Cmd::ACT, 1, "matrix" );
  Command::AddCmd( new Action_MinImage(),      Cmd::ACT, 1, "minimage" );
  Command::AddCmd( new Action_Molsurf(),       Cmd::ACT, 1, "molsurf" );
  Command::AddCmd( new Action_MultiDihedral(), Cmd::ACT, 1, "multidihedral" );
  Command::AddCmd( new Action_MultiVector(),   Cmd::ACT, 1, "multivector" );
  Command::AddCmd( new Action_NAstruct(),      Cmd::ACT, 1, "nastruct" );
  Command::AddCmd( new Action_NativeContacts(),Cmd::ACT, 1, "nativecontacts" );
  Command::AddCmd( new Action_NMRrst(),        Cmd::ACT, 1, "nmrrst" ); // HIDDEN
  Command::AddCmd( new Action_Outtraj(),       Cmd::ACT, 1, "outtraj" );
  Command::AddCmd( new Action_PairDist(),      Cmd::ACT, 1, "pairdist" );
  Command::AddCmd( new Action_Pairwise(),      Cmd::ACT, 1, "pairwise" );
  Command::AddCmd( new Action_Principal(),     Cmd::ACT, 1, "principal" );
  Command::AddCmd( new Action_Projection(),    Cmd::ACT, 1, "projection" );
  Command::AddCmd( new Action_Pucker(),        Cmd::ACT, 1, "pucker" );
  Command::AddCmd( new Action_Radgyr(),        Cmd::ACT, 2, "radgyr", "rog" );
  Command::AddCmd( new Action_Radial(),        Cmd::ACT, 2, "radial", "rdf" );
  Command::AddCmd( new Action_RandomizeIons(), Cmd::ACT, 1, "randomizeions" );
  Command::AddCmd( new Action_Remap(),         Cmd::ACT, 1, "remap" );
  Command::AddCmd( new Action_ReplicateCell(), Cmd::ACT, 1, "replicatecell" );
  Command::AddCmd( new Action_Rmsd(),          Cmd::ACT, 2, "rms", "rmsd" );
  Command::AddCmd( new Action_Rotate(),        Cmd::ACT, 1, "rotate" );
  Command::AddCmd( new Action_RunningAvg(),    Cmd::ACT, 2, "runavg", "runningaverage" );
  Command::AddCmd( new Action_Scale(),         Cmd::ACT, 1, "scale" );
  Command::AddCmd( new Action_SetVelocity(),   Cmd::ACT, 1, "setvelocity" );
  Command::AddCmd( new Action_Spam(),          Cmd::ACT, 1, "spam" ); // HIDDEN
  Command::AddCmd( new Action_STFC_Diffusion(),Cmd::ACT, 1, "stfcdiffusion" );
  Command::AddCmd( new Action_Strip(),         Cmd::ACT, 1, "strip" );
  Command::AddCmd( new Action_Surf(),          Cmd::ACT, 1, "surf" );
  Command::AddCmd( new Action_SymmetricRmsd(), Cmd::ACT, 1, "symmrmsd" );
  Command::AddCmd( new Action_Temperature(),   Cmd::ACT, 1, "temperature" );
  Command::AddCmd( new Action_Translate(),     Cmd::ACT, 2, "trans", "translate" );
  Command::AddCmd( new Action_Unstrip(),       Cmd::ACT, 1, "unstrip" );
  Command::AddCmd( new Action_Unwrap(),        Cmd::ACT, 1, "unwrap" );
  Command::AddCmd( new Action_Vector(),        Cmd::ACT, 1, "vector" );
  Command::AddCmd( new Action_VelocityAutoCorr(),Cmd::ACT,1,"velocityautocorr" );
  Command::AddCmd( new Action_Volmap(),        Cmd::ACT, 1, "volmap" );
  Command::AddCmd( new Action_Volume(),        Cmd::ACT, 1, "volume" );
  Command::AddCmd( new Action_Watershell(),    Cmd::ACT, 1, "watershell" );
  // ANALYSIS
  Command::AddCmd( new Analysis_AmdBias(),     Cmd::ANA, 1, "amdbias" ); // HIDDEN
  Command::AddCmd( new Analysis_AutoCorr(),    Cmd::ANA, 1, "autocorr" );
  Command::AddCmd( new Analysis_Average(),     Cmd::ANA, 1, "avg" );
  Command::AddCmd( new Analysis_State(),       Cmd::ANA, 1, "calcstate" );
  Command::AddCmd( new Analysis_Clustering(),  Cmd::ANA, 1, "cluster" );
  Command::AddCmd( new Analysis_Corr(),        Cmd::ANA, 2, "corr", "correlationcoe" );
  Command::AddCmd( new Analysis_ConstantPHStats,Cmd::ANA,1, "cphstats" );
  Command::AddCmd( new Analysis_CrankShaft(),  Cmd::ANA, 2, "crank", "crankshaft" );
  Command::AddCmd( new Analysis_CrdFluct(),    Cmd::ANA, 1, "crdfluct" );
  Command::AddCmd( new Analysis_CrossCorr(),   Cmd::ANA, 1, "crosscorr" );
  Command::AddCmd( new Analysis_CurveFit(),    Cmd::ANA, 1, "curvefit" );
  Command::AddCmd( new Analysis_Matrix(),      Cmd::ANA, 2, "diagmatrix", "matrix" );
  Command::AddCmd( new Analysis_Divergence(),  Cmd::ANA, 1, "divergence" );
  Command::AddCmd( new Analysis_FFT(),         Cmd::ANA, 1, "fft" );
  Command::AddCmd( new Analysis_Hist(),        Cmd::ANA, 2, "hist", "histogram" );
  Command::AddCmd( new Analysis_Integrate(),   Cmd::ANA, 1, "integrate" );
  Command::AddCmd( new Analysis_IRED(),        Cmd::ANA, 1, "ired" );
  Command::AddCmd( new Analysis_KDE(),         Cmd::ANA, 1, "kde" );
  Command::AddCmd( new Analysis_Lifetime(),    Cmd::ANA, 1, "lifetime" );
  Command::AddCmd( new Analysis_LowestCurve(), Cmd::ANA, 1, "lowestcurve" );
  Command::AddCmd( new Analysis_MeltCurve(),   Cmd::ANA, 1, "meltcurve" );
  Command::AddCmd( new Analysis_Modes(),       Cmd::ANA, 1, "modes" );
  Command::AddCmd( new Analysis_Multicurve(),  Cmd::ANA, 1, "multicurve" );
  Command::AddCmd( new Analysis_MultiHist(),   Cmd::ANA, 1, "multihist" );
  Command::AddCmd( new Analysis_Overlap(),     Cmd::ANA, 1, "overlap" ); // HIDDEN
  Command::AddCmd( new Analysis_PhiPsi(),      Cmd::ANA, 1, "phipsi" );
  Command::AddCmd( new Analysis_Regression(),  Cmd::ANA, 1, "regress" );
  Command::AddCmd( new Analysis_RemLog(),      Cmd::ANA, 1, "remlog" );
  Command::AddCmd( new Analysis_Rms2d(),       Cmd::ANA, 2, "2drms", "rms2d" );
  Command::AddCmd( new Analysis_RmsAvgCorr(),  Cmd::ANA, 1, "rmsavgcorr" );
  Command::AddCmd( new Analysis_Rotdif(),      Cmd::ANA, 1, "rotdif" );
  Command::AddCmd( new Analysis_RunningAvg(),  Cmd::ANA, 1, "runningavg" );
  Command::AddCmd( new Analysis_Spline(),      Cmd::ANA, 1, "spline" );
  Command::AddCmd( new Analysis_Statistics(),  Cmd::ANA, 2, "stat", "statistics" );
  Command::AddCmd( new Analysis_TI(),          Cmd::ANA, 1, "ti" );
  Command::AddCmd( new Analysis_Timecorr(),    Cmd::ANA, 1, "timecorr" );
  Command::AddCmd( new Analysis_VectorMath(),  Cmd::ANA, 1, "vectormath" );
  Command::AddCmd( new Analysis_Wavelet(),     Cmd::ANA, 1, "wavelet" );
  // CONTROL STRUCTURES
  Command::AddCmd( new ControlBlock_For(),     Cmd::BLK, 1, "for" );
  Command::AddCmd( new Control_Set(),          Cmd::CTL, 1, "set" );
  Command::AddCmd( new Control_Show(),         Cmd::CTL, 1, "show" );
  // DEPRECATED COMMANDS
  Command::AddCmd( new Deprecated_AvgCoord(),    Cmd::DEP, 1, "avgcoord" );
  Command::AddCmd( new Deprecated_DihScan(),     Cmd::DEP, 1, "dihedralscan" );
  Command::AddCmd( new Deprecated_Hbond(),       Cmd::DEP, 2, "acceptor", "donor" );
  Command::AddCmd( new Deprecated_MinDist(),     Cmd::DEP, 2, "mindist", "maxdist" );
  Command::AddCmd( new Deprecated_ParmBondInfo(),Cmd::DEP, 1, "parmbondinfo" );
  Command::AddCmd( new Deprecated_ParmMolInfo(), Cmd::DEP, 1, "parmmolinfo" );
  Command::AddCmd( new Deprecated_ParmResInfo(), Cmd::DEP, 1, "parmresinfo" );
  Command::AddCmd( new Deprecated_TopSearch(),   Cmd::DEP, 4, "bondsearch", "molsearch", "nobondsearch", "nomolsearch" );

  // Add null ptr to indicate end of command key addresses for ReadLine
  names_.push_back( 0 );
}

/** Clear any existing control blocks. */
void Command::ClearControlBlocks() {
  for (CtlArray::iterator it = control_.begin(); it != control_.end(); ++it)
    delete *it;
  control_.clear();
  ctlidx_ = -1;
}

/** Free all commands. Should only be called just before program exit. Also
  * remove any remaining control blocks.
  */
void Command::Free() {
  commands_.Clear();
  ClearControlBlocks();
}

/** \param oIn Pointer to DispatchObject to add as command.
  * \param dIn Command destination
  * \param nKeys Number of command keywords associated with this command.
  * The remaining arguments are the nKeys command keywords.
  */
void Command::AddCmd(DispatchObject* oIn, Cmd::DestType dIn, int nKeys, ...) {
  Cmd::Sarray keys;
  va_list args;
  va_start(args, nKeys);
  for (int nk = 0; nk < nKeys; nk++) {
    char* key = va_arg(args, char*);
    keys.push_back( std::string(key) );
  }
  va_end(args);
  commands_.Add( Cmd(oIn, keys, dIn) );
  // Store memory addresses of command keys for ReadLine
  for (Cmd::key_iterator key = commands_.Back().keysBegin();
                         key != commands_.Back().keysEnd(); ++key)
    names_.push_back( key->c_str() );
}

/** Search Commands list for command with given keyword and object type. */
Cmd const& Command::SearchTokenType(DispatchObject::Otype catIn, const char* keyIn,
                                    bool silent)
{
  for (CmdList::const_iterator cmd = commands_.begin(); cmd != commands_.end(); ++cmd)
  {
    if (catIn != cmd->Obj().Type()) continue;
    if (cmd->KeyMatches(keyIn)) return *cmd;
  }
  if (!silent) mprinterr("'%s': Command not found.\n", keyIn);
  return EMPTY_;
}

/** Search Commands list for command with given keyword and object type. */
Cmd const& Command::SearchTokenType(DispatchObject::Otype catIn, const char* keyIn) {
  return SearchTokenType(catIn, keyIn, false);
}

/** Search the Commands list for given command.
  * \return the token if found, 0 if not.
  */
Cmd const& Command::SearchToken(ArgList& argIn) {
  // Search for command.
  for (CmdList::const_iterator cmd = commands_.begin(); cmd != commands_.end(); ++cmd)
  {
    if ( cmd->KeyMatches( argIn.Command() ) )
      return *cmd;
  }
  //mprinterr("'%s': Command not found.\n", argIn.Command());
  return EMPTY_;
}

/** First list the command category, then the commands for that category
  * in alphabetical order. Should not be called with NONE, HIDDEN, or
  * DEPRECATED.
  */
void Command::ListCommandsForType(DispatchObject::Otype typeIn) {
  std::vector< std::string > command_keys;
  mprintf("%s Commands:\n", DispatchObject::ObjKeyword(typeIn));
  for (CmdList::const_iterator cmd = commands_.begin(); cmd != commands_.end(); ++cmd)
  {
    if (cmd->Obj().Type() == typeIn)
      for (Cmd::key_iterator key = cmd->keysBegin(); key != cmd->keysEnd(); ++key)
        command_keys.push_back( *key );
  }
  std::sort( command_keys.begin(), command_keys.end() );
  std::string Line("        ");
  for (std::vector< std::string >::const_iterator key = command_keys.begin();
                                                  key != command_keys.end(); ++key)
  {
    if ( Line.size() + key->size() + 1 > 80 ) {
      mprintf("%s\n", Line.c_str());
      Line.assign("        ");
    }
    Line.append( *key + " " );
  }
  if (!Line.empty()) // TODO is it ever empty?
    mprintf("%s\n", Line.c_str());
}
    
/** List all commands of the given type, or all commands if type
  * is DispatchObject::NONE.
  */
void Command::ListCommands(DispatchObject::Otype typeIn) {
  if (typeIn == DispatchObject::NONE) {
    for (int idx = 1; idx != DispatchObject::HIDDEN; idx++)
      ListCommandsForType( (DispatchObject::Otype)idx );
  } else
    ListCommandsForType( typeIn );
}

/** \return true if any control blocks remain. */
bool Command::UnterminatedControl() {
  if (!control_.empty()) {
    mprinterr("Error: %u unterminated control block(s) detected.\n", ctlidx_+1);
    for (int i = 0; i <= ctlidx_; i++)
      mprinterr("Error:   %i : %s\n", i, control_[i]->Description().c_str());
    return true;
  }
  return false;
}

/** Create new control block with given Block. */
int Command::AddControlBlock(ControlBlock* ctl, CpptrajState& State, ArgList& cmdArg) {
  if ( ctl->SetupBlock( State, cmdArg ) )
    return 1;
  if (ctlidx_ == -1) mprintf("CONTROL: Starting control block.\n");
  control_.push_back( ctl );
  ctlidx_++;
  mprintf("  BLOCK %2i: ", ctlidx_);
  for (int i = 0; i < ctlidx_; i++)
    mprintf("  ");
  mprintf("%s\n", ctl->Description().c_str());
  //mprintf("DEBUG: Begin control block %i\n", ctlidx_);
  return 0;
}

#define NEW_BLOCK "__NEW_BLOCK__"
/** Execute the specified control block. */
int Command::ExecuteControlBlock(int block, CpptrajState& State)
{
  control_[block]->Start();
  ControlBlock::DoneType ret = control_[block]->CheckDone(CurrentVars_);
  if (State.Debug() > 0) {
    mprintf("DEBUG: Start: CurrentVars:");
    CurrentVars_.PrintVariables();
  }
  while (ret == ControlBlock::NOT_DONE) {
    for (ControlBlock::const_iterator it = control_[block]->begin();
                                      it != control_[block]->end(); ++it)
    {
      if (it->CommandIs(NEW_BLOCK)) {
        // Execute next control block
        if (ExecuteControlBlock(block+1, State)) return 1;
      } else {
        for (int i = 0; i < block; i++) mprintf("  ");
        // Execute command
        if ( ExecuteCommand(State, *it) != CpptrajState::OK ) return 1;
      }
    }
    ret = control_[block]->CheckDone(CurrentVars_);
  }
  if (ret == ControlBlock::ERROR) return 1;
  return 0;
}

/** Handle the given command. If inside a control block, if the command is
  * a control command a new block will be created, otherwise the command will
  * be added to the current control block. Once all control blocks are
  * complete they will be executed. If not inside a control block, just
  * execute the command.
  */
CpptrajState::RetType Command::Dispatch(CpptrajState& State, std::string const& commandIn)
{
  ArgList cmdArg( commandIn );
  cmdArg.MarkArg(0); // Always mark the first arg as the command
  // Check for control block
  if (!control_.empty()) {
    mprintf("  [%s]\n", cmdArg.ArgLine());
    // In control block.
    if ( control_[ctlidx_]->EndBlock( cmdArg ) ) {
      // End the current control block.
      //mprintf("DEBUG: End control block %i.\n", ctlidx_);
      mprintf("  BLOCK %2i: ", ctlidx_);
      for (int i = 0; i < ctlidx_; i++)
        mprintf("  ");
      mprintf("END\n");
      ctlidx_--;
      if (ctlidx_ < 0) {
        // Outermost control structure is ended. Execute control block(s).
        mprintf("CONTROL: Executing %u control block(s).\n", control_.size());
        if (State.QuietBlocks()) SetWorldSilent(true);
        int cbret = ExecuteControlBlock(0, State);
        ClearControlBlocks();
        if (State.QuietBlocks()) SetWorldSilent(false);
        if (cbret != 0) return CpptrajState::ERR;
        mprintf("CONTROL: Control block finished.\n\n");
      }
    } else {
      // Check if this is another control block statement (silently)
      Cmd const& ctlCmd = SearchTokenType(DispatchObject::CONTROL, cmdArg.Command(), true);
      if (ctlCmd.Empty() || ctlCmd.Destination() != Cmd::BLK) { // TODO just check Destination?
        // Add this command to current control block.
        control_[ctlidx_]->AddCommand( cmdArg );
        mprintf("\tAdded command '%s' to control block %i.\n", cmdArg.Command(), ctlidx_);
      } else {
        // Tell current block that a new block is being created
        control_[ctlidx_]->AddCommand(NEW_BLOCK);
        // Create new control block
        DispatchObject* obj = ctlCmd.Alloc();
        if (AddControlBlock( (ControlBlock*)obj, State, cmdArg )) {
          delete obj;
          ClearControlBlocks();
          return CpptrajState::ERR;
        }
      }
    }
    return CpptrajState::OK;
  }
  return ExecuteCommand( State, cmdArg );
}

#undef NEW_BLOCK

/** Search for the given command and execute it. Replace any variables in
  * command with their values. EXE and CTL commands are executed immediately
  * and then freed. ACT and ANA commands are sent to the CpptrajState for later
  * execution. BLK commands set up control blocks which will be executed when
  * the outer control block is completed.
  */
CpptrajState::RetType Command::ExecuteCommand( CpptrajState& State, ArgList const& cmdArgIn ) {
  // Replace variable names in command with entries from CurrentVars
  ArgList cmdArg = CurrentVars_.ReplaceVariables( cmdArgIn, State.DSL(), State.Debug() );
  if (cmdArg.empty()) return CpptrajState::ERR;
  // Print modified command
  mprintf("  [%s]\n", cmdArg.ArgLine());
  // Look for command in command list.
  Cmd const& cmd = SearchToken( cmdArg );
  CpptrajState::RetType ret_val = CpptrajState::OK;
  if (cmd.Empty()) {
    // Not a command. Try to evaluate the expression
    RPNcalc calc;
    calc.SetDebug( State.Debug() );
    if (calc.ProcessExpression( cmdArg.ArgLineStr() ))
      ret_val = CpptrajState::ERR;
    else {
      if (calc.Evaluate( State.DSL() ))
        ret_val = CpptrajState::ERR;
    }
    if (ret_val == CpptrajState::ERR)
      mprinterr("'%s': Invalid command or expression.\n", cmdArg.ArgLine());
  } else {
    DispatchObject* obj = cmd.Alloc();
    switch (cmd.Destination()) {
      case Cmd::CTL:
        ret_val = ((Control*)obj)->SetupControl(State, cmdArg, CurrentVars_);
        delete obj;
        break;
      case Cmd::BLK:
        if (AddControlBlock( (ControlBlock*)obj, State, cmdArg )) {
          delete obj;
          return CpptrajState::ERR;
        }
        break;
      case Cmd::EXE:
#       ifdef MPI
        if (Parallel::TrajComm().Master())
          ret_val = ((Exec*)obj)->Execute( State, cmdArg );
        else
          ret_val = CpptrajState::OK;
        if (Parallel::World().CheckError( (int)ret_val ) != 0)
          ret_val = CpptrajState::ERR;
#       else
        ret_val = ((Exec*)obj)->Execute( State, cmdArg );
#       endif
        delete obj;
        break;
      case Cmd::ACT: ret_val = State.AddToActionQueue( (Action*)obj, cmdArg ); break;
      case Cmd::ANA: ret_val = State.AddToAnalysisQueue( (Analysis*)obj, cmdArg ); break;
      case Cmd::DEP:
        mprinterr("Error: '%s' is deprecated.\n", cmdArg.Command());
        cmd.Help();
        break;
    }
  }
  return ret_val;
}

/** Read command input from file. */
CpptrajState::RetType Command::ProcessInput(CpptrajState& State, std::string const& inputFilename)
{
  BufferedLine infile;
  if (infile.OpenFileRead( inputFilename )) {
    if (!inputFilename.empty())
      mprinterr("Error: Could not open input file '%s'\n", inputFilename.c_str());
    return CpptrajState::ERR;
  }
  mprintf("INPUT: Reading input from '%s'\n", infile.Filename().full());
  // Read in each line of input.
  int nInputErrors = 0;
  CpptrajState::RetType cmode = CpptrajState::OK;
  CmdInput input;
  const char* ptr = infile.Line();
  while (ptr != 0) {
    bool moreInput = input.AddInput( ptr );
    while (moreInput) {
      ptr = infile.Line();
      moreInput = input.AddInput( ptr );
    }
    // Only attempt to execute if the command is not blank.
    if (!input.Empty()) {
#     ifdef TIMER
      Timer time_cmd; // DEBUG
      time_cmd.Start(); // DEBUG
#     endif
      // Call Dispatch to convert input to ArgList and process.
      cmode = Command::Dispatch(State, input.Str());
#     ifdef TIMER
      time_cmd.Stop(); // DEBUG
      time_cmd.WriteTiming(0," Command time: "); // DEBUG
#     endif
      if (cmode == CpptrajState::ERR) {
        nInputErrors++;
        if (State.ExitOnError()) break;
      } else if (cmode == CpptrajState::QUIT)
        break;
    }
    // Reset Input line
    input.Clear();
    ptr = infile.Line();
  }
  infile.CloseFile();
  if (nInputErrors > 0) {
    mprinterr("\t%i errors encountered reading input.\n", nInputErrors);
    return CpptrajState::ERR;
  }
  return cmode;
}
