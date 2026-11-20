#include "Analysis_MDANCE.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"
#include "DataSet_integer_mem.h" // for JSON output and MDANCE results
#include "ProgressBar.h"
#include "StringRoutines.h" // integerToString
#include "Trajout_Single.h"
#include <algorithm> // std::sort
#ifdef HAS_EIGEN
# include "Mdance/KMeans.h"
# include "Mdance/helm.h"
# include "Mdance/result_utils.h"
using namespace Cpptraj::Mdance;
#endif

/** CONSTRUCTOR */
Analysis_MDANCE::Analysis_MDANCE() :
  debug_(0),
  coords_(0),
  kClusters_(0),
  kSeed_(0),
  percentage_(0),
  vthresh_(0),
  metric_(ExtendedSimilarity::NO_METRIC),
#ifdef HAS_EIGEN
  kinit_(MD::KinitType::StratAll),
#endif
  cnumvtime_(0),
  centers_(0),
  clusterfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  centerfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  infofile_(0),
  summaryfile_(0),
  jsonfile_(0),
  sort_(true)
{}

/** DESTRUCTOR */
Analysis_MDANCE::~Analysis_MDANCE() {
}

const char* Analysis_MDANCE::kinitKeys_[] = {
  "all",
  "reduced",
  "compsim",
  "divselect",
  "kmeanspp",
  "random",
  "vanillakmpp",
  0
};

const char* Analysis_MDANCE::kinitStr_[] = {
  "All",
  "Reduced",
  "CompSim",
  "DivSelect",
  "Kmeans++",
  "Random",
  "Vanilla Kmeans++"
};

#ifdef HAS_EIGEN
const Cpptraj::Mdance::MD::KinitType Analysis_MDANCE::kinitTypes_[] = {
  MD::KinitType::StratAll,
  MD::KinitType::StratReduced,
  MD::KinitType::CompSim,
  MD::KinitType::DivSelect,
  MD::KinitType::KmeansPP,
  MD::KinitType::Random,
  MD::KinitType::VanillaKmeansPP
};
#endif

// Analysis_MDANCE::Help()
void Analysis_MDANCE::Help() const {
# ifdef HAS_EIGEN
  mprintf("\t[crdset <COORDS set>] clusters <#> [mask <mask>]\n"
          "\t[metric <metric>] [vthresh <vectthreshold>]\n"
          "\t[kinit <init>] [pct <percentage>] [kseed <#>]\n"
          "\t[name <set name>] [out <cnumvtime file>]\n"
          "\t[clusterout <trajfileprefix> [clusterfmt <trajformat>]]\n"
          "\t[centerout <trajfilename> [centerfmt <trajformat>]]\n"
          "\t[info <infofile>] [summary <summaryfile>] [nosort]\n"
          "\t[json <jsonfile>]\n"
         );
  mprintf("  <metric> = %s\n", ExtendedSimilarity::MetricKeys().c_str());
  mprintf("  <init>   =");
  for (int i = 0; kinitKeys_[i] != 0; i++)
    mprintf(" %s", kinitKeys_[i]);
  mprintf("\n");
# else
  mprintf("CPPTRAJ was compiled without Eigen - MDANCE is disabled.\n");
# endif
}


/** Get arguments related to writing cluster data to trajectories.
  * Copied from Cluster/Results_Coords.
  */
void Analysis_MDANCE::getClusterTrajArgs(ArgList& argIn,
                                         const char* trajKey, const char* fmtKey,
                                         std::string& trajName,
                                         TrajectoryFile::TrajFormatType& fmt) const
{
  trajName = argIn.GetStringKey( trajKey );
  fmt = TrajectoryFile::WriteFormatFromString( argIn.GetStringKey(fmtKey), fmt );
  // If file name specified but not format, try to guess from name
  if (!trajName.empty() && fmt == TrajectoryFile::UNKNOWN_TRAJ)
    fmt = TrajectoryFile::WriteFormatFromFname( trajName, TrajectoryFile::AMBERTRAJ );
}


// Analysis_MDANCE::Setup()
Analysis::RetType Analysis_MDANCE::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
# ifdef HAS_EIGEN
  debug_ = debugIn;
  sort_ = !analyzeArgs.hasKey("nosort");
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: Could not locate COORDS set corresponding to '%s'\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
  // Target # of clusters, percentage, threshold
  kClusters_ = analyzeArgs.getKeyInt("clusters", 0);
  kSeed_ = analyzeArgs.getKeyInt("kseed", 0);
  percentage_ = analyzeArgs.getKeyInt("pct", 10);
  vthresh_ = analyzeArgs.getKeyInt("vthresh", 16);
  // Metric
  std::string mstr = analyzeArgs.GetStringKey("metric");
  metric_ = ExtendedSimilarity::NO_METRIC;
  if (!mstr.empty()) {
    metric_ = ExtendedSimilarity::TypeFromKeyword( mstr );
    if (metric_ == ExtendedSimilarity::NO_METRIC) {
      mprinterr("Error: Metric '%s' not recognized.\n", mstr.c_str());
      return Analysis::ERR;
    }
  } else {
    metric_ = ExtendedSimilarity::MSD;
  }
  // Init strategy
  std::string kstr = analyzeArgs.GetStringKey("kinit");
  int iKinit = -1;
  if (!kstr.empty()) {
    for (int i = 0; kinitKeys_[i] != 0; i++) {
      const char* key = kinitKeys_[i];
      if (key != 0 && kstr == std::string(key)) {
        iKinit = i;
        break;
      }
    }
    if (iKinit < 0) {
      mprinterr("Error: Unrecognized keyword for 'kinit': %s\n", kstr.c_str());
      return Analysis::ERR;
    }
    kinit_ = kinitTypes_[iKinit];
  } else {
    kinit_ = MD::KinitType::StratAll;
    iKinit = 0;
  }
  // Atom mask
  std::string maskstr = analyzeArgs.GetStringKey("mask");
  if (maskstr.empty())
    maskstr = "*";
  if (mask_.SetMaskString(maskstr)) {
    mprinterr("Error: Could not set mask string '%s'\n", maskstr.c_str());
    return Analysis::ERR;
  }
  // Set up results that depend on COORDS DataSet
  getClusterTrajArgs(analyzeArgs, "clusterout",   "clusterfmt",   clusterfile_,  clusterfmt_);
  getClusterTrajArgs(analyzeArgs, "centerout",    "centerfmt",    centerfile_,   centerfmt_);
  // Output files/data
  DataFile* cnumvtimefile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  infofile_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("info"), "MDANCE cluster info",
                                         DataFileList::TEXT, true);
  summaryfile_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("summary"), "MDANCE cluster summary",
                                        DataFileList::TEXT, true);
  if (infofile_ == 0 || summaryfile_ == 0) {
    mprinterr("Error: Could not allocate cluster info/summary file.\n");
    return Analysis::ERR;
  }
  jsonfile_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("json"), "MDANCE JSON file",
                                         DataFileList::TEXT, false);
  // Overall set name extracted here. All other arguments should already be processed. 
  std::string dsname = analyzeArgs.GetStringKey("name");
  if (dsname.empty())
    dsname = setup.DSL().GenerateDefaultName("MDANCE");
  // ---------------------------------------------
    
  // Cluster number vs time data set
  cnumvtime_ = setup.DSL().AddSet(DataSet::INTEGER, dsname, "Cnum");
  if (cnumvtime_ == 0) return Analysis::ERR;
  if (cnumvtimefile != 0) cnumvtimefile->AddDataSet( cnumvtime_ );
  // Cluster centers data set
  centers_ = (DataSet_Coords*)setup.DSL().AddSet(DataSet::COORDS, MetaData(dsname, "centers"));

  // Check input
  if (kClusters_ < 1) {
    mprinterr("Error: 'clusters' must be > 0.\n");
    return Analysis::ERR;
  }

  // Info
  mprintf("    MDANCE:\n");
  mprintf("\t----==== Clustering Options ====----\n");
  mprintf("\t  COORDS set       : %s\n", coords_->legend());
  mprintf("\t  # clusters       : %i\n", kClusters_);
  mprintf("\t  KMeans seed      : %i\n", kSeed_);
  mprintf("\t  Metric           : %s\n", ExtendedSimilarity::metricStr(metric_));
  mprintf("\t  Init. Strat.     : %s\n", kinitStr_[iKinit]);
  mprintf("\t  Percentage       : %i%%\n", percentage_);
  mprintf("\t  Vect. threshhold : %i\n", vthresh_);
  mprintf("\t  Atom selection   : %s\n", mask_.MaskString());

  mprintf("\n");
  mprintf("\t----==== Output Options ====----\n");
  //mprintf("\tData set name          : %s\n", dsname.c_str());
  mprintf("\tCluster # vs time set  : %s\n", cnumvtime_->Meta().PrintName().c_str());
  mprintf("\tCluster centers set    : %s\n", centers_->Meta().PrintName().c_str());
  if (sort_)
    mprintf("\tWill sort clusters by population.\n");
  else
    mprintf("\tNot sorting clusters by population.\n");
  mprintf("\tInfo file              : %s\n", infofile_->Filename().full());
  mprintf("\tSummary output file    : %s\n", summaryfile_->Filename().full());
  if (cnumvtimefile != 0)
    mprintf("\tCluster # vs time file : %s\n", cnumvtimefile->DataFilename().full());
  if (jsonfile_ != 0)
    mprintf("\tMDANCE JSON file       : %s\n", jsonfile_->Filename().full());
  if (!clusterfile_.empty())
    mprintf("\tCluster trajectories will be written to %s.cX, format %s\n",
            clusterfile_.c_str(), TrajectoryFile::FormatString(clusterfmt_));
  if (!centerfile_.empty())
    mprintf("\tCluster centers will be written to %s, format %s\n",
            centerfile_.c_str(), TrajectoryFile::FormatString(centerfmt_));

  mprintf("# Citation: Lexin Chen, Daniel R. Roe, Matthew Kochert, Carlos Simmerling,\n"
          "#           Ramón Alain Miranda-Quintana;\n"
          "#           k‑Means NANI: An Improved Clustering Algorithm for Molecular Dynamics Simulations.\n"
          "#           J. Chem. Theory Comput. 9 July 2024; 20 (13): 5583–5597.\n"
          "#            https://doi.org/10.1021/acs.jctc.4c00308\n");

  return Analysis::OK;
# else /* HAS_EIGEN */
  mprintf("CPPTRAJ was compiled without Eigen - MDANCE is disabled.\n");
  return Analysis::ERR;
# endif /* HAS_EIGEN */
}

/** Write frames in each cluster to a trajectory file.  */
void Analysis_MDANCE::writeClusterTraj(ClusterArray const& Clusters) const {
  Topology* clusterparm = coords_->TopPtr();
  // Loop over all clusters
  for (unsigned int cidx = 0; cidx != Clusters.size(); cidx++)
  {
    Iarray const& cluster = Clusters[cidx].Frames();
    // Create filename based on cluster number.
    std::string cfilename =  clusterfile_ + ".c" + integerToString( cidx );
    // Set up trajectory file 
    Trajout_Single clusterout;
    if (clusterout.PrepareTrajWrite(cfilename, ArgList(), DataSetList(), clusterparm,
                                    coords_->CoordsInfo(), cluster.size(),
                                    clusterfmt_))
    {
      mprinterr("Error: Could not set up cluster trajectory %s for write.\n",
                cfilename.c_str());  
      return;
    } 
    // Loop over all frames in cluster
    unsigned int set = 0;
    Frame clusterframe = coords_->AllocateFrame();
    for (Iarray::const_iterator fnum = cluster.begin(); fnum != cluster.end(); ++fnum)
    {
      coords_->GetFrame( *fnum, clusterframe );
      clusterout.WriteSingle(set++, clusterframe);
    }
    // Close traj
    clusterout.EndTraj();
  }
}

/** Write cluster centers to a trajectory file.  */
void Analysis_MDANCE::writeCenterTraj(ClusterArray const& Clusters) const {
  // Set up trajectory file 
  Trajout_Single clusterout;
  if (clusterout.PrepareTrajWrite(centerfile_, ArgList(), DataSetList(), centers_->TopPtr(),
                                  centers_->CoordsInfo(), centers_->Size(),
                                  centerfmt_))
  {
    mprinterr("Error: Could not set up cluster centers trajectory %s for write.\n",
              centerfile_.c_str());  
    return;
  }
  // Loop over all clusters
  for (unsigned int cidx = 0; cidx != Clusters.size(); cidx++)
  {
    Frame const& centerframe = Clusters[cidx].Ctr();
    clusterout.WriteSingle(cidx, centerframe);
  }
  // Close traj
  clusterout.EndTraj();
}

/** Write summary to given file */
void Analysis_MDANCE::writeSummary(CpptrajFile& outfile, ClusterArray const& Clusters) const {
  outfile.Printf("%-8s %8s %8s %8s %8s\n","#Cluster","Frames","Frac","MSD","Rep");
  for (ClusterArray::const_iterator clust = Clusters.begin(); clust != Clusters.end(); ++clust)
  {
    double frac = (double)clust->size() / (double)Clusters.Nframes();
    outfile.Printf("%8li %8u %8.3f %8.3f %8i\n", clust-Clusters.begin(), clust->size(), frac, clust->MSD(), clust->Rep()+1);
  }
}

/** Write info to given file */
void Analysis_MDANCE::writeInfo(CpptrajFile& outfile, ClusterArray const& Clusters)
const
{
  std::string buffer;
  
  outfile.Printf("#Clustering: %u clusters %u frames\n",
                 Clusters.size(), Clusters.Nframes());
  // DBI
  outfile.Printf("#DBI: %f\n", Clusters.DBI());
  // Pseudo-F
  if (Clusters.size() > 1) {
    //double SSRSST = 0.0;
    //double pseudof = clusters.ComputePseudoF( SSRSST, metricIn );
    outfile.Printf("#pSF: %f\n", Clusters.PSF());
    //outfile.Printf("#SSR/SST: %f\n", SSRSST);
  } else
    mprintf("Warning: Fewer than 2 clusters. Not calculating pseudo-F.\n");

  // Do not print trajectory stuff if no filename given (i.e. STDOUT output)
  if (!outfile.IsStream()) {
    for (ClusterArray::const_iterator C1 = Clusters.begin(); C1 != Clusters.end(); ++C1)
    {
      buffer.clear();
      buffer.resize(Clusters.Nframes(), '.');
      for (Iarray::const_iterator f1 = C1->Frames().begin(); f1 != C1->Frames().end(); ++f1)
        buffer[ *f1 ] = 'X';
      buffer += '\n';
      outfile.Write((void*)buffer.c_str(), buffer.size());
    }
  }
  // Print representative frame numbers
  outfile.Printf("#Representative frames:");
  for (ClusterArray::const_iterator C1 = Clusters.begin(); C1 != Clusters.end(); ++C1)
    outfile.Printf(" %i", C1->Rep()+1);
  outfile.Printf("\n");
}

/** Write cluster results in MDANCE Json format */
void Analysis_MDANCE::writeJson(CpptrajFile& outfile, ClusterArray const& Clusters)
const
{
  outfile.Printf("{\n");
  outfile.Printf("  \"algorithm\": \"kmeans\",\n"); // FIXME
  outfile.Printf("  \"nFrames\": %u,\n", Clusters.Nframes());
  outfile.Printf("  \"nClusters\": %u,\n", Clusters.size());
  outfile.Printf("  \"labels\": [");
  std::string buffer;
  DataSet_integer const& CVT = static_cast<DataSet_integer const&>( *cnumvtime_ );
  for (unsigned int idx = 0; idx != CVT.Size(); idx++) {
    if (idx > 0)
      buffer.append(", " + integerToString(CVT[idx]));
    else
      buffer.append(integerToString(CVT[idx]));
  }
  outfile.Write((void*)buffer.c_str(), buffer.size());
  outfile.Printf("],\n");
  outfile.Printf("  \"clusterSizes\": [");
  for (unsigned int cnum = 0; cnum != Clusters.size(); cnum++) {
    if (cnum > 0)
      outfile.Printf(", %u", Clusters[cnum].size());
    else
      outfile.Printf("%u", Clusters[cnum].size());
  }
  outfile.Printf("],\n");
  outfile.Printf("  \"representatives\": [");
  for (unsigned int cnum = 0; cnum != Clusters.size(); cnum++) {
    if (cnum > 0)
      outfile.Printf(", %i", Clusters[cnum].Rep());
    else
      outfile.Printf("%i", Clusters[cnum].Rep());
  }
  outfile.Printf("],\n");
  outfile.Printf("  \"clusterMSD\": [");
  for (unsigned int cnum = 0; cnum != Clusters.size(); cnum++) {
    if (cnum > 0)
      outfile.Printf(", %.10g", Clusters[cnum].MSD());
    else
      outfile.Printf("%.10g", Clusters[cnum].MSD());
  }
  outfile.Printf("],\n");
  outfile.Printf("  \"scores\": {\n");
  outfile.Printf("    \"calinskiHarabasz\": %.10g,\n", Clusters.PSF());
  outfile.Printf("    \"daviesBouldin\": %.10g\n", Clusters.DBI());
  outfile.Printf("  }\n");
  outfile.Printf("}\n");
}

// Analysis_MDANCE::Analyze()
Analysis::RetType Analysis_MDANCE::Analyze() {
# ifdef HAS_EIGEN
  mprintf("    MDANCE: Starting MDANCE.\n");
  if (coords_ == 0) {
    mprinterr("Error: COORDS are null.\n");
    return Analysis::ERR;
  }
  DataSet_Coords& CRD = static_cast<DataSet_Coords&>( *coords_ );
  // Get the atom selection
  if (CRD.Top().SetupIntegerMask( mask_ )) {
    mprinterr("Error: Could not set up atom mask.\n");
    return Analysis::ERR;
  }
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: No atoms selected.\n");
    return Analysis::ERR;
  }
  // Set the topology/COORDS set for centers
  CoordinateInfo ctrInfo; // Coordinates only
  if (mask_.Nselected() < CRD.Top().Natom()) {
    // Strip top to match clustered coords
    Topology* clusterTop = CRD.Top().modifyStateByMask( mask_ );
    if (clusterTop == 0) {
      mprinterr("Error: Could not create topology for cluster centers.\n");
      return Analysis::ERR;
    }
    clusterTop->Brief("Topology for cluster centers");
    if (centers_->CoordsSetup( *clusterTop, ctrInfo )) {
      mprinterr("Error: Could not set up COORDS set for cluster centers.\n");
      return Analysis::ERR;
    }
    delete clusterTop;
  } else {
    if (centers_->CoordsSetup( CRD.Top(), ctrInfo )) {
      mprinterr("Error: Could not set up COORDS set for cluster centers.\n");
      return Analysis::ERR;
    }
  }
  centers_->Allocate(DataSet::SizeArray(1, kClusters_));
  // This is an Eigen matrix. Each row is a frame, each column is a coordinate.
  unsigned int nSelectedAtoms = mask_.Nselected();
  unsigned int ncoords = nSelectedAtoms * 3;
  ArrayXXd data( CRD.Size(), ncoords );
  mprintf("\tSaving internal Eigen matrix (%zd rows/frames, %zd cols/coords).\n", data.rows(), data.cols());
  ProgressBar progress(CRD.Size());
  Frame frmIn = CRD.AllocateFrame();
  for (unsigned int idx = 0; idx != CRD.Size(); idx++)
  {
    progress.Update(idx);
    CRD.GetFrame(idx, frmIn);
    unsigned int icrd = 0;
    for (unsigned int iat = 0; iat < nSelectedAtoms; iat++)
    {
      const double* XYZ = frmIn.XYZ(mask_[iat]);
      data( idx, icrd   ) = XYZ[0];
      data( idx, icrd+1 ) = XYZ[1];
      data( idx, icrd+2 ) = XYZ[2];
      icrd += 3;
    }
  }

  // Convert metric to internal MDANCE
  MD::Metric mt = MD::Metric::MSD;
  switch(metric_) {
    case ExtendedSimilarity::MSD : mt = MD::Metric::MSD; break;
    case ExtendedSimilarity::BUB : mt = MD::Metric::BUB; break;
    case ExtendedSimilarity::FAI : mt = MD::Metric::Fai; break;
    case ExtendedSimilarity::GLE : mt = MD::Metric::Gle; break;
    case ExtendedSimilarity::JA  : mt = MD::Metric::Ja; break;
    case ExtendedSimilarity::JT  : mt = MD::Metric::JT; break;
    case ExtendedSimilarity::RT  : mt = MD::Metric::RT; break;
    case ExtendedSimilarity::RR  : mt = MD::Metric::RR; break;
    case ExtendedSimilarity::SM  : mt = MD::Metric::SM; break;
    case ExtendedSimilarity::SS1 : mt = MD::Metric::SS1; break;
    case ExtendedSimilarity::SS2 : mt = MD::Metric::SS2; break;
    case ExtendedSimilarity::NO_METRIC :
      mprinterr("Internal Error: Analysis_MDANCE::Analyze(): No metric.\n");
      return Analysis::ERR;
  }
  // Initialize and run Kmeans
  KmeansNANI kmeans(data, kClusters_, mt, kinit_, nSelectedAtoms, percentage_, vthresh_, kSeed_);

  // Results
  // First check the clustering assignments.
  // MDANCE labels each frame with the cluster number
  ClusterArray Clusters( coords_->Size() );
  Clusters.resize( kClusters_ );
  Veci cluster_of_frame = kmeans.getLabels();
  cnumvtime_->Allocate(DataSet::SizeArray(1, cluster_of_frame.size()));
  for (int ifrm = 0; ifrm < cluster_of_frame.size(); ifrm++) {
    int cnum = cluster_of_frame[ifrm];
    if (cnum < 0 || cnum >= kClusters_) {
      mprinterr("Error: Cluster of frame %i is out of bounds: %i\n", ifrm+1, cnum);
    }
    Clusters[cnum].push_back( ifrm );
    cnumvtime_->Add(ifrm, &cnum);
    if (debug_ > 1) mprintf("DEBUG: Frame %8i Cluster %8i\n", ifrm+1, cnum);
  }
  // Get the pseudo-F (Calinski-Harabasz) and DBI scores
  std::pair<double,double> scores = kmeans.computeScores();
  Clusters.SetDBI( scores.second );
  Clusters.SetPSF( scores.first );
  //mprintf("\tDBI      : %f\n", scores.second);
  //mprintf("\tpseudo-F : %f\n", scores.first);
  std::vector<int> reps = Cpptraj::Mdance::computeRepresentatives( data,
                                                                   (((DataSet_integer_mem*)cnumvtime_)->Array()),
                                                                   Clusters.size(),
                                                                   nSelectedAtoms,
                                                                   mt );
  for (std::vector<int>::const_iterator it = reps.begin(); it != reps.end(); ++it)
    Clusters[it-reps.begin()].SetRep( *it );
  std::vector<double> MSDS = Cpptraj::Mdance::computeClusterMSD( data,
                                                                 (((DataSet_integer_mem*)cnumvtime_)->Array()),
                                                                  Clusters.size(),
                                                                  nSelectedAtoms );
  for (std::vector<double>::const_iterator it = MSDS.begin(); it != MSDS.end(); ++it)
    Clusters[it-MSDS.begin()].SetMSD( *it );

  // Get centers
  Frame ctrFrame = centers_->AllocateFrame();
  Mat clusterCenters = kmeans.getCenters();
  if (debug_ > 0) mprintf("DEBUG: centers rows %zd, cols %zd\n", clusterCenters.rows(), clusterCenters.cols());
  for (int iclust = 0; iclust != kClusters_; iclust++) {
    unsigned int icrd = 0;
    ctrFrame.ClearAtoms();
    for (unsigned int iat = 0; iat < nSelectedAtoms; iat++) {
      double XYZ[3];
      XYZ[0] = clusterCenters( iclust, icrd   );
      XYZ[1] = clusterCenters( iclust, icrd+1 );
      XYZ[2] = clusterCenters( iclust, icrd+2 );
      ctrFrame.AddXYZ( XYZ );
      icrd += 3;
    }
    //centers_->AddFrame( ctrFrame );
    Clusters[iclust].SetCtr( ctrFrame );
  }

  // Write JSON before sorting since MDANCE does not sort
  if (jsonfile_ != 0)
    writeJson(*jsonfile_, Clusters);

  // Sort if needed
  if (sort_)
    std::sort(Clusters.begin(), Clusters.end());

  // Add centers to the centers DataSet
  for (int iclust = 0; iclust != kClusters_; iclust++)
    centers_->AddFrame( Clusters[iclust].Ctr() );

  // Write info
  writeInfo(*infofile_, Clusters);
  // Write summary
  writeSummary(*summaryfile_, Clusters);
  // Write cluster trajectories
  if (!clusterfile_.empty())
    writeClusterTraj( Clusters );
  if (!centerfile_.empty())
    writeCenterTraj( Clusters );

  return Analysis::OK; // DEBUG
# else /* HAS_EIGEN */
  return Analysis::ERR;
# endif /* HAS_EIGEN */
}
