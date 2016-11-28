// Action_Spam
#include <cmath> // sqrt
#include <cstdio> // sscanf, sprintf
#include "Action_Spam.h"
#include "Box.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "Constants.h" // ELECTOAMBER
#include "DataSet_integer.h"
#include "DistRoutines.h"
#include "StringRoutines.h" // integerToString
#include "KDE.h"

// CONSTRUCTOR
Action_Spam::Action_Spam() : Action(HIDDEN),
  bulk_(0.0),
  purewater_(false),
  reorder_(false),
  calcEnergy_(false),
  cut2_(144.0),
  onecut2_(1.0 / 144.0),
  doublecut_(24.0),
  infofile_(0),
  site_size_(2.5),
  sphere_(false),
  Nframes_(0),
  overflow_(false)
  ,set_counter_(0) // DEBUG
{ }

void Action_Spam::Help() const {
  mprintf("\t<filename> [solv <solvname>] [reorder] [name <name>] [bulk <value>]\n"
          "\t[purewater] [cut <cut>] [info <infofile>] [summary <summary>]\n"
          "\t[site_size <size>] [sphere] [out <datafile>]\n\n"
          "    <filename> : File with the peak locations present (XYZ- format)\n"
          "    <solvname> : Name of the solvent residues\n"
          "    <cut>      : Non-bonded cutoff for energy evaluation\n"
          "    <value>    : SPAM free energy of the bulk solvent\n"
          "    <infofile> : File with stats about which sites are occupied when.\n"
          "    <size>     : Size of the water site around each density peak.\n"
          "    [sphere]   : Treat each site like a sphere.\n"
          "    [purewater]: The system is pure water---used to parametrize the bulk values.\n"
          "    [reorder]  : The solvent should be re-ordered so the same solvent molecule\n"
          "                 is always in the same site.\n"
          "    <summary>  : File with the summary of all SPAM results. If not specified,\n"
          "                 no SPAM energies will be calculated.\n"
          "    <datafile> : Data file with all SPAM energies for each snapshot.\n");
}

// Action_Spam::Init()
Action::RetType Action_Spam::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'spam' action does not work with > 1 thread (%i threads currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  // Always use imaged distances
  image_.InitImaging(true);
  // This is needed everywhere in this function scope
  FileName filename;

  // See if we're doing pure water. If so, we don't need a peak file
  purewater_ = actionArgs.hasKey("purewater");

  // Get data set name.
  std::string ds_name = actionArgs.GetStringKey("name");
  if (ds_name.empty())
    ds_name = init.DSL().GenerateDefaultName("SPAM");

  // Get output data file
  DataFile* datafile = init.DFL().AddDataFile(actionArgs.GetStringKey("out"), actionArgs);

  // Solvent residue name
  solvname_ = actionArgs.GetStringKey("solv");
  if (solvname_.empty())
    solvname_ = std::string("WAT");

  if (purewater_) {
    // We still need the cutoff
    double cut = actionArgs.getKeyDouble("cut", 12.0);
    cut2_ = cut * cut;
    doublecut_ = 2 * cut;
    onecut2_ = 1 / cut2_;
    // We only have one data set averaging over every water. Add it here
    DataSet* ds = init.DSL().AddSet(DataSet::DOUBLE, MetaData(ds_name));
    if (ds == 0) return Action::ERR;
    if (datafile != 0) datafile->AddDataSet( ds );
    myDSL_.push_back( ds );
  } else {
    // Get the file name with the peaks defined in it
    filename.SetFileName( actionArgs.GetStringNext() );
    if (filename.empty()) {
      mprinterr("Error: No Peak file specified.\n");
      return Action::ERR;
    } else if (!File::Exists(filename)) {
      File::ErrorMsg( filename.full() );
      return Action::ERR;
    }
    // Get the remaining optional arguments
    reorder_ = actionArgs.hasKey("reorder");
    bulk_ = actionArgs.getKeyDouble("bulk", 0.0);
    double cut = actionArgs.getKeyDouble("cut", 12.0);
    cut2_ = cut * cut;
    doublecut_ = 2 * cut;
    onecut2_ = 1 / cut2_;
    std::string infoname = actionArgs.GetStringKey("info");
    if (infoname.empty())
      infoname = std::string("spam.info");
    infofile_ = init.DFL().AddCpptrajFile(infoname, "SPAM info");
    if (infofile_ == 0) return Action::ERR;
    // The default maskstr is the Oxygen atom of the solvent
    summaryfile_ = actionArgs.GetStringKey("summary");
    // Divide site size by 2 to make it half the edge length (or radius)
    site_size_ = actionArgs.getKeyDouble("site_size", 2.5) / 2.0;
    sphere_ = actionArgs.hasKey("sphere");
    // If it's a sphere, square the radius to compare with
    if (sphere_)
      site_size_ *= site_size_;

    // Parse through the peaks file and extract the peaks
    CpptrajFile peakfile;
    if (peakfile.OpenRead(filename)) {
      mprinterr("SPAM: Error: Could not open %s for reading!\n", filename.full());
      return Action::ERR;
    }
    std::string line = peakfile.GetLine();
    int npeaks = 0;
    while (!line.empty()) {
      if (sscanf(line.c_str(), "%d", &npeaks) != 1) {
        line = peakfile.GetLine();
        continue;
      }
      line = peakfile.GetLine();
      break;
    }
    while (!line.empty()) {
      double x, y, z, dens;
      if (sscanf(line.c_str(), "C %lg %lg %lg %lg", &x, &y, &z, &dens) != 4) {
        line = peakfile.GetLine();
        continue;
      }
      line = peakfile.GetLine();
      peaks_.push_back(Vec3(x, y, z));
    }
    peakfile.CloseFile();
    // Check that our initial number of peaks matches our parsed peaks. Warn
    // otherwise
    if (npeaks != (int)peaks_.size())
      mprinterr("SPAM: Warning: %s claims to have %d peaks, but really has %d!\n",
                filename.full(), npeaks, peaks_.size());
    // Now add all of the data sets
    MetaData md(ds_name);
    for (int i = 0; i < (int)peaks_.size(); i++) {
      md.SetAspect( integerToString(i+1) ); // TODO: Should this be Idx?
      DataSet* ds = init.DSL().AddSet(DataSet::DOUBLE, md);
      if (ds == 0) return Action::ERR;
      myDSL_.push_back( ds );
      if (datafile != 0) datafile->AddDataSet( ds );
    }
    // peakFrameData will keep track of omitted frames for each peak.
    peakFrameData_.clear();
    peakFrameData_.resize( peaks_.size() );
  }
  // Determine if energy calculation needs to happen
  calcEnergy_ = (!summaryfile_.empty() || datafile != 0);
  // Set function for determining if water is inside peak
  if (sphere_)
    Inside_ = &Action_Spam::inside_sphere;
  else
    Inside_ = &Action_Spam::inside_box;

  // Print info now
  if (purewater_) {
    mprintf("SPAM: Calculating bulk value for pure solvent\n");
    if (datafile != 0)
      mprintf("SPAM: Printing solvent energies to %s\n", datafile->DataFilename().full());
    mprintf("SPAM: Using a %.2f Angstrom non-bonded cutoff with shifted EEL.\n",
            sqrt(cut2_));
    if (reorder_)
      mprintf("SPAM: Warning: Re-ordering makes no sense for pure solvent.\n");
    if (!summaryfile_.empty())
      mprintf("SPAM: Printing solvent SPAM summary to %s\n", summaryfile_.c_str());
  }else {
    mprintf("SPAM: Solvent [%s] density peaks taken from %s.\n",
            solvname_.c_str(), filename.base());
    mprintf("SPAM: %d density peaks will be analyzed from %s.\n",
            peaks_.size(), filename.base());
    mprintf("SPAM: Occupation information printed to %s.\n", infofile_->Filename().full());
    mprintf("SPAM: Sites are ");
    if (sphere_)
      mprintf("spheres with diameter %.3lf\n", site_size_);
    else
      mprintf("boxes with edge length %.3lf\n", site_size_);
    if (reorder_) {
      mprintf("SPAM: Re-ordering trajectory so each site always has ");
      mprintf("the same water molecule.\n");
    }
    if (!calcEnergy_) {
      if (!reorder_) {
        mprinterr("SPAM: Error: Not re-ordering trajectory or calculating energies. ");
        mprinterr("Nothing to do!\n");
        return Action::ERR;
      }
      mprintf("SPAM: Not calculating any SPAM energies\n");
    }else {
      mprintf("SPAM: Using a non-bonded cutoff of %.2lf Ang. with a EEL shifting function.\n",
              sqrt(cut2_));
      mprintf("SPAM: Bulk solvent SPAM energy taken as %.3lf kcal/mol\n", bulk_);
    }
  }
  mprintf("#Citation: Cui, G.; Swails, J.M.; Manas, E.S.; \"SPAM: A Simple Approach\n"
          "#          for Profiling Bound Water Molecules\"\n"
          "#          J. Chem. Theory Comput., 2013, 9 (12), pp 5539–5549.\n");

  return Action::OK;
}

// Action_Spam::Setup()
Action::RetType Action_Spam::Setup(ActionSetup& setup) {
  // We need box info
  Box const& currentBox = setup.CoordInfo().TrajBox();
  if (currentBox.Type() == Box::NOBOX) {
    mprinterr("Error: SPAM: Must have explicit solvent with periodic boundaries!\n");
    return Action::ERR;
  }

  // See if our box dimensions are too small for our cutoff...
  if (currentBox.BoxX() < doublecut_ ||
      currentBox.BoxY() < doublecut_ ||
      currentBox.BoxZ() < doublecut_)
  {
    mprinterr("Error: SPAM: The box appears to be too small for your cutoff!\n");
    return Action::ERR;
  }
  // Set up the solvent_residues_ vector
  int resnum = 0;
  for (Topology::res_iterator res = setup.Top().ResStart();
                              res != setup.Top().ResEnd(); res++)
  {
    if (res->Name().Truncated() == solvname_) {
      solvent_residues_.push_back(*res);
      // Tabulate COM
      double mass = 0.0;
      for (int i = res->FirstAtom(); i < res->LastAtom(); i++)
        mass += setup.Top()[i].Mass();
    }
    resnum++;
  }
  resPeakNum_.reserve( solvent_residues_.size() );
  comlist_.reserve( solvent_residues_.size() );

  mprintf("\tFound %zu solvent residues [%s]\n", solvent_residues_.size(),
          solvname_.c_str());

  // Set up the charge array and check that we have enough info
  if (SetupParms(setup.Top())) return Action::ERR;

  // Back up the parm
  // NOTE: This is a full copy - use reference instead?
  CurrentParm_ = setup.TopAddress();

  return Action::OK;
}

// Action_Spam::SetupParms
/** Sets the temporary charge array and makes sure that we have the necessary
  * parameters in our topology to calculate nonbonded energy terms
  */
int Action_Spam::SetupParms(Topology const& ParmIn) {
  // Store the charges
  atom_charge_.clear();
  atom_charge_.reserve( ParmIn.Natom() );
  for (Topology::atom_iterator atom = ParmIn.begin(); atom != ParmIn.end(); ++atom)
    atom_charge_.push_back( atom->Charge() * Constants::ELECTOAMBER );
  if (!ParmIn.Nonbond().HasNonbond()) {
    mprinterr("Error: SPAM: Parm does not have LJ information.\n");
    return 1;
  }
  return 0;
}

// Action_Spam::Calculate_Energy()
double Action_Spam::Calculate_Energy(Frame const& frameIn, Residue const& res) {
  // The first atom of the solvent residue we want the energy from
  double result = 0;
  Matrix_3x3 ucell, recip;
  if ( image_.ImageType() == NONORTHO )
    frameIn.BoxCrd().ToRecip(ucell, recip);
  /* Now loop through all atoms in the residue and loop through the pairlist to
   * get the energies
   */
  for (int i = res.FirstAtom(); i < res.LastAtom(); i++) {
    Vec3 atm1 = Vec3(frameIn.XYZ(i));
    for (int j = 0; j < CurrentParm_->Natom(); j++) {
      if (j >= res.FirstAtom() && j < res.LastAtom()) continue;
      Vec3 atm2 = Vec3(frameIn.XYZ(j));
      double dist2;
      // Get imaged distance
      switch( image_.ImageType() ) {
        case NONORTHO : dist2 = DIST2_ImageNonOrtho(atm1, atm2, ucell, recip); break;
        case ORTHO    : dist2 = DIST2_ImageOrtho(atm1, atm2, frameIn.BoxCrd()); break;
        default       : dist2 = DIST2_NoImage(atm1, atm2); break;
      }
      if (dist2 < cut2_) {
        double qiqj = atom_charge_[i] * atom_charge_[j];
        NonbondType const& LJ = CurrentParm_->GetLJparam(i, j);
        double r2 = 1 / dist2;
        double r6 = r2 * r2 * r2;
                  // Shifted electrostatics: qiqj/r * (1-r/rcut)^2 + VDW
        double shift = (1 - dist2 * onecut2_);
        result += qiqj / sqrt(dist2) * shift * shift + LJ.A() * r6 * r6 - LJ.B() * r6;
      }
    }
  }
  return result;
}

// Action_Spam::DoAction()
Action::RetType Action_Spam::DoAction(int frameNum, ActionFrame& frm) {
  Nframes_++;

  // Check that our box is still big enough...
  overflow_ = overflow_ || frm.Frm().BoxCrd().BoxX() < doublecut_ ||
                           frm.Frm().BoxCrd().BoxY() < doublecut_ ||
                           frm.Frm().BoxCrd().BoxZ() < doublecut_;
  if (purewater_)
    return DoPureWater(frameNum, frm.Frm());
  else
    return DoSPAM(frameNum, frm.ModifyFrm());
}

// Action_Spam::DoPureWater
/** Carries out SPAM analysis for pure water to parametrize bulk */
  /* This is relatively simple... We have to loop through every water molecule
   * for every frame, calculate the energy of that solvent molecule, and add
   * that to our one data set. Therefore we will have NFRAMES * NWATER data
   * points
   */
Action::RetType Action_Spam::DoPureWater(int frameNum, Frame const& frameIn)
{
  int wat = 0;
  int basenum = frameNum * solvent_residues_.size();
  for (std::vector<Residue>::const_iterator res = solvent_residues_.begin();
        res != solvent_residues_.end(); res++) {
    double ene = Calculate_Energy(frameIn, *res);
    myDSL_[0]->Add(basenum + wat, &ene);
    wat++;
  }
  return Action::OK;
}

// Action_Spam::inside_box()
bool Action_Spam::inside_box(Vec3 gp, Vec3 pt, double edge) const {
  return (gp[0] + edge > pt[0] && gp[0] - edge < pt[0] &&
          gp[1] + edge > pt[1] && gp[1] - edge < pt[1] &&
          gp[2] + edge > pt[2] && gp[2] - edge < pt[2]);
}

// Action_Spam::inside_sphere()
bool Action_Spam::inside_sphere(Vec3 gp, Vec3 pt, double rad2) const {
  return ( (gp[0]-pt[0])*(gp[0]-pt[0]) + (gp[1]-pt[1])*(gp[1]-pt[1]) +
           (gp[2]-pt[2])*(gp[2]-pt[2]) < rad2 );
}

// Action_Spam::DoSPAM
/** Carries out SPAM analysis on a typical system */
Action::RetType Action_Spam::DoSPAM(int frameNum, Frame& frameIn) {
  /* A list of all solvent residues and the sites that they are reserved for. An
   * unreserved solvent residue has an index -1. At the end, we will go through
   * and re-order the frame if requested.
   */
  resPeakNum_.assign(solvent_residues_.size(), -1);
  // Tabulate all of the COMs
  comlist_.clear();
  for (std::vector<Residue>::const_iterator res = solvent_residues_.begin();
                                            res != solvent_residues_.end(); res++)
    comlist_.push_back(frameIn.VCenterOfMass(res->FirstAtom(), res->LastAtom()));
  // Loop through each peak and then scan through every residue, and assign a
  // solvent residue to each peak
  int pknum = 0;
  for (Varray::const_iterator pk = peaks_.begin(); pk != peaks_.end(); ++pk, ++pknum)
  {
    for (unsigned int resnum = 0; resnum != comlist_.size(); resnum++)
    {
      // If we're inside, make sure this residue is not already `claimed'. If it
      // is, assign it to the closer peak center
      if ((this->*Inside_)(*pk, comlist_[resnum], site_size_)) {
        if (resPeakNum_[resnum] > 0) {
          Vec3 diff1 = comlist_[resnum] - *pk;
          Vec3 diff2 = comlist_[resnum] - peaks_[ resPeakNum_[resnum] ];
          // If we are closer, update. Otherwise do nothing
          if (diff1.Magnitude2() < diff2.Magnitude2())
            resPeakNum_[resnum] = pknum;
        } else
          resPeakNum_[resnum] = pknum;
      }
    }
  }

  /* Now we have a vector of reservations. We want to make sure that each site
   * is occupied once and only once. If a site is unoccupied, add frameNum to
   * this peak's data set in peakFrameData_. If a site is double-occupied, add
   * -frameNum to this peak's data set in peakFrameData_.
   */
  std::vector<bool> occupied(peaks_.size(), false);
  std::vector<bool> doubled(peaks_.size(), false); // to avoid double-additions
  for (Iarray::const_iterator it = resPeakNum_.begin();
                              it != resPeakNum_.end(); it++)
  {
    if (*it > -1) {
      if (!occupied[*it])
        occupied[*it] = true;
      else if (!doubled[*it]) {
        peakFrameData_[*it].push_back(-frameNum); // double-occupied, frameNum will be ignored
        doubled[*it] = true;
      }
    }
  }
  // Now loop through and add all non-occupied sites
  for (unsigned int i = 0; i < peaks_.size(); i++)
    if (!occupied[i]) 
      peakFrameData_[i].push_back(frameNum);
  // Now adjust the occupied vectors to only contain 'true' for sites we need to
  // analyze (i.e., make all doubled points 'unoccupied')
  for (unsigned int i = 0; i < peaks_.size(); i++)
    if (doubled[i])
      occupied[i] = false;

  // If we have to calculate energies, do that here
  if (calcEnergy_) {
    /* Loop through every peak, then loop through the water molecules to find
     * which one is in that site, and calculate the LJ and EEL energies for that
     * water molecule within a given cutoff.
     */
    for (int peak = 0; peak < (int)peaks_.size(); peak++) {
      // Skip unoccupied peaks
      if (!occupied[peak]) {
        double val = 0;
        myDSL_[peak]->Add(frameNum, &val);
        continue;
      }
      for (unsigned int i = 0; i < resPeakNum_.size(); i++) {
        if (resPeakNum_[i] != peak) continue;
        /* Now we have our residue number. Create a pairlist for each solvent
         * molecule that can be used for each atom in that residue. Should
         * provide some time savings.
         */
        double ene = Calculate_Energy(frameIn, solvent_residues_[i]);
        myDSL_[peak]->Add(frameNum, &ene);
        break;
      }
    }
  }

  // If we have to re-order trajectories, do that here
  if (reorder_) {
    /* Loop over every occupied site and swap the atoms so the same solvent
     * residue is always in the same site
     */
    for (int i = 0; i < (int)peaks_.size(); i++) {
      // Skip unoccupied sites
      if (!occupied[i]) continue;
      for (unsigned int j = 0; j < solvent_residues_.size(); j++) {
        // This is the solvent residue in our site
        if (resPeakNum_[j] == i) {
          for (int k = 0; k < solvent_residues_[j].NumAtoms(); k++)
            frameIn.SwapAtoms(solvent_residues_[i].FirstAtom()+k,
                              solvent_residues_[j].FirstAtom()+k);
          // Since we swapped solvent_residues_ of 2 solvent atoms, we also have
          // to swap reservations[i] and reservations[j]...
          int tmp = resPeakNum_[j];
          resPeakNum_[j] = resPeakNum_[i];
          resPeakNum_[i] = tmp;
        }
      }
    }
    return MODIFY_COORDS;
  }

  return Action::OK;
}

/** Calculate the DELTA G of an individual water site */
void Action_Spam::Calc_G_Wat(DataSet* dsIn, Iarray const& SkipFrames,
                             double& dg_avg, double& dg_std, double& dh_avg,
                             double& dh_std, double& ntds)
{
  DataSet_1D const& dataIn = static_cast<DataSet_1D const&>( *dsIn );
  // Create energy vector containing only frames that are singly-occupied
  DataSet_double enevec;
  Iarray::const_iterator fnum = SkipFrames.begin();
  int frm = 0;
  for (; frm != (int)dataIn.Size(); frm++) {
    if (fnum == SkipFrames.end()) break;
    int skippedFrame = *fnum;
    if (skippedFrame < 0) skippedFrame = -skippedFrame;
    if (frm == skippedFrame)
      ++fnum;
    else
      enevec.AddElement( dataIn.Dval(frm) );
  }
  for (; frm != (int)dataIn.Size(); frm++)
    enevec.AddElement( dataIn.Dval(frm) );
  if (enevec.Size() < 1) {
    mprintf("Warning: No energy values for peak. Skipping.\n");
    return;
  }

  DataSet_double kde1;
  //global DG_BULK, DH_BULK, BW
  KDE gkde;
  if (gkde.CalcKDE( kde1, enevec )) {
    mprinterr("Error: Could not calculate E KDE histogram.\n");
    return;
  }
  // DEBUG
  DataFile rawout;
  rawout.SetupDatafile( FileName("dbgraw." + integerToString(set_counter_)   + ".dat"), 0 );
  rawout.AddDataSet( dsIn );
  rawout.WriteDataOut();
  DataFile kdeout;
  kdeout.SetupDatafile( FileName("dbgkde." + integerToString(set_counter_++) + ".dat"), 0 );
  kdeout.AddDataSet( &kde1 );
  kdeout.WriteDataOut();

  // def calc_g_wat(enevec, sample_size, sample_num)
/*
  // Set up defaults
  unsigned int sample_size = 0;
  unsigned int sample_num = 0;
  if (sample_size < 1) sample_size = enevec.Size();
  sample_size = std::min(enevec.Size(), sample_size);
  unsigned int sample_num = std::max(1, sample_num);

  std::vector<double> enthalpy(sample_num, 0.0);
  std::vector<double> free_energy(sample_num, 0.0);
  // We want to do "sample_num" subsamples with "sample_size" elements in each
  // subsample.  Loop over those now
  for (unsigned int ii = 0; ii != sample_num; ii++) {
    // If sample_num is 1, then don't resample
      if (sample_num == 1)
         skde = gaussian_kde(enevec);
      else
         // Generate the subsample kernel
         skde = gaussian_kde(kde1.resample(sample_size))
      // Determine the widths of the bins from the covariance factor method
      binwidth = skde.covariance_factor()
      // Determine the number of bins over the range of the data set
      nbins = int(((skde.dataset.max() - skde.dataset.min()) / binwidth) + 0.5) \
            + 100
      // Make a series of points from the minimum to the maximum
      pts = np.zeros(nbins)
      for i in range(nbins): pts[i] = skde.dataset.min() + (i-50) * binwidth
      // Evaluate the KDE at all of those points
      kdevals = skde.evaluate(pts)
      // Get the enthalpy and free energy
      enthalpy[ii] = np.sum(skde.dataset) / sample_size
      free_energy[ii] = -1.373 * log10(binwidth *
                                      np.sum(kdevals * np.exp(-pts / 0.596)))

   // Our subsampling is over, now find the average and standard deviation
   dg_avg = np.sum(free_energy) / sample_num - DG_BULK
   dg_std = free_energy.std()
   dh_avg = np.sum(enthalpy) / sample_num - DH_BULK
   dh_std = enthalpy.std()
   ntds = dg_avg - dh_avg
*/
}

// Action_Spam::Print()
void Action_Spam::Print() {
  // Print the spam info file if we didn't do pure water
  if (!purewater_) {
    // Warn about any overflows
    if (overflow_)
      mprinterr("Warning: SPAM: Some frames had a box too small for the cutoff.\n");

    // Print information about each missing peak
    infofile_->Printf("# There are %d density peaks and %d frames\n\n",
                (int)peaks_.size(), Nframes_);
    // Loop over every Data set
    for (unsigned int i = 0; i < peakFrameData_.size(); i++) {
      // Skip peaks with 0 unoccupied sites
      if (peakFrameData_[i].size() == 0) continue;
      // Find out how many double-occupied frames there are
      int ndouble = 0;
      for (unsigned int j = 0; j < peakFrameData_[i].size(); j++)
        if (peakFrameData_[i][j] < 0) ndouble++;
      infofile_->Printf("# Peak %u has %d omitted frames (%d double-occupied)\n",
                        i, (int)peakFrameData_[i].size(), ndouble);
      for (unsigned int j = 0; j < peakFrameData_[i].size(); j++) {
        if (j > 0 && j % 10 == 0) infofile_->Printf("\n");
        infofile_->Printf("%7d ", peakFrameData_[i][j]);
      }
      infofile_->Printf("\n\n");
    }

    double dg_avg, dg_std, dh_avg, dh_std, ntds;
    unsigned int p = 0;
    for (std::vector<DataSet*>::const_iterator ds = myDSL_.begin(); ds != myDSL_.end(); ++ds, ++p)
      Calc_G_Wat( *ds, peakFrameData_[p], dg_avg, dg_std, dh_avg, dh_std, ntds );
  }

  // Print the summary file with the calculated SPAM energies
  if (!summaryfile_.empty()) {
    // Not enabled yet -- just print out the data files.
    mprinterr("Warning: SPAM: SPAM calculation not yet enabled.\n");
    //if (datafile_.empty()) datafile_ = summaryfile_;
  }
}
