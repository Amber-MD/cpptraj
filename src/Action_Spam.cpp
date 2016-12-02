// Action_Spam
#include <cmath> // sqrt
#include <cstdio> // sscanf
#include "Action_Spam.h"
#include "Box.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "Constants.h" // ELECTOAMBER
#include "DistRoutines.h"
#include "StringRoutines.h" // integerToString
#include "KDE.h"
#include "OnlineVarT.h" // Stats
#include "DataSet_Mesh.h"

// CONSTRUCTOR
Action_Spam::Action_Spam() : Action(HIDDEN),
  debug_(0),
  DG_BULK_(-30.3), // Free energy of bulk SPCE water
  DH_BULK_(-22.2), // Enthalpy of bulk SPCE water
  temperature_(300.0),
  purewater_(false),
  reorder_(false),
  calcEnergy_(false),
  cut2_(144.0),
  onecut2_(1.0 / 144.0),
  doublecut_(24.0),
  infofile_(0),
  site_size_(2.5),
  sphere_(false),
  ds_dg_(0),
  ds_dh_(0),
  ds_ds_(0),
  Nframes_(0),
  overflow_(false)
{ }

void Action_Spam::Help() const {
  mprintf("\t<filename> [solv <solvname>] [reorder] [name <name>]\n"
          "\t[purewater] [cut <cut>] [info <infofile>] [summary <summary>]\n"
          "\t[site_size <size>] [sphere] [out <datafile>]\n"
          "\t[dgbulk <dgbulk>] [dhbulk <dhbulk>] [temperature <T>]\n"
          "    <filename> : File with the peak locations present (XYZ- format)\n"
          "    <solvname> : Name of the solvent residues\n"
          "    <cut>      : Non-bonded cutoff for energy evaluation\n"
          "    <dgbulk>   : SPAM free energy of the bulk solvent in kcal/mol\n"
          "    <dhbulk>   : SPAM enthalpy of the bulk solvent in kcal/mol\n"
          "    <T>        : Temperature at which SPAM calculation was run.\n"
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
  trajComm_ = init.TrajComm();
# endif
  debug_ = debugIn;
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

  // Get energy cutoff
  double cut = actionArgs.getKeyDouble("cut", 12.0);
  cut2_ = cut * cut;
  doublecut_ = 2 * cut;
  onecut2_ = 1 / cut2_;

  if (purewater_) {
    // We only have one data set averaging over every water. Add it here
    DataSet* ds = init.DSL().AddSet(DataSet::DOUBLE, MetaData(ds_name));
    if (ds == 0) return Action::ERR;
    if (datafile != 0) datafile->AddDataSet( ds );
    myDSL_.push_back( ds );
    DG_BULK_ = 0.0;
    DH_BULK_ = 0.0;
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
    DG_BULK_ = actionArgs.getKeyDouble("dgbulk", -30.3);
    if (!actionArgs.Contains("dgbulk"))
      mprintf("Warning: 'dgbulk' not specified; using default for SPC/E water.\n");
    DH_BULK_ = actionArgs.getKeyDouble("dhbulk", -22.2);
    if (!actionArgs.Contains("dhbulk"))
      mprintf("Warning: 'dhbulk' not specified; using default for SPC/E water.\n");
    temperature_ = actionArgs.getKeyDouble("temperature", 300.0);
    std::string infoname = actionArgs.GetStringKey("info");
    if (infoname.empty())
      infoname = std::string("spam.info");
    infofile_ = init.DFL().AddCpptrajFile(infoname, "SPAM info");
    if (infofile_ == 0) return Action::ERR;
    DataFile* summaryfile = init.DFL().AddDataFile(actionArgs.GetStringKey("summary"), actionArgs);
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
    // Now add all of the individual peak energy data sets
    for (int i = 0; i < (int)peaks_.size(); i++) {
      DataSet* ds = init.DSL().AddSet(DataSet::DOUBLE, MetaData(ds_name,i+1));
      if (ds == 0) return Action::ERR;
      myDSL_.push_back( ds );
      if (datafile != 0) datafile->AddDataSet( ds );
    }
    // Make the peak overall energy sets Mesh so we can skip unoccupied peaks
    Dimension Pdim( 1.0, 0.0, "Peak" );
    ds_dg_ = init.DSL().AddSet(DataSet::XYMESH, MetaData(ds_name,"DG"));
    ds_dh_ = init.DSL().AddSet(DataSet::XYMESH, MetaData(ds_name,"DH"));
    ds_ds_ = init.DSL().AddSet(DataSet::XYMESH, MetaData(ds_name,"-TDS"));
    if (ds_dg_==0 || ds_dh_==0 || ds_ds_==0) return Action::ERR;
    ds_dg_->SetDim(Dimension::X, Pdim);
    ds_dh_->SetDim(Dimension::X, Pdim);
    ds_ds_->SetDim(Dimension::X, Pdim);
    if (summaryfile != 0) {
      summaryfile->AddDataSet( ds_dg_ );
      summaryfile->AddDataSet( ds_dh_ );
      summaryfile->AddDataSet( ds_ds_ );
    }
#   ifdef MPI
    ds_dg_->SetNeedsSync(false);
    ds_dh_->SetNeedsSync(false);
    ds_ds_->SetNeedsSync(false);
#   endif
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
  mprintf("    SPAM:\n");
  if (purewater_) {
    mprintf("\tCalculating bulk value for pure solvent\n");
    if (datafile != 0)
      mprintf("\tPrinting solvent energies to %s\n", datafile->DataFilename().full());
    mprintf("\tUsing a %.2f Angstrom non-bonded cutoff with shifted EEL.\n",
            sqrt(cut2_));
    if (reorder_)
      mprintf("\tWarning: Re-ordering makes no sense for pure solvent.\n");
    if (!summaryfile_.empty())
      mprintf("\tPrinting solvent SPAM summary to %s\n", summaryfile_.c_str());
  } else {
    mprintf("\tSolvent [%s] density peaks taken from %s.\n",
            solvname_.c_str(), filename.base());
    mprintf("\t%zu density peaks will be analyzed from %s.\n",
            peaks_.size(), filename.base());
    mprintf("\tOccupation information printed to %s.\n", infofile_->Filename().full());
    mprintf("\tSites are ");
    if (sphere_)
      mprintf("spheres with diameter %.3f\n", site_size_);
    else
      mprintf("boxes with edge length %.3f\n", site_size_);
    if (reorder_)
      mprintf("\tRe-ordering trajectory so each site always has the same water molecule.\n");
    if (!calcEnergy_) {
      if (!reorder_) {
        mprinterr("Error: Not re-ordering trajectory or calculating energies. Nothing to do!\n");
        return Action::ERR;
      }
      mprintf("\tNot calculating any SPAM energies\n");
    } else {
      mprintf("\tUsing a non-bonded cutoff of %.2f Ang. with a EEL shifting function.\n",
              sqrt(cut2_));
      mprintf("\tBulk solvent SPAM free energy: %.3f kcal/mol\n", DG_BULK_);
      mprintf("\tBulk solvent SPAM enthalpy: %.3f kcal/mol\n", DH_BULK_);
      mprintf("\tTemperature: %.3f K\n", temperature_);
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
  if (solvent_residues_.empty()) {
    mprinterr("Error: No solvent residues found with name '%s'\n", solvname_.c_str());
    return Action::ERR;
  }
  resPeakNum_.reserve( solvent_residues_.size() );
  comlist_.reserve( solvent_residues_.size() );

  mprintf("\tFound %zu solvent residues [%s]\n", solvent_residues_.size(),
          solvname_.c_str());

  // Set up the charge array and check that we have enough info
  if (SetupParms(setup.Top())) return Action::ERR;

  // Save topology address so we can get NB params during energy calc. 
  CurrentParm_ = setup.TopAddress();

  return Action::OK;
}

// Action_Spam::SetupParms
/** Sets the temporary charge array and makes sure that we have the necessary
  * parameters in our topology to calculate nonbonded energy terms
  */
int Action_Spam::SetupParms(Topology const& ParmIn) {
  // Store the charges, convert to Amber style so we get energies in kcal/mol
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

  // Now loop through all atoms in the residue and loop through the pairlist to
  // get the energies
  for (int i = res.FirstAtom(); i < res.LastAtom(); i++) {
    Vec3 atm1 = Vec3(frameIn.XYZ(i));
    for (int j = 0; j < CurrentParm_->Natom(); j++) {
      if (j >= res.FirstAtom() && j < res.LastAtom()) continue;
      Vec3 atm2 = Vec3(frameIn.XYZ(j));
      double dist2;
      // Get imaged distance
      switch( image_.ImageType() ) {
        case NONORTHO : dist2 = DIST2_ImageNonOrtho(atm1, atm2, ucell_, recip_); break;
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
  // Calculate unit cell and fractional matrices for non-orthorhombic system
  if ( image_.ImageType() == NONORTHO )
    frm.Frm().BoxCrd().ToRecip(ucell_, recip_);
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
  t_action_.Start();
  int wat = 0;
  int maxwat = (int)solvent_residues_.size();
  int basenum = frameNum * solvent_residues_.size();
  DataSet_double& evals = static_cast<DataSet_double&>( *myDSL_[0] );
  // Make room for each solvent residue energy this frame.
  evals.Resize( evals.Size() + solvent_residues_.size() );
  t_energy_.Start();
# ifdef _OPENMP
# pragma omp parallel private(wat)
  {
# pragma omp for
# endif
  for (wat = 0; wat < maxwat; wat++)
    evals[basenum + wat] = Calculate_Energy(frameIn, solvent_residues_[wat]);
# ifdef _OPENMP
  }
# endif
  t_energy_.Stop();
  t_action_.Stop();
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
  t_action_.Start();
  t_resCom_.Start();
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
  t_resCom_.Stop();
  t_assign_.Start();
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
  t_assign_.Stop();
  t_occupy_.Start();
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
  t_occupy_.Stop();
  t_energy_.Start();
  // If we have to calculate energies, do that here
  if (calcEnergy_) {
    int peak;
    int npeaks = (int)peaks_.size();
    const double ZERO = 0.0;
    // Loop through every peak, then loop through the water molecules to find
    // which one is in that site, and calculate the LJ and EEL energies for that
    // water molecule within a given cutoff.
#   ifdef _OPENMP
#   pragma omp parallel private(peak)
    {
#   pragma omp for schedule(dynamic)
#   endif
    for (peak = 0; peak < npeaks; peak++)
    {
      if (occupied[peak]) {
        for (unsigned int i = 0; i < resPeakNum_.size(); i++)
          if (resPeakNum_[i] == peak) {
            // Now we have our residue number. Create a pairlist for each solvent
            // molecule that can be used for each atom in that residue. Should
            // provide some time savings.
            double ene = Calculate_Energy(frameIn, solvent_residues_[i]);
            myDSL_[peak]->Add(frameNum, &ene);
            break;
          }
      } else
        myDSL_[peak]->Add(frameNum, &ZERO);
    }
#   ifdef _OPENMP
    }
#   endif
  }
  t_energy_.Stop();
  t_reordr_.Start();
  // If we have to re-order trajectories, do that here
  Action::RetType ret = Action::OK;
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
    ret = MODIFY_COORDS;
  }
  t_reordr_.Stop();
  t_action_.Stop();

  return ret;
}

static inline int absval(int i) { if (i < 0) return -i; else return i; }

/** Calculate the DELTA G of an individual water site */
int Action_Spam::Calc_G_Wat(DataSet* dsIn, unsigned int peaknum)
{
  DataSet_1D const& dataIn = static_cast<DataSet_1D const&>( *dsIn );
  // Create energy vector containing only frames that are singly-occupied.
  // Calculate the mean (enthalpy) while doing this.
  DataSet_double enevec;
  Stats<double> Havg;
  double min = 0.0, max = 0.0;
  if (!peakFrameData_.empty()) {
    Iarray const& SkipFrames = peakFrameData_[peaknum];
    Iarray::const_iterator fnum = SkipFrames.begin();
    for (int frm = 0; frm != (int)dataIn.Size(); frm++) {
      bool frameIsSkipped = (fnum != SkipFrames.end() && absval(*fnum) == frm);
      if (frameIsSkipped)
        ++fnum;
      else {
        double ene = dataIn.Dval(frm);
        if (enevec.Size() < 1) {
          min = ene; 
          max = ene;
        } else {
          min = std::min(min, ene);
          max = std::max(max, ene);
        }
        enevec.AddElement( ene );
        Havg.accumulate( ene );
      }
    }
  } else {
    min = dataIn.Dval(0);
    max = dataIn.Dval(0);
    for (unsigned int frm = 0; frm != dataIn.Size(); frm++) {
      double ene = dataIn.Dval(frm);
      min = std::min(min, ene);
      max = std::max(max, ene);
      enevec.AddElement( ene );
      Havg.accumulate( ene );
    }
  }
  if (enevec.Size() < 1)
    return 1;
  // Calculate distribution of energy values using KDE. Get the bandwidth
  // factor here since we already know the SD.
  double BWfac = KDE::BandwidthFactor( enevec.Size() );
  //mprintf("DEBUG:\tNvals=%zu min=%g max=%g BWfac=%g\n", enevec.Size(), min, max, BWfac);
  // Estimate number of bins the same way spamstats.py does.
  int nbins = (int)(((max - min) / BWfac) + 0.5) + 100;

  HistBin Xdim(nbins, min - (50*BWfac), BWfac, "P(Ewat)");
  //Xdim.CalcBinsOrStep(min - Havg.variance(), max + Havg.variance(), 0.0, nbins, "P(Ewat)");
  if (debug_ > 0) Xdim.PrintHistBin();
  DataSet_double kde1;
  KDE gkde;
  if (gkde.CalcKDE( kde1, enevec, Xdim, 1.06 * sqrt(Havg.variance()) * BWfac )) {
    mprinterr("Error: Could not calculate E KDE histogram.\n");
    return -1;
  }
  kde1.SetupFormat() = TextFormat(TextFormat::GDOUBLE, 12, 5);
  // Determine SUM[ P(Ewat) * exp(-Ewat / RT) ]
  double RT = Constants::GASK_KCAL * temperature_;
  double KB = 1.0 / RT;
  double sumQ = 0.0;
  for (unsigned int i = 0; i != kde1.Size(); i++) {
    double Ewat = kde1.Xcrd(i);
    double PEwat = kde1.Dval(i);
    sumQ += (PEwat * exp( -Ewat * KB ));
  }
  //mprintf("DEBUG: sumQ= %20.10E\n", sumQ);
  double DG = -RT * log(BWfac * sumQ);

  double adjustedDG = DG - DG_BULK_;
  double adjustedDH = Havg.mean() - DH_BULK_;
  double ntds = adjustedDG - adjustedDH;
  if (ds_dg_ == 0) {
    mprintf("\tSPAM bulk energy values:\n");
    mprintf("\t  <G>= %g, <H>= %g +/- %g, -TdS= %g\n", adjustedDG, adjustedDH,
            sqrt(Havg.variance()), ntds);
  } else {
    ((DataSet_Mesh*)ds_dg_)->AddXY(peaknum+1, adjustedDG);
    ((DataSet_Mesh*)ds_dh_)->AddXY(peaknum+1, adjustedDH);
    ((DataSet_Mesh*)ds_ds_)->AddXY(peaknum+1, ntds);
  }
  
  // DEBUG
  if (debug_ > 1) {
    FileName rawname("dbgraw." + integerToString(peaknum+1) + ".dat");
    FileName kdename("dbgkde." + integerToString(peaknum+1) + ".dat");
    mprintf("DEBUG: Writing peak %u raw energy values to '%s', KDE histogram to '%s'\n",
            peaknum+1, rawname.full(), kdename.full());
    DataFile rawout;
    rawout.SetupDatafile( rawname, 0 );
    rawout.AddDataSet( &enevec );
    rawout.WriteDataOut();
    DataFile kdeout;
    kdeout.SetupDatafile( kdename, 0 );
    kdeout.AddDataSet( &kde1 );
    kdeout.WriteDataOut();
  }

  return 0;
}

#ifdef MPI
int Action_Spam::SyncAction() {
  // Get total number of frames.
  std::vector<int> frames_on_rank( trajComm_.Size() );
  int myframes = Nframes_;
  trajComm_.GatherMaster( &myframes, 1, MPI_INT, &frames_on_rank[0] );
  if (trajComm_.Master())
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      Nframes_ += frames_on_rank[rank];

  // Sync peakFrameData_
  std::vector<int> size_on_rank( trajComm_.Size() );
  for (unsigned int i = 0; i != peakFrameData_.size(); i++)
  {
    Iarray& Data = peakFrameData_[i];
    int mysize = (int)Data.size();
    trajComm_.GatherMaster( &mysize, 1, MPI_INT, &size_on_rank[0] );
    if (trajComm_.Master()) {
      int total = size_on_rank[0];
      for (int rank = 1; rank < trajComm_.Size(); rank++)
        total += size_on_rank[rank];
      Data.resize( total );
      int* endptr = &(Data[0]) + size_on_rank[0];
      // Receive data from each rank
      for (int rank = 1; rank < trajComm_.Size(); rank++) {
        trajComm_.SendMaster( endptr, size_on_rank[rank], rank, MPI_INT );
        // Properly offset the frame numbers
        for (int j = 0; j != size_on_rank[rank]; j++, endptr++)
          if (*endptr < 0)
            *endptr -= frames_on_rank[rank-1];
          else
            *endptr += frames_on_rank[rank-1];
      }
    } else // Send data to master
      trajComm_.SendMaster( &(Data[0]), Data.size(), trajComm_.Rank(), MPI_INT );
  }
  return 0;
}
#endif

// Action_Spam::Print()
void Action_Spam::Print() {
  // Timings
  mprintf("\tSPAM timing data:\n");
  t_resCom_.WriteTiming(2, "Residue c.o.m. calc:", t_action_.Total());
  t_assign_.WriteTiming(2, "Peak assignment    :", t_action_.Total());
  t_occupy_.WriteTiming(2, "Occupancy calc.    :", t_action_.Total());
  t_energy_.WriteTiming(2, "Energy calc        :", t_action_.Total());
  t_reordr_.WriteTiming(2, "Residue reordering :", t_action_.Total());
  t_action_.WriteTiming(1, "SPAM Action Total:");
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
                        i+1, (int)peakFrameData_[i].size(), ndouble);
      for (unsigned int j = 0; j < peakFrameData_[i].size(); j++) {
        if (j > 0 && j % 10 == 0) infofile_->Printf("\n");
        // Adjust frame number.
        int fnum = peakFrameData_[i][j];
        if (fnum < 0)
         fnum--;
        else
          fnum++;
        infofile_->Printf(" %7d", fnum);
      }
      infofile_->Printf("\n\n");
    }

    unsigned int p = 0;
    int n_peaks_no_energy = 0;
    for (std::vector<DataSet*>::const_iterator ds = myDSL_.begin(); ds != myDSL_.end(); ++ds, ++p)
    {
      int err = Calc_G_Wat( *ds, p );
      if (err == 1)
        n_peaks_no_energy++;
      else if (err == -1)
        mprintf("Warning: Error calculating SPAM energies for peak %zu\n", ds - myDSL_.begin());
    }
    if (n_peaks_no_energy > 0)
      mprintf("Warning: No energies for %i peaks.\n", n_peaks_no_energy);
  } else
    Calc_G_Wat( myDSL_[0], 0 );
}
