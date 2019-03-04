#include "Action_XtalSymm.h"
#include "CpptrajStdio.h"
#include "Matrix.h"
#include "Molecule.h"
#include "Topology.h"
#include "SpaceGroup.h"
#include <cmath> // floor

const int Action_XtalSymm::IASU_GRID_BINS_ = 192;

const double Action_XtalSymm::DASU_GRID_BINS_ = 192;

Action_XtalSymm::Action_XtalSymm() :
  nops_(0),
  nCopyA_(0),
  nCopyB_(0),
  nCopyC_(0),
  sgID_(0),
  refType_(NONE),
  nMolecule_(0),
  allToFirstASU_(false),
  molCentToASU_(false)
{}

// Action_XtalSymm::Help()
void Action_XtalSymm::Help() const {
  mprintf("\t<mask> group <space group> [collect [centroid]]\n"
          "\t[ first | %s ]\n"
          "\t[na <na>] [nb <nb>] [nc <nc>]\n",
          DataSetList::RefArgs);
  mprintf("  Calculate the optimal approach for superimposing symmetry-related subunits\n"
          "  of the simulation back onto one another.  This modifies the coordinate state\n"
          "  for all future actions.\n");
}

// Action_XtalSymm::Init()
Action::RetType Action_XtalSymm::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'xtalsymm' action does not work with > 1 process (%i processes currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  // Get the space group (user must supply this) and fill out the symmetry operations
  spaceGrp_ = actionArgs.GetStringKey("group");
  nCopyA_ = actionArgs.getKeyInt("na", 1);
  nCopyB_ = actionArgs.getKeyInt("nb", 1);
  nCopyC_ = actionArgs.getKeyInt("nc", 1);
  allToFirstASU_ = actionArgs.hasKey("collect");  
  molCentToASU_ = actionArgs.hasKey("centroid");
  // Allocate space to hold all of the symmetry operations.
  // No space group has more than 96.
  R_.assign(   96 * nCopyA_ * nCopyB_ * nCopyC_, Matrix_3x3());
  T_.assign(   96 * nCopyA_ * nCopyB_ * nCopyC_, Vec3());
  RefT_.assign(96 * nCopyA_ * nCopyB_ * nCopyC_, Vec3());
  // Load symmetry operations for space group
  SpaceGroup SG;
  sgID_ = SG.ID( spaceGrp_ );
  if (sgID_ == -1) {
    mprinterr("Error: Could not ID space group string '%s'\n", spaceGrp_.c_str());
    return Action::ERR;
  }
  //mprintf("DEBUG: space group '%s', ID= %i\n", spaceGrp_.c_str(), sgID_);
  nops_ = SG.LoadSymmOps(nCopyA_, nCopyB_, nCopyC_, R_, T_);
  //mprintf("DEBUG: nops= %i\n", nops_);
  if (nops_ < 1) {
    mprinterr("Internal Error: No operations for space group '%s', ID %i\n", spaceGrp_.c_str(), sgID_);
    return Action::ERR;
  }
  // Prepare a table of transposed (inverse) rotation matrices
  Rinv_.clear();
  Rinv_.reserve( 96 * nCopyA_ * nCopyB_ * nCopyC_ );
  for (int i = 0; i < nops_; i++)
    Rinv_.push_back( R_[i].Transposed() );

  // Get the reference frame if possible (if there is no reference, the first
  // frame of the trajectory will fill in for it much later, in DoAction)
  refType_ = NONE;
  if (actionArgs.hasKey("first"))
    refType_ = FIRST;
  else {
    REF_ = init.DSL().GetReferenceFrame( actionArgs );
    if (REF_.error()) return Action::ERR;
    if (REF_.empty()) {
      mprintf("Warning: No reference structure specified. Defaulting to first.\n");
      refType_ = FIRST;
    } else
      refType_ = SPECIFIED;
  }
  
  // Set the masks for all symmetry-related subunits
  Masks_.assign(nops_, AtomMask());
  subunitOpID_.reserve(nops_);
  std::string mask = actionArgs.GetMaskNext();
  if (mask.empty()) {
    mprintf("Error.  A mask for the asymmetric unit must be specified.\n");
    return Action::ERR;
  }
  Masks_[0].SetMaskString(mask);
  // Allocate memory for subunits
  other_.assign( nops_, Frame() );

  mprintf("    XTALSYMM: Mask is %s\n", Masks_[0].MaskString());
  if (refType_ == FIRST)
    mprintf("\tReference is first frame.\n");
  else
    mprintf("\tReference is %s\n", REF_.refName());
  mprintf("\tSpace group '%s'\n", spaceGrp_.c_str());
  mprintf("\tNA= %8i  NB= %8i  NC= %8i  Nops= %8i\n", nCopyA_, nCopyB_, nCopyC_, nops_);
  if (allToFirstASU_) {
    mprintf("\tReimaging");
    if (molCentToASU_)
      mprintf(" by molecule\n");
    else
      mprintf(" by atom\n");
  }

  return Action::OK;
}

//  Action_XtalSymm::BestOrigin()
/** Find the best location for the origin at which to make the rotation that will take one
  * asymmetric unit back onto the original.  Return the result as a three-element vector.
  * The two frames must have already been translated as expected relative to one another for
  * this to work.
  * \param orig   The original frame
  * \param othr   The other frames to superimpose.  This is an array because the problem may need
  *               to simultaneously deal with multiple operations at once.
  * \param operID The list of IDs of the rotation for each symmetry operation.
  */
Vec3 Action_XtalSymm::BestOrigin(Frame& orig, Frame* othr, std::vector<int>& operID)
const
{
  int stride = 3 * orig.Natom();
  
  Matrix<double> A;
  std::vector<double> b;
  A.resize(3, stride * operID.capacity());
  b.resize(stride * operID.capacity());
  int i, j;
  for (i = 0; i < (int)operID.capacity(); i++) {
    int opIDi = operID[i];
    for (j = 0; j < orig.Natom(); j++) {
      A.setElement(0, i*stride + 3*j    , 1.0 - Rinv_[opIDi][0]);
      A.setElement(1, i*stride + 3*j    ,     - Rinv_[opIDi][1]);
      A.setElement(2, i*stride + 3*j    ,     - Rinv_[opIDi][2]);
      A.setElement(0, i*stride + 3*j + 1,     - Rinv_[opIDi][3]);
      A.setElement(1, i*stride + 3*j + 1, 1.0 - Rinv_[opIDi][4]);
      A.setElement(2, i*stride + 3*j + 1,     - Rinv_[opIDi][5]);
      A.setElement(0, i*stride + 3*j + 2,     - Rinv_[opIDi][6]);
      A.setElement(1, i*stride + 3*j + 2,     - Rinv_[opIDi][7]);
      A.setElement(2, i*stride + 3*j + 2, 1.0 - Rinv_[opIDi][8]);

      // CHECK
      if (opIDi < 0 || opIDi >= nops_) {
        printf("opIDi = %d out of %d ops\n", opIDi, nops_);
      }
      if (j >= othr[i].Natom() || j < 0) {
        printf("j = %d out of %d atoms, i = %d out of %zu capacity.\n", j, othr[i].Natom(),
        i, operID.capacity());
      }
      // END CHECK
      
      b[i*stride + 3*j    ] = (*orig.CRD(3*j    )) -
                              (Rinv_[opIDi][0] * (*othr[i].CRD(3*j    )) +
                               Rinv_[opIDi][1] * (*othr[i].CRD(3*j + 1)) +
                               Rinv_[opIDi][2] * (*othr[i].CRD(3*j + 2)));
      b[i*stride + 3*j + 1] = (*orig.CRD(3*j + 1)) -
                              (Rinv_[opIDi][3] * (*othr[i].CRD(3*j    )) +
                               Rinv_[opIDi][4] * (*othr[i].CRD(3*j + 1)) +
                               Rinv_[opIDi][5] * (*othr[i].CRD(3*j + 2)));
      b[i*stride + 3*j + 2] = (*orig.CRD(3*j + 2)) -
                              (Rinv_[opIDi][6] * (*othr[i].CRD(3*j    )) +
                               Rinv_[opIDi][7] * (*othr[i].CRD(3*j + 1)) +
                               Rinv_[opIDi][8] * (*othr[i].CRD(3*j + 2)));
    }
  }
  A.LinearLeastSquares(b.data());

  return Vec3(b[0], b[1], b[2]);
}

// Action_XtalSymm::BestSuperposition()
/** Take a symmetry operation and find how well it could possibly make a superposition of the
  * given asymmetric unit back onto the original.  Return the result as mass-weighted atom
  * positional RMSD plus a three-element vector describing the optimal displacement between
  * centers of mass for the two subunits.
  * \param  maskID The number of the mask for the asymmetric unit to superimpose, >= 1
  * \param  operID The number of the operation to try, >= 1
  * \param  leads  A growing array of leads that will ultimately produce an origin for the
  *                transformation of the entire system.
  * \param  nLead  The number of leads found thus far.
  */
void Action_XtalSymm::BestSuperposition(int maskID, int operID, XtalDock* leads, int& nLead)
const
{
  // Set up frames for the original and symmetry-related subunits
  Frame orig, othr;
  orig = Frame(Masks_[0].Nselected());
  othr = Frame(Masks_[maskID].Nselected());
  orig.SetCoordinates(RefFrame_, Masks_[0]);
  othr.SetCoordinates(RefFrame_, Masks_[maskID]);
  Vec3 corig = orig.VCenterOfMass(0, orig.Natom());
  Vec3 cothr = othr.VCenterOfMass(0, othr.Natom());
  Vec3 cdiff = cothr - corig;
  Matrix_3x3 U, invU;
  U = RefFrame_.BoxCrd().UnitCell(1.0);
  RefFrame_.BoxCrd().ToRecip(U, invU);
  cdiff = invU * cdiff;
  const double nearTwo = 1.999999999;  
  double minx = floor(cdiff[0] - nearTwo);
  double miny = floor(cdiff[1] - nearTwo);
  double minz = floor(cdiff[2] - nearTwo);
  double dx, dy, dz, di, dj, dk;
  std::vector<int> operIDv(1);
  operIDv[0] = operID;
  for (di = 0.0; di <= 4.0; di++) {
    for (dj = 0.0; dj <= 4.0; dj++) {
      for (dk = 0.0; dk <= 4.0; dk++) {
        dx = minx + di;
        dy = miny + dj;
        dz = minz + dk;
        
        // Translate to test a new starting location, then reverse the symmetry operation.
        // This involves computation of a new origin that will take both subunits into a
        // unique frame of reference.
        Vec3 Tvec = Vec3(-dx - T_[operID][0], -dy - T_[operID][1], -dz - T_[operID][2]);
        Tvec = U * Tvec;
        othr.SetCoordinates(RefFrame_, Masks_[maskID]);
        othr.Translate(Tvec);
        Vec3 Ovec = BestOrigin(orig, &othr, operIDv);
        orig.NegTranslate(Ovec);
        othr.NegTranslate(Ovec);
        othr.Rotate(Rinv_[operID]);
        
        // Compute the root-mean squared deviation of the two subunits, then
        // put the original subunit back in its official frame of reference.
        orig.Translate(Ovec);
        othr.Translate(Ovec);
        leads[nLead].subunit_ = maskID;
        leads[nLead].opID_    = operID;
        leads[nLead].rmsd_    = orig.RMSD_NoFit(othr, false);
        leads[nLead].origin_  = Ovec;
        Tvec = Vec3(cdiff[0] - dx, cdiff[1] - dy, cdiff[2] - dz);
        leads[nLead].displc_  = Tvec;
        if (leads[nLead].rmsd_ < 10.0) {
          nLead++;
        }
      }
    }
  }
}

// Action_XtalSymm::DetectAsuResidence
/** Detect which ASU the given point lies in given the fractional coordinates of the point.
  */
Action_XtalSymm::TransOp Action_XtalSymm::DetectAsuResidence(double x, double y, double z)
const
{
  int i, j, k, m;
  Vec3 pt;
  TransOp amove;
  
  for (i = -1; i <= 1; i++) {
    for (j = -1; j <= 1; j++) {
      for (k = -1; k <= 1; k++) {
        for (m = 0; m < nops_; m++) {
          pt = Vec3(x + (double)i, y + (double)j, z + (double)k);
          pt = pt - T_[m];
          pt = Rinv_[m] * pt;
          if (PointInPrimaryASU(pt[0], pt[1], pt[2])) {

            // This is a solution, break out of all loops
            amove.opID_ = m;
            amove.tr_x_ = (double)i;
            amove.tr_y_ = (double)j;
            amove.tr_z_ = (double)k;
            return amove;
          }
        }
      }
    }
  }

  // Raise a warning if this didn't fall into any ASU
  amove.opID_ = 0;
  amove.tr_x_ = 0.0;
  amove.tr_y_ = 0.0;
  amove.tr_z_ = 0.0;
  mprintf("Warning: point %9.4lf %9.4lf %9.4lf did not fall into any asymmteric unit.\n",
          x, y, z);

  return amove;
}

// Action_XtalSymm::FindPrevious()
void Action_XtalSymm::FindPrevious(int& prevOpID, double& prevTx, double& prevTy, double& prevTz,
                                   double ptx, double pty, double ptz) const
{
  for (int ii = -1; ii <= 1; ii++) {
    for (int jj = -1; jj <= 1; jj++) {
      for (int kk = -1; kk <= 1; kk++) {
        for (int m = 0; m < nops_; m++) {
          Vec3 pt(ptx + (double)ii, pty + (double)jj, ptz + (double)kk);
          pt = pt - T_[m];
          pt = Rinv_[m] * pt;
          if (PointInPrimaryASU(pt[0], pt[1], pt[2])) {

            // This is a solution, break out of all loops
            prevOpID = m;
            prevTx = (double)ii;
            prevTy = (double)jj;
            prevTz = (double)kk;
            return;
          }
        }
      }
    }
  }
}

// Action_XtalSymm::BuildAsuGrid()
/** Build a grid spanning the unit cell to indicate the approximate extent of
  * each asymmetric unit's volume.  Grid bins that do not fall entirely within
  * one asymmetric unit will be labelled as wildcards and any coordinates that
  * fall in those bins will have to be checked against all asymmetric units.
  */
void Action_XtalSymm::BuildAsuGrid()
{
  int i, j, k, ii, jj, kk;

  int prevOpID = 0;
  double prevTx = 0.0;
  double prevTy = 0.0;
  double prevTz = 0.0;
  AsuGrid_.assign(IASU_GRID_BINS_ * IASU_GRID_BINS_ * IASU_GRID_BINS_, TransOp());
  for (i = 0; i < IASU_GRID_BINS_; i++) {
    for (j = 0; j < IASU_GRID_BINS_; j++) {
      for (k = 0; k < IASU_GRID_BINS_; k++) {

        // First, select a point in the middle of the grid bin and find its ASU.
        double ptx = ((double)i + 0.5) / DASU_GRID_BINS_;
        double pty = ((double)j + 0.5) / DASU_GRID_BINS_;
        double ptz = ((double)k + 0.5) / DASU_GRID_BINS_;
        Vec3 pt = Vec3(ptx + prevTx, pty + prevTy, ptz + prevTz);
        pt = pt - T_[prevOpID];
        pt = Rinv_[prevOpID] * pt;
        bool success = PointInPrimaryASU(pt[0], pt[1], pt[2]);
        if (!success) {
          FindPrevious(prevOpID, prevTx, prevTy, prevTz, ptx, pty, ptz);
        }
        bool complete = true;
        for (ii = 0; ii < 2; ii++) {
          for (jj = 0; jj < 2; jj++) {
            for (kk = 0; kk < 2; kk++) {
              pt = Vec3(ptx + (double)ii/DASU_GRID_BINS_ + prevTx,
                        pty + (double)jj/DASU_GRID_BINS_ + prevTy,
                        ptz + (double)kk/DASU_GRID_BINS_ + prevTz);
              pt = pt - T_[prevOpID];
              pt = Rinv_[prevOpID] * pt;
              if (PointInPrimaryASU(pt[0], pt[1], pt[2]) == false) {
                complete = false;
                ii = 2;
                jj = 2;
                kk = 2;
              }
            }
          }
        }

        // Record the result
        int idx = (i*IASU_GRID_BINS_ + j)*IASU_GRID_BINS_ + k;
        if (complete) {
          AsuGrid_[idx].opID_ = prevOpID;
          AsuGrid_[idx].tr_x_ = prevTx;
          AsuGrid_[idx].tr_y_ = prevTy;
          AsuGrid_[idx].tr_z_ = prevTz;
        }
        else {
          AsuGrid_[idx].opID_ = -1;
          AsuGrid_[idx].tr_x_ = 0.0;
          AsuGrid_[idx].tr_y_ = 0.0;
          AsuGrid_[idx].tr_z_ = 0.0;
        }
      }
    }
  }
}

// Action_XtalSymm::Setup()
Action::RetType Action_XtalSymm::Setup(ActionSetup& setup)
{
  int i, j, k;
  
  // Checks on the sanity of the system and suitability for crystal symmetry operations
  if (setup.CoordInfo().HasCrd() == false || setup.CoordInfo().HasBox() == false) {
    return Action::ERR;
  }
  
  // Set up the integer mask for the first, required, mask string
  if (setup.Top().SetupIntegerMask(Masks_[0])) {
    return Action::ERR;
  }

  // Set up the mask for the entire topology, for making the RefFrame_ clone of it later
  std::string str;
  str.assign(":*");
  tgtMask_.SetMaskString(str);
  if (setup.Top().SetupIntegerMask(tgtMask_)) {
    return Action::ERR;
  }
  if (refType_ == SPECIFIED) {
    // The way masks are currently used, the reference topology must match
    // current topology. For now just check same # atoms.
    if (REF_.Parm().Natom() != setup.Top().Natom())
      mprintf("Warning: Reference '%s' # atoms %i does not match current # atoms %i.\n",
              REF_.refName(), REF_.Parm().Natom(), setup.Top().Natom());

    RefFrame_.SetupFrame( REF_.Parm().Natom() );
    RefFrame_.SetCoordinates( REF_.Coord(), tgtMask_ );
  }

  // If there are not enough masks specified (it's tedious to do, and would make
  // the action command lengthy), then generate all the rest based on the first.
  int nMaskAtom = Masks_[0].Nselected();
  int nTopolAtom = setup.Top().Natom();
  std::vector<int> baseMask = Masks_[0].Selected();
  std::vector<int> newMask = Masks_[0].Selected();

  // Create an array to hold whether each atom has been included in a mask
  std::vector<int> occupancy(nTopolAtom, 0);
  for (i = 0; i < nMaskAtom; i++) {
    occupancy[baseMask[i]] = 1;
  }

  // For every other mask that is needed, loop over the topology and find
  // equivalent atoms.  If all atoms can be matched in a sequence and spacing
  // identical to the original mask, without stepping on atoms that are already
  // taken, then this is a new, equivalent, asymmetric unit.
  int startpos = 0;
  int maskwidth = 0;
  for (i = 0; i < nMaskAtom; i++) {
    if (baseMask[i] - baseMask[0] > maskwidth) {
      maskwidth = baseMask[i] - baseMask[0];
    }
  }
  for (i = 1; i < nops_; i++) {
    for (j = startpos; j < nTopolAtom - maskwidth; j++) {
      int nmatched = 0;
      for (k = 0; k < nMaskAtom; k++) {
        int basepos = baseMask[k];
        int candpos = j + baseMask[k] - baseMask[0];
        if (occupancy[candpos] == 0 && setup.Top().TruncResNameAtomName(candpos) ==
                                       setup.Top().TruncResNameAtomName(basepos)) {
          nmatched++;
        }
        else {
          break;
        }
      }
        
      // If there are as many matches as there are atoms in the mask, this was a success.
      // Otherwise, increment the starting position so that future searches do not run
      // back over the same atoms again.
      if (nmatched == nMaskAtom) {
        Masks_[i].ClearSelected();
        for (k = 0; k < nMaskAtom; k++) {
          newMask[k] = j + baseMask[k] - baseMask[0];
          occupancy[j + baseMask[k] - baseMask[0]] = 1;
        }
        Masks_[i].AddAtoms(newMask);
        break;
      }
      else {
        startpos++;
      }
    }
  }

  // Determine which rotations are identity matrices
  rotIdentity_.clear();
  rotIdentity_.reserve( nops_ );
  for (i = 0; i < nops_; i++) {
    if (fabs(R_[i].Row1()[0] - 1.0) < 1.0e-6 && fabs(R_[i].Row2()[1] - 1.0) < 1.0e-6 &&
        fabs(R_[i].Row3()[2] - 1.0) < 1.0e-6 && fabs(R_[i].Row1()[1]) < 1.0e-6 &&
        fabs(R_[i].Row1()[2]) < 1.0e-6 && fabs(R_[i].Row2()[0]) < 1.0e-6 &&
        fabs(R_[i].Row2()[2]) < 1.0e-6 && fabs(R_[i].Row3()[0]) < 1.0e-6 &&
        fabs(R_[i].Row3()[1]) < 1.0e-6) {
      rotIdentity_.push_back( true );
    }
    else {
      rotIdentity_.push_back( false );
    }
  }
  
  // Create a lookup table to help determine ASUs for loose particles,
  // then a mask to detail all particles not in one of the designated
  // asymmetric units.  These 'loose' particles, which are probably
  // solvent, do not have an assigned ASU but nonetheless must fall
  // into one at any given time.
  if (allToFirstASU_) {
    BuildAsuGrid();
    std::vector<int> LoneAtoms( setup.Top().Natom() );
    std::vector<int> MoleAtoms( setup.Top().Natom() );
    for (i = 0; i < setup.Top().Natom(); i++) {
      LoneAtoms[i] = 1;
      MoleAtoms[i] = 0;
    }
    for (i = 0; i < nops_; i++) {
      for (j = 0; j < Masks_[i].Nselected(); j++) {
        LoneAtoms[Masks_[i].Selected()[j]] = 0;
      }
    }
    int nnonasu = 0;
    for (i = 0; i < setup.Top().Natom(); i++) {
      if (LoneAtoms[i] == 1) {
        nnonasu++;
      }
    }
    if (molCentToASU_) {

      // If solvent molecules are to be re-imaged whole into the primary ASU, mark them
      // separately and remove them from the list of individual atoms to remove.
      nMolecule_ = setup.Top().Nmol();
      molLimits_.clear();
      molLimits_.reserve(2 * nMolecule_);
      molInSolvent_.assign(nMolecule_, true);
      for (i = 0; i < nMolecule_; i++) {
        molLimits_.push_back( setup.Top().Mol(i).BeginAtom() );
        molLimits_.push_back( setup.Top().Mol(i).EndAtom()   );
        for (j = molLimits_[2*i]; j < molLimits_[2*i + 1]; j++) {
          if (LoneAtoms[j] == 0) {
            molInSolvent_[i] = false;
          }
        }
        if (molInSolvent_[i]) {
          for (j = molLimits_[2*i]; j < molLimits_[2*i + 1]; j++) {
            LoneAtoms[j] = 0;
            MoleAtoms[j] = 1;
          }
        }
      }
    }
    else {
      nMolecule_ = 0;
    }

    std::vector<int> SolventList;
    SolventList.reserve(nnonasu);
    for (i = 0; i < setup.Top().Natom(); i++) {
      if (LoneAtoms[i] == 1) {
        SolventList.push_back(i);
      }
    }
    SolventParticles_ = AtomMask(SolventList, nnonasu);
    if (molCentToASU_) {
      SolventList.clear();
      SolventList.reserve(nnonasu);
      for (i = 0; i < setup.Top().Natom(); i++) {
        if (MoleAtoms[i] == 1) {
          SolventList.push_back(i);
        }
      }
      SolventMolecules_ = AtomMask(SolventList, nnonasu);
    }
  }

  // Allocate space for subunit frames
  for (i = 0; i < nops_; i++)
    other_[i].SetupFrame(Masks_[i].Nselected());

  return Action::OK;
}

// Action_XtalSymm::OperationAvailable()
/** Test whether a given symmetry operation is available for use in an approach to reconstruct
  * the unit cell.
  * \param leads:          The list of leads, each specifying an operation that will take one subunit
                           back onto the original subunit given a properly imaged displacement and
                           an origin.
  * \param HowToGetThere:  List of leads accumulated thus far.
  * \param ncurr:          The position to add the candidate lead to the list.
  */
bool Action_XtalSymm::OperationAvailable(XtalDock* leads, std::vector<int> const& HowToGetThere, int ncurr)
{
  int i;

  for (i = 0; i < ncurr; i++) {
    if (leads[HowToGetThere[i]].opID_ == leads[HowToGetThere[ncurr]].opID_) {
      return false;
    }
  }

  return true;
}

// Action_XtalSymm::OriginsAlign()
/** Test whether the origin of a particular lead will work in the context of the others
  * accumulated thus far.
  */
bool Action_XtalSymm::OriginsAlign(XtalDock* leads, std::vector<int> const& HowToGetThere, int ncurr)
const
{
  int i;
  double origx, origy, origz;
  
  // First, find a symmetry operation that involves an
  // actual rotation, where the origin would matter
  bool oxfound = false, oyfound = false, ozfound = false;
  for (i = 0; i <= ncurr; i++) {
    int thisSU = leads[HowToGetThere[i]].subunit_;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    if (fabs(R_[thisSU].Row1()[0] - 1.0) >= 1.0e-6 || fabs(R_[thisSU].Row2()[0]) >= 1.0e-6 ||
        fabs(R_[thisSU].Row3()[0]) >= 1.0e-6) {
      if (oxfound == false) {
        origx = leads[HowToGetThere[i]].origin_[0];
        oxfound = true;
      }
      else {
        dx = origx - leads[HowToGetThere[i]].origin_[0];
      }
    }
    if (fabs(R_[thisSU].Row1()[1]) >= 1.0e-6 || fabs(R_[thisSU].Row2()[1] - 1.0) >= 1.0e-6 ||
        fabs(R_[thisSU].Row3()[1]) >= 1.0e-6) {
      if (oyfound == false) {
        origy = leads[HowToGetThere[i]].origin_[1];
        oyfound = true;
      }
      else {
        dy = origy - leads[HowToGetThere[i]].origin_[1];
      }
    }
    if (fabs(R_[thisSU].Row1()[2]) >= 1.0e-6 || fabs(R_[thisSU].Row2()[2]) >= 1.0e-6 ||
        fabs(R_[thisSU].Row3()[2] - 1.0) >= 1.0e-6) {
      if (ozfound == false) {
        origz = leads[HowToGetThere[i]].origin_[2];
        ozfound = true;
      }
      else {
        dz = origz - leads[HowToGetThere[i]].origin_[2];
      }
    }
    if (dx*dx + dy*dy + dz*dz >= 100.0) {
      return false;
    }
  }

  return true;
}

// Action_XtalSymm::DoAction()
Action::RetType Action_XtalSymm::DoAction(int frameNum, ActionFrame& frm)
{
  int i, j;
  Frame orig;
  Matrix_3x3 U, invU;

  // Allocate space for the subunit frames
  orig = Frame(Masks_[0].Nselected());

  // Determine the optimal strategy for superimposing subunits.  This has to be
  // done here, rather than in Setup, because the reference frame for determining
  // the strategy may have to be the first frame.
  if (refType_ != NONE) {
    // Use the reference if supplied.  Otherwise use the first frame.
    if (refType_ == FIRST)
    {
      RefFrame_.SetupFrame(frm.Frm().Natom());
      RefFrame_.SetCoordinates(frm.Frm(), tgtMask_);
    }
    // Ensure this is only done once FIXME ok?
    refType_ = NONE;

    U = RefFrame_.BoxCrd().UnitCell(1.0);
    RefFrame_.BoxCrd().ToRecip(U, invU);
    XtalDock* leads = new XtalDock[nops_ * nops_ * 125]; // TODO vector?
    int nLead = 0;
    for (i = 0; i < nops_; i++) {
      for (j = 0; j < nops_; j++) {
        BestSuperposition(i, j, leads, nLead);
      }
    }

    // Obtain clusters of origins, then see whether any of them holds all symmetry operations
    // and subunits, which of that group contains the tightest cluster of origins, and which
    // contains the lowest overall atom positional RMSD.  
    int fnidXfrm = -1;
    for (i = 0; i < nops_; i++) {
      if (rotIdentity_[i] == false) {
        fnidXfrm = i;
        break;
      }
    }
    if (fnidXfrm == -1) {

      // This was a P1 crystal.  Any origin will work for applying the symmetry operations.
      for (i = 0; i < nops_; i++) {
        for (j = 0; j < nLead; j++) {
          if (leads[j].subunit_ == i) {
            subunitOpID_[i] = leads[j].opID_;
            RefT_[i] = leads[j].displc_;
            break;
          }
        }
      }
    }
    else {
      
      // This is not a P1 crystal, and there are many possible combinations of the discovered
      // leads that could paint the correct picture of how to reassemble the unit (super)
      // cell.
      std::vector<int> HowToGetThere(nops_, 0);
      double bestRmsd = 1.0e8;
      double bestOrig = 1.0e8;
      i = 0;
      while (HowToGetThere[0] < nLead) {
        while (i < nops_) {
          while (HowToGetThere[0] < nLead &&
                 (leads[HowToGetThere[i]].subunit_ != i ||
                  OperationAvailable(leads, HowToGetThere, i) == false ||
                  OriginsAlign(leads, HowToGetThere, i) == false)) {
            HowToGetThere[i] += 1;
            while (HowToGetThere[i] == nLead && i > 0) {
              HowToGetThere[i] = 0;
              i--;
              HowToGetThere[i] += 1;
            }
          }
          i++;
          if (i == nops_) {
            
            // Check the RMSD that would result from this situation
            orig.SetCoordinates(RefFrame_, Masks_[0]);
            Vec3 corig = orig.VCenterOfMass(0, orig.Natom());
            std::vector<int> trialOpID;
            trialOpID.reserve(nops_);
            for (j = 0; j < nops_; j++) {
              other_[j].SetCoordinates(RefFrame_, Masks_[leads[HowToGetThere[j]].subunit_]);
              Vec3 cothr = other_[j].VCenterOfMass(0, other_[j].Natom());
              Vec3 cdiff = cothr - corig;
              cdiff = invU * cdiff;
              Vec3 cmove = leads[HowToGetThere[j]].displc_ - cdiff;
              cmove[0] = round(cmove[0]);
              cmove[1] = round(cmove[1]);
              cmove[2] = round(cmove[2]);
              cmove = (U * cmove) - (U * T_[leads[HowToGetThere[j]].opID_]);
              other_[j].Translate(cmove);
              trialOpID[j] = leads[HowToGetThere[j]].opID_;
            }
            Vec3 trOvec = BestOrigin(orig, &other_[0], trialOpID);
            orig.NegTranslate(trOvec);
            for (j = 0; j < nops_; j++) {
              other_[j].NegTranslate(trOvec);
              other_[j].Rotate(Rinv_[trialOpID[j]]);
            }
            double trmsd = 0.0;
            for (j = 0; j < nops_; j++) {
              trmsd += orig.RMSD_NoFit(other_[j], false);
            }
            if (trmsd < bestRmsd + 0.1) {
              double torig = 0.0;
              for (j = 0; j < nops_; j++) {
                Vec3 dsp = U * leads[HowToGetThere[j]].displc_;
                torig += dsp[0]*dsp[0] + dsp[1]*dsp[1] + dsp[2]*dsp[2];
              }
              torig += trOvec[0]*trOvec[0] + trOvec[1]*trOvec[1]+ trOvec[2]*trOvec[2];
              if (torig < bestOrig) {
                bestRmsd = trmsd;
                bestOrig = torig;
                for (j = 0; j < nops_; j++) {
                  subunitOpID_[j] = trialOpID[j];
                  RefT_[j] = leads[HowToGetThere[j]].displc_;
                }
              }
            }
      
            // Exit if the standard solution works well enough
            if (bestRmsd < 1.0) {
              bool stdworks = true;
              for (j = 0; j < nops_; j++) {
                if (trialOpID[j] != j) {
                  stdworks = false;
                }
              }
              if (stdworks) {
                HowToGetThere[0] = nLead;
              }
            }
      
            // If there is more to do, decrement i and keep on going
            if (HowToGetThere[0] < nLead && i > 0) {
              i--;
              HowToGetThere[i] += 1;
              while (HowToGetThere[i] == nLead && i > 0) {
                HowToGetThere[i] = 0;
                i--;
                HowToGetThere[i] += 1;
              }
            }
          }
        }
      }
    }
          
    // Free allocated memory
    delete[] leads;
  }
  
  // Loop over all subunits, set them according to the correct displacements from
  // the original subunits (imaging considerations), and apply the transformations.
  U = frm.Frm().BoxCrd().UnitCell(1.0);
  frm.Frm().BoxCrd().ToRecip(U, invU);
  orig = Frame(Masks_[0].Nselected());

  orig.SetCoordinates(frm.Frm(), Masks_[0]);
  for (i = 0; i < nops_; i++) {
    int opID = subunitOpID_[i];
    
    // Get each subunit and the box transformations
    other_[i].SetCoordinates(frm.Frm(), Masks_[i]);

    // Get the displacement between the two subunits' centers of mass
    Vec3 corig = orig.VCenterOfMass(0, orig.Natom());
    Vec3 cothr = other_[i].VCenterOfMass(0, other_[i].Natom());
    Vec3 cdiff = cothr - corig;
    cdiff = invU * cdiff;

    // How much does the displacment need to change in order to set the second
    // subunit in the correct place for the symmetry operation to make any sense?
    Vec3 cmove = RefT_[opID] - cdiff;
    cmove[0] = round(cmove[0]);
    cmove[1] = round(cmove[1]);
    cmove[2] = round(cmove[2]);
    cmove = (U * cmove) - (U * T_[opID]);
    
    // Translate the second subunit as needed.  This is the first part of the
    // reverse symmetry operation.  After all subunits have had their translations
    // removed, a consensus origin for the whole system will be computed in order
    // to apply the reverse rotations.
    other_[i].Translate(cmove);
    frm.ModifyFrm().Translate(cmove, Masks_[i]);
  }
  Vec3 Ovec = BestOrigin(orig, &other_[0], subunitOpID_);

  // Apply the final result, the origin at which all of these transformations are valid
  frm.ModifyFrm().NegTranslate(Ovec);  
  for (i = 0; i < nops_; i++) {
    int opID = subunitOpID_[i];
    frm.ModifyFrm().Rotate(Rinv_[opID], Masks_[i]);
  }

  // It is now possible to re-image all solvent molecules (those not in one of the
  // specified asymmetric units), and then determine to which asymmetric unit they
  // belong.
  if (allToFirstASU_) {
    frm.ModifyFrm().Rotate(invU, SolventParticles_);
    if (molCentToASU_) {
      frm.ModifyFrm().Rotate(invU, SolventMolecules_);
      for (i = 0; i < nMolecule_; i++) {
        if (molInSolvent_[i]) {

          // Re-image the entire molecule
          double x = 0.0;
          double y = 0.0;
          double z = 0.0;
          for (j = molLimits_[2*i]; j < molLimits_[2*i + 1]; j++) {
            x += *frm.Frm().CRD(3*j);
            y += *frm.Frm().CRD(3*j + 1);
            z += *frm.Frm().CRD(3*j + 2);
          }
          double dfac = 1.0 / (double)(molLimits_[2*i + 1] - molLimits_[2*i]);
          x *= dfac;
          y *= dfac;
          z *= dfac;
          frm.ModifyFrm().Translate(Vec3(0.5 - round(x), 0.5 - round(y), 0.5 - round(z)),
                                    molLimits_[2*i], molLimits_[2*i + 1]);
          x += 0.5 - round(x);
          y += 0.5 - round(y);
          z += 0.5 - round(z);
          
          // Use the grid to determine the asymmetric unit (if the grid says "operation -1",
          // an intensive search will be done to find the correct ASU)
          int gidx = (int)(x * DASU_GRID_BINS_);
          int gidy = (int)(y * DASU_GRID_BINS_);
          int gidz = (int)(z * DASU_GRID_BINS_);
          gidx = (gidx*IASU_GRID_BINS_ + gidy)*IASU_GRID_BINS_ + gidz;
          TransOp Vm = (AsuGrid_[gidx].opID_ == -1) ? DetectAsuResidence(x, y, z) :
                                                    AsuGrid_[gidx];
          frm.ModifyFrm().Translate(Vec3(Vm.tr_x_, Vm.tr_y_, Vm.tr_z_) - T_[Vm.opID_],
                                    molLimits_[2*i], molLimits_[2*i + 1]);
          frm.ModifyFrm().Rotate(Rinv_[Vm.opID_], molLimits_[2*i], molLimits_[2*i + 1]);
        }
      }
      frm.ModifyFrm().Rotate(U, SolventMolecules_);
    }
    for (i = 0; i < SolventParticles_.Nselected(); i++) {
      int iatm = SolventParticles_.Selected()[i];

      // Re-image the ith solvent atom
      double x = *frm.Frm().CRD(3*iatm);
      double y = *frm.Frm().CRD(3*iatm + 1);
      double z = *frm.Frm().CRD(3*iatm + 2);
      frm.ModifyFrm().Translate(Vec3(0.5 - round(x), 0.5 - round(y), 0.5 - round(z)), iatm);
      x += 0.5 - round(x);
      y += 0.5 - round(y);
      z += 0.5 - round(z);

      // Use the grid to determine the asymmetric unit (if the grid says "operation -1",
      // an intensive search will be done to find the correct ASU)
      int gidx = (int)(x * DASU_GRID_BINS_);
      int gidy = (int)(y * DASU_GRID_BINS_);
      int gidz = (int)(z * DASU_GRID_BINS_);
      gidx = (gidx*IASU_GRID_BINS_ + gidy)*IASU_GRID_BINS_ + gidz;
      TransOp Vm = (AsuGrid_[gidx].opID_ == -1) ? DetectAsuResidence(x, y, z) : AsuGrid_[gidx];
      frm.ModifyFrm().Translate(Vec3(Vm.tr_x_, Vm.tr_y_, Vm.tr_z_) - T_[Vm.opID_], iatm);
      frm.ModifyFrm().Rotate(Rinv_[Vm.opID_], iatm);
    }
    frm.ModifyFrm().Rotate(U, SolventParticles_);
  }
  
  return Action::MODIFY_COORDS;
}

//---------------------------------------------------------------------------------------------
// dmin: Various functions for overloading the min function from the std library
//---------------------------------------------------------------------------------------------
double Action_XtalSymm::dmin(double a, double b)
{
  if (a < b) {
    return a;
  }
  else {
    return b;
  }
}

double Action_XtalSymm::dmin(double a, double b, double c)
{
  if (a < b) {
    if (a < c) {
      return a;
    }
    else {
      return c;
    }
  }
  else {
    if (b < c) {
      return b;
    }
    else {
      return c;
    }
  }
}

double Action_XtalSymm::dmin(double a, double b, double c, double d)
{
  double ab = dmin(a, b);
  double cd = dmin(c, d);

  return dmin(ab, cd);
}

//---------------------------------------------------------------------------------------------
// dmax: Various functions for overloading the min function from the std library
//---------------------------------------------------------------------------------------------
double Action_XtalSymm::dmax(double a, double b)
{
  if (a > b) {
    return a;
  }
  else {
    return b;
  }
}

double Action_XtalSymm::dmax(double a, double b, double c)
{
  if (a > b) {
    if (a > c) {
      return a;
    }
    else {
      return c;
    }
  }
  else {
    if (b > c) {
      return b;
    }
    else {
      return c;
    }
  }
}

double Action_XtalSymm::dmax(double a, double b, double c, double d)
{
  double ab = dmax(a, b);
  double cd = dmax(c, d);

  return dmax(ab, cd);
}

// Action_XtalSymm::PointInPrimaryASU()
/** Given the fractional coordinates of a point, determine whether it lies 
  * within the primary asymmetric unit as defined by the space group's geometry.
  */
bool Action_XtalSymm::PointInPrimaryASU(double x, double y, double z)
const
{
  const double twothr  = 0.66666666666667;
  const double half    = 0.5;
  const double threig  = 0.375;
  const double third   = 0.33333333333333;
  const double fourth  = 0.25;
  const double sixth   = 0.16666666666667;
  const double eighth  = 0.125;
  const double twelfth = 0.83333333333333;
  bool result = true;

  // Scale the coordinates by the number of times the unit cell has been replicated.
  x *= (double)nCopyA_;
  y *= (double)nCopyB_;
  z *= (double)nCopyC_;

  // Giant case switch to test each space group
  switch(sgID_) {
    case 0:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > 1.0) return false;
      break;
    case 1:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > 1.0) return false;
      break;
    case 2:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > half) return false;
      break;
    case 3:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > half) return false;
      break;
    case 4:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 5:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 6:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 7:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 8:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 9:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 10:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 11:
      if (0.0 > x || x > half || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 12:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > half) return false;
      break;
    case 13:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 14:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 15:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 16:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 17:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 18:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
    case 19:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 20:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 21:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 22:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 23:
      if (0.0 > x || x > fourth || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 24:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 25:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 26:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 27:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 28:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 29:
      if (0.0 > x || x > fourth || 0.0 > y || y > 1.0 || 0.0 > z || z > 1.0) return false;
      break;
    case 30:
      if (0.0 > x || x > fourth || 0.0 > y || y > 1.0 || 0.0 > z || z > 1.0) return false;
      break;
    case 31:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > half) return false;
      break;
    case 32:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 33:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 34:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 35:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 36:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 37:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 38:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 39:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 40:
      if (0.0 > x || x > half || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 41:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 42:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 43:
      if (0.0 > x || x > fourth || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 44:
      if (0.0 > x || x > fourth || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 45:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 46:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 47:
      if (0.0 > x || x > fourth || 0.0 > y || y > 1.0 || 0.0 > z || z > half) return false;
      break;
    case 48:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 49:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 50:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 51:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 52:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 53:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > fourth || 0.0 > z || z > half) return false;
      break;
    case 54:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > fourth) return false;
      break;
    case 55:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 56:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 57:
      if (0.0 > x || x > fourth || 0.0 > y || y > 1.0 || 0.0 > z || z > half) return false;
      break;
    case 58:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > fourth) return false;
      break;
    case 59:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 60:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 61:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 62:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 63:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 64:
      if (0.0 > x || x > half || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 65:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 66:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 67:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 68:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 69:
      if (0.0 > x || x > half || 0.0 > y || y > fourth || 0.0 > z || z > half) return false;
      break;
    case 70:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 71:
      if (0.0 > x || x > fourth || 0.0 > y || y > fourth || 0.0 > z || z > half) return false;
      break;
    case 72:
      if (0.0 > x || x > eighth || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 73:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 74:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 75:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 76:
      if (0.0 > x || x > fourth || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 77:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 78:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 79:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 80:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 81:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 82:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > fourth) return false;
      break;
    case 83:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0) return false;
      break;
    case 84:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 85:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 86:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 87:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 88:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > fourth) return false;
      break;
    case 89:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 90:
      if (0.0 > x || x > fourth || 0.0 > y || y > fourth || 0.0 > z || z > 1.0) return false;
      break;
    case 91:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 92:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 93:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > eighth) return false;
      break;
    case 94:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > eighth) return false;
      break;
    case 95:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > fourth) return false;
      break;
    case 96:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 97:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > eighth) return false;
      break;
    case 98:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > eighth) return false;
      break;
    case 99:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 100:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > eighth) return false;
      break;
    case 101:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0 || x > y) {
        return false;
      }
      break;
    case 102:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0 || y > half-x) {
        return false;
      }
      break;
    case 103:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0 || x > y) {
        return false;
      }
      break;
    case 104:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0 || x > y) {
        return false;
      }
      break;
    case 105:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 106:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 107:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 108:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 109:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || x > y) {
        return false;
      }
      break;
    case 110:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || y > half-x) {
        return false;
      }
      break;
    case 111:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 112:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 113:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0 || x > y) {
        return false;
      }
      break;
    case 114:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 115:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > 1.0 || y > half-x) {
        return false;
      }
      break;
    case 116:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 117:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 118:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > fourth) return false;
      break;
    case 119:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half) return false;
      break;
    case 120:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > fourth) return false;
      break;
    case 121:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 122:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 123:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || x > y) {
        return false;
      }
      break;
    case 124:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > eighth) return false;
      break;
    case 125:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || x > y) {
        return false;
      }
      break;
    case 126:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 127:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || y > half-x) {
        return false;
      }
      break;
    case 128:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 129:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || y > half-x) {
        return false;
      }
      break;
    case 130:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 131:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || y > half-x) {
        return false;
      }
      break;
    case 132:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 133:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 134:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || x > y) {
        return false;
      }
      break;
    case 135:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 136:
      if (0.0 > x || x > half || 0.0 > y || y > 1.0 || 0.0 > z || z > fourth || x > y ||
          y > 1.0-x) {
        return false;
      }
      break;
    case 137:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 138:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || x > y) {
        return false;
      }
      break;
    case 139:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth) return false;
      break;
    case 140:
      if (0.0 > x || x > fourth || 0.0 > y || y > half || 0.0 > z || z > 1.0 || x > y ||
          y > half-x) {
        return false;
      }
      break;
    case 141:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth || x > y) {
        return false;
      }
      break;
    case 142:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth ||
          y > half-x) {
        return false;
      }
      break;
    case 143:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > eighth) return false;
      break;
    case 144:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > eighth) return false;
      break;
    case 145:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > 1.0 ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 146:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > third) return false;
      break;
    case 147:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > third) return false;
      break;
    case 148:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > third ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 149:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 150:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > sixth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 151:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 152:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 153:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > sixth) return false;
      break;
    case 154:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > sixth) return false;
      break;
    case 155:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > sixth) return false;
      break;
    case 156:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > sixth) return false;
      break;
    case 157:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > sixth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 158:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > 1.0 ||
          x > 2.0*y || y > dmin(1.0-x, 2.0*x)) {
        return false;
      }
      break;
    case 159:
      if (0.0 > x || x > twothr || 0.0 > y || y > half || 0.0 > z || z > 1.0 ||
          x > (y+1.0)/2.0 || y > dmin(1.0-x, x)) {
        return false;
      }
      break;
    case 160:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 161:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 162:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > third ||
          x > 2.0*y || y > dmin(1.0-x, 2.0*x)) {
        return false;
      }
      break;
    case 163:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > sixth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 164:
      if (0.0 > x || x > twothr || 0.0 > y || y > half || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, x)) {
        return false;
      }
      break;
    case 165:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > fourth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 166:
      if (0.0 > x || x > twothr || 0.0 > y || y > third || 0.0 > z || z > 1.0 ||
          x > (1.0+y)/2.0 || y > x/2.0) {
        return false;
      }
      break;
    case 167:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > fourth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 168:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > sixth ||
          x > 2.0*y || y > dmin(1.0-x, 2.0*x)) {
        return false;
      }
      break;
    case 169:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > twelfth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 170:
      if (0.0 > x || x > twothr || 0.0 > y || y > half || 0.0 > z || z > 1.0 ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, x)) {
        return false;
      }
      break;
    case 171:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > sixth) return false;
      break;
    case 172:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > sixth) return false;
      break;
    case 173:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > third || y > x) {
        return false;
      }
      break;
    case 174:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > third || y > x) {
        return false;
      }
      break;
    case 175:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 176:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 177:
      if (0.0 > x || x > twothr || 0.0 > y || y > half || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, x)) {
        return false;
      }
      break;
    case 178:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > fourth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 179:
      if (0.0 > x || x > twothr || 0.0 > y || y > half || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, x)) {
        return false;
      }
      break;
    case 180:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > twelfth) return false;
      break;
    case 181:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > twelfth) return false;
      break;
    case 182:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > sixth || y > x) {
        return false;
      }
      break;
    case 183:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > sixth || y > x) {
        return false;
      }
      break;
    case 184:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > fourth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 185:
      if (0.0 > x || x > twothr || 0.0 > y || y > third || 0.0 > z || z > 1.0 ||
          x > (1.0+y)/2.0 || y > x/2.0) {
        return false;
      }
      break;
    case 186:
      if (0.0 > x || x > twothr || 0.0 > y || y > half || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, x)) {
        return false;
      }
      break;
    case 187:
      if (0.0 > x || x > twothr || 0.0 > y || y > half || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, x)) {
        return false;
      }
      break;
    case 188:
      if (0.0 > x || x > twothr || 0.0 > y || y > third || 0.0 > z || z > 1.0 ||
          x > (1.0+y)/2.0 || y > x/2.0) {
        return false;
      }
      break;
    case 189:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > half ||
          x > 2.0*y || y > dmin(1.0-x, 2.0*x)) {
        return false;
      }
      break;
    case 190:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > fourth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 191:
      if (0.0 > x || x > twothr || 0.0 > y || y > half || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, x)) {
        return false;
      }
      break;
    case 192:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > fourth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, (1.0+x)/2.0)) {
        return false;
      }
      break;
    case 193:
      if (0.0 > x || x > twothr || 0.0 > y || y > third || 0.0 > z || z > half ||
          x > (1.0+y)/2.0 || y > x/2.0) {
        return false;
      }
      break;
    case 194:
      if (0.0 > x || x > twothr || 0.0 > y || y > half || 0.0 > z || z > fourth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, x)) {
        return false;
      }
      break;
    case 195:
      if (0.0 > x || x > twothr || 0.0 > y || y > half || 0.0 > z || z > fourth ||
          x > (1.0+y)/2.0 || y > dmin(1.0-x, x)) {
        return false;
      }
      break;
    case 196:
      if (0.0 > x || x > twothr || 0.0 > y || y > twothr || 0.0 > z || z > fourth ||
          x > 2.0*y || y > dmin(1.0-x, 2.0*x)) {
        return false;
      }
      break;
    case 197:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > 1.0 || 0.0 > z || z > half || y > 1.0-x ||
          z > dmin(x, y)) {
        return false;
      }
      break;
    case 198:
      if (0.0 > x || x > half || 0.0 > y || y > half || -fourth > z || z > fourth ||
          y > x || dmax(x-half, -y) > z || z > dmin(half-x, y)) {
        return false;
      }
      break;
    case 199:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > half || 0.0 > z || z > half ||
          y > dmin(x, 1.0-x) || z > y) {
        return false;
      }
      break;
    case 200:
      if (0.0 > x || x > half || 0.0 > y || y > half || -half > z || z > half ||
          dmax(x-half, -y) > z || z > dmin(x, y)) {
        return false;
      }
      break;
    case 201:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half ||
          z > dmin(x, y)) {
        return false;
      }
      break;
    case 202:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half ||
          z > dmin(x, y)) {
        return false;
      }
      break;
    case 203:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > half || 0.0 > z || z > half ||
          y > dmin(x, 1.0-x) || z > y) {
        return false;
      }
      break;
    case 204:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth || y > x ||
          z > dmin(half-x, y)) {
        return false;
      }
      break;
    case 205:
      if (0.0 > x || x > half || 0.0 > y || y > fourth || -fourth > z || z > fourth ||
          y > dmin(x, half-x) || -y > z || z > y) {
        return false;
      }
      break;
    case 206:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || y > x ||
          z > y) {
        return false;
      }
      break;
    case 207:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half ||
          z > dmin(x, y)) {
        return false;
      }
      break;
    case 208:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth ||
          z > dmin(x, half-x, half-y)) {
        return false;
      }
      break;
    case 209:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > half || 0.0 > z || z > half ||
          y > dmin(x, 1.0-x) || z > y) {
        return false;
      }
      break;
    case 210:
      if (0.0 > x || x > half || 0.0 > y || y > half || -fourth > z || z > fourth ||
          dmax(-x, x-half, -y, y-half) > z || z > dmin(x, half-x, y, half-y)) {
        return false;
      }
      break;
    case 211:
      if (0.0 > x || x > half || 0.0 > y || y > fourth || -fourth > z || z > fourth ||
          y > dmin(x, half-x) || -y > z || z > y) {
        return false;
      }
      break;
    case 212:
      if (0.0 > x || x > half || -eighth > y || y > eighth || -eighth > z || z > eighth ||
          y > dmin(x, half-x) || -y > z || z > dmin(x, half-x)) {
        return false;
      }
      break;
    case 213:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth ||
          z > dmin(x, half-x, y, half-y)) {
        return false;
      }
      break;
    case 214:
      if (0.0 > x || x > half || 0.0 > y || y > 3.0/4 || -half > z || z > fourth ||
          dmax(-y, x-half) > z || z > dmin(half-y, 2.0*x-y, 2.0*y-x, y-2.0*x+half)) {
        return false;
      }
      break;
    case 215:
      if (-fourth > x || x > half || 0.0 > y || y > 3.0/4 || 0.0 > z || z > half ||
          x > y || y > x+half || (y-x)/2.0 > z ||
          z > dmin(y, (-4.0*x-2.0*y+3.0)/2.0, (3.0-2.0*x-2.0*y)/4)) {
        return false;
      }
      break;
    case 216:
      if (-threig > x || x > eighth || -eighth > y || y > eighth || -eighth > z ||
          z > threig || dmax(x, y, y-x-eighth) > z || z > y+fourth) {
        return false;
      }
      break;
    case 217:
      if (0.0 > x || x > 1.0 || 0.0 > y || y > half || 0.0 > z || z > half ||
          y > dmin(x, 1.0-x) || z > y) {
        return false;
      }
      break;
    case 218:
      if (0.0 > x || x > half || 0.0 > y || y > fourth || -fourth > z || z > fourth ||
          y > dmin(x, half-x) || -y > z || z > y) {
        return false;
      }
      break;
    case 219:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || y > x ||
          z > y) {
        return false;
      }
      break;
    case 220:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half ||
          z > dmin(x, y)) {
        return false;
      }
      break;
    case 221:
      if (0.0 > x || x > half || 0.0 > y || y > fourth || -fourth > z || z > fourth ||
          y > dmin(x, half-x) || -y > z || z > y) {
        return false;
      }
      break;
    case 222:
      if (fourth > x || x > half || fourth > y || y > half || 0.0 > z || z > half ||
          z > dmin(x, y)) {
        return false;
      }
      break;
    case 223:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || y > x ||
          z > y) {
        return false;
      }
      break;
    case 224:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > half || y > x ||
          z > y) {
        return false;
      }
      break;
    case 225:
      if (0.0 > x || x > half || 0.0 > y || y > half || 0.0 > z || z > fourth ||
          z > dmin(x, half-x, y, half-y)) {
        return false;
      }
      break;
    case 226:
      if (0.0 > x || x > half || 0.0 > y || y > half || -fourth > z || z > fourth ||
          y > x || dmax(x-half, -y) > z || z > dmin(half-x, y)) {
        return false;
      }
      break;
    case 227:
      if (0.0 > z || z > half || 0.0 > y || y > half || 0.0 > x || x > 1.0) return false;
    default:
      break;
  }
  
  return result;
}
