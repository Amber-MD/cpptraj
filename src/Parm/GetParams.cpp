#include "GetParams.h"
#include "DihedralParmSet.h"
#include "ParameterSet.h"
#include "../Atom.h"
#include "../AtomType.h"
#include "../CmapParmHolder.h"
#include "../CpptrajStdio.h"
#include "../ParameterTypes.h"
#include "../Topology.h"
#include "../TypeNameHolder.h"

using namespace Cpptraj::Parm;

/** CONSTRUCTOR */
GetParams::GetParams() :
  debug_(0)
{}

/** Set debug level */
void GetParams::SetDebug(int debugIn) {
  debug_ = debugIn;
}

static void paramOverwriteWarning(const char* type) {
  mprintf("Warning: An existing %s parameter would have been overwritten. This\n"
          "Warning:  usually means that the atom type information in the Topology is\n"
          "Warning:  incomplete. This can happen for example with Chamber topologies\n"
          "Warning:  if the original atom type names were > 4 characters.\n", type);
  mprintf("Warning: The %s parameters in this topology may now be incorrect.\n", type);
}

// GetBondParams()
void GetParams::GetBondParams(ParmHolder<BondParmType>& BP, std::vector<Atom> const& atoms, BondArray const& bonds, BondParmArray const& bpa) {
  for (BondArray::const_iterator b = bonds.begin(); b != bonds.end(); ++b)
  {
    if (b->Idx() != -1) {
      TypeNameHolder types(2);
      types.AddName( atoms[b->A1()].Type() );
      types.AddName( atoms[b->A2()].Type() );
      Cpptraj::Parm::RetType ret = BP.AddParm( types, bpa[b->Idx()], false );
      if (ret == Cpptraj::Parm::ERR)
        paramOverwriteWarning("bond");
    }
  }
}

// GetAngleParams()
void GetParams::GetAngleParams(ParmHolder<AngleParmType>& AP, std::vector<Atom> const& atoms, AngleArray const& angles, AngleParmArray const& apa) {
  for (AngleArray::const_iterator b = angles.begin(); b != angles.end(); ++b)
  {
    if (b->Idx() != -1) {
      TypeNameHolder types(3);
      types.AddName( atoms[b->A1()].Type() );
      types.AddName( atoms[b->A2()].Type() );
      types.AddName( atoms[b->A3()].Type() );
      Cpptraj::Parm::RetType ret = AP.AddParm( types, apa[b->Idx()], false );
      if (ret == Cpptraj::Parm::ERR)
        paramOverwriteWarning("angle");
    }
  }
}

// GetImproperParams()
void GetParams::GetImproperParams(ImproperParmHolder& IP, std::vector<Atom> const& atoms, DihedralArray const& imp, DihedralParmArray const& ipa) {
  IP.SetRequireExactMatch(true);
  DihedralParmSet impPrm;
  for (DihedralArray::const_iterator b = imp.begin(); b != imp.end(); ++b)
  {
    if (b->Idx() != -1) {
      TypeNameHolder types(4);
      types.AddName( atoms[b->A1()].Type() );
      types.AddName( atoms[b->A2()].Type() );
      types.AddName( atoms[b->A3()].Type() );
      types.AddName( atoms[b->A4()].Type() );
      Cpptraj::Parm::RetType ret = impPrm.AddDihParm( types, ipa[b->Idx()], false );
      if (ret == Cpptraj::Parm::ERR)
        paramOverwriteWarning("improper");
    }
  }
  if (impPrm.ToImpParm(IP)) {
    mprinterr("Internal Error: GetParams::GetImproperParams(): Could not transfer parameters from set to holder.\n");
  }
}

// GetDihedralParams()
void GetParams::GetDihedralParams(DihedralParmHolder& DP, ImproperParmHolder& IP, std::vector<Atom> const& atoms, DihedralArray const& dih, DihedralParmArray const& dpa) {
  IP.SetRequireExactMatch(true);
  DihedralParmSet dihPrm;
  DihedralParmSet impPrm;
  for (DihedralArray::const_iterator b = dih.begin(); b != dih.end(); ++b)
  {
    if (b->Idx() != -1) {
      TypeNameHolder types(4);
      types.AddName( atoms[b->A1()].Type() );
      types.AddName( atoms[b->A2()].Type() );
      types.AddName( atoms[b->A3()].Type() );
      types.AddName( atoms[b->A4()].Type() );
      //mprintf("DEBUG: dihedral %li ( %i %i %i %i )\n", b - dih.begin() + 1, b->A1()+1, b->A2()+1, b->A3()+1, b->A4()+1);
      //mprintf("DEBUG: dihedral %li %s %s %s %s idx=%i type=%i PK=%g PN=%g Phase=%g SCEE=%g SCNB=%g\n", b - dih.begin() + 1,
      //        *(types[0]), *(types[1]), *(types[2]), *(types[3]), b->Idx(), (int)b->Type(),
      //        dpa[b->Idx()].Pk(), dpa[b->Idx()].Pn(), dpa[b->Idx()].Phase(), dpa[b->Idx()].SCEE(), dpa[b->Idx()].SCNB());
      Cpptraj::Parm::RetType ret;
      if (b->IsImproper()) {
        ret = impPrm.AddDihParm( types, dpa[b->Idx()], false );
      } else {
        ret = dihPrm.AddDihParm( types, dpa[b->Idx()], false );
      }
      // DEBUG
      //if (ret == Cpptraj::Parm::ADDED) {
      //  mprintf("DEBUG: Added %s %s %s %s idx=%i isImproper=%i\n", *(types[0]), *(types[1]), *(types[2]), *(types[3]), b->Idx(), (int)b->IsImproper());
      //}
      if (ret == Cpptraj::Parm::ERR) {
        paramOverwriteWarning("dihedral");
        mprintf("Warning: Dihedral %s %s %s %s PK=%g PN=%g Phase=%g SCEE=%g SCNB=%g\n",
                *(types[0]), *(types[1]), *(types[2]), *(types[3]),
                dpa[b->Idx()].Pk(), dpa[b->Idx()].Pn(), dpa[b->Idx()].Phase(), dpa[b->Idx()].SCEE(), dpa[b->Idx()].SCNB());
        //bool found;
        //DihedralParmArray dpa = DP.FindParam(types, found);
        //mprintf("Warning: Existing params:\n");
        //for (DihedralParmArray::const_iterator d = dpa.begin(); d != dpa.end(); ++d)
        //  mprintf("Warning:\t\tPK=%g PN=%g Phase=%g SCEE=%g SCNB=%g\n",
        //          d->Pk(), d->Pn(), d->Phase(), d->SCEE(), d->SCNB());
      }
    }
  }
  if (impPrm.ToImpParm(IP)) {
    mprinterr("Internal Error: GetParams::GetImproperParams(): Could not transfer parameters from set to holder.\n");
  }
  if (dihPrm.ToDihParm(DP)) {
    mprinterr("Internal Error: GetParams::GetDihedralParams(): Could not transfer parameters from set to holder.\n");
  }
}

/// \return an error if cmap atom names do not match
static inline int check_cmap_atom_name(NameType const& n0, NameType const& n1)
{
  if (n0 != n1) {
    mprinterr("Error: CMAP term atom name %s does not match expected CMAP term atom name %s\n",
              *n1, *n0);
    return 1;
  }
  return 0;
}

/** Get existing CMAP parameters.
  * Unlike other parameters, CMAPs are big and are uniquely identified by
  * a combination of residue and 5 atom names. Depending on if they were
  * read in from a parameter file or a topology, they may or may not
  * have residue/atom name information, which is needed for assignment.
  * This info needs to be generated if it is missing.
  */
int GetParams::GetCmapParams(CmapParmHolder& cmapParm, CmapArray const& cmapTerms,
                             CmapGridArray const& cmapGrids,
                             std::vector<Atom> const& atoms, std::vector<Residue> const& residues)
const
{
  if (cmapGrids.empty() || cmapTerms.empty()) {
    //mprintf("DEBUG: CMAP grids/terms are empty. No parameters to get.\n");
    return 0;
  }
  // Check if we need to generate residue/atom information for grids.
  bool needResAtomInfo = false;
  for (CmapGridArray::const_iterator it = cmapGrids.begin(); it != cmapGrids.end(); ++it)
  {
    if (it->ResNames().empty() || it->AtomNames().empty()) {
      needResAtomInfo = true;
      break;
    }
  }
  // Mark off which terms need to be added.
  std::vector<bool> addGrid( cmapGrids.size(), false );
  for (CmapArray::const_iterator cm = cmapTerms.begin(); cm != cmapTerms.end(); ++cm) {
    if (cm->Idx() != -1)
      addGrid[cm->Idx()] = true;
  }
  // Add grids, adding res/atom info if needed
  if (needResAtomInfo) {
    mprintf("CMAP terms need residue/atom info.\n");
    for (unsigned int idx = 0; idx != cmapGrids.size(); idx++) {
      if (addGrid[idx]) {
        // Figure out residue names this cmap applies to.
        std::set<NameType> resNames;
        std::vector<NameType> atomNames;
        atomNames.reserve(5);
        for (CmapArray::const_iterator cm = cmapTerms.begin(); cm != cmapTerms.end(); ++cm) {
          if (cm->Idx() == (int)idx) {
            int resnum = atoms[cm->A2()].ResNum();
            resNames.insert( residues[resnum].Name() );
            if (atomNames.empty()) {
              atomNames.push_back( atoms[cm->A1()].Name() );
              atomNames.push_back( atoms[cm->A2()].Name() );
              atomNames.push_back( atoms[cm->A3()].Name() );
              atomNames.push_back( atoms[cm->A4()].Name() );
              atomNames.push_back( atoms[cm->A5()].Name() );
            } else {
              // Check atom names
              if (check_cmap_atom_name(atomNames[0], atoms[cm->A1()].Name())) return 1;
              if (check_cmap_atom_name(atomNames[1], atoms[cm->A2()].Name())) return 1;
              if (check_cmap_atom_name(atomNames[2], atoms[cm->A3()].Name())) return 1;
              if (check_cmap_atom_name(atomNames[3], atoms[cm->A4()].Name())) return 1;
              if (check_cmap_atom_name(atomNames[4], atoms[cm->A5()].Name())) return 1;
            }
          }
        }
        //mprintf("DEBUG: Cmap term %u residues", idx);
        //for (std::set<NameType>::const_iterator it = resNames.begin(); it != resNames.end(); ++it)
        //  mprintf(" %s", *(*it));
        //mprintf(" Atoms={");
        //for (std::vector<NameType>::const_iterator it = atomNames.begin(); it != atomNames.end(); ++it)
        //  mprintf(" %s", *(*it));
        //mprintf(" }\n");
        CmapGridType newGrid = cmapGrids[idx];
        // Add the atom/res info to the grid
        newGrid.SetNumCmapRes( resNames.size() );
        for (std::set<NameType>::const_iterator it = resNames.begin(); it != resNames.end(); ++it)
          newGrid.AddResName( it->Truncated() ); // FIXME just use NameType
        for (std::vector<NameType>::const_iterator it = atomNames.begin(); it != atomNames.end(); ++it)
          newGrid.AddAtomName( it->Truncated() ); // FIXME just use NameType
        if (!newGrid.CmapIsValid()) {
          mprinterr("Error: CMAP is not valid, could not get parameter.\n");
          mprinterr("Error:   Term %u residues", idx);
          for (std::set<NameType>::const_iterator it = resNames.begin(); it != resNames.end(); ++it)
            mprinterr(" %s", *(*it));
          mprinterr(" Atoms={");
          for (std::vector<NameType>::const_iterator it = atomNames.begin(); it != atomNames.end(); ++it)
            mprinterr(" %s", *(*it));
          mprinterr(" }\n");
        }
        // Set a default title if needed
        if (newGrid.Title().empty())
          newGrid.SetTitle( "CMAP for " + newGrid.ResNames().front() );
        Cpptraj::Parm::RetType ret = cmapParm.AddParm( newGrid, false, debug_ );
        if (ret == Cpptraj::Parm::ERR)
          paramOverwriteWarning("CMAP");
      } // END if adding existing grid to parms
    } // END loop over existing grids
  } else {
    mprintf("CMAP terms have residue/atom info.\n");
    for (unsigned int idx = 0; idx != cmapGrids.size(); idx++) {
      if (addGrid[idx]) {
        Cpptraj::Parm::RetType ret = cmapParm.AddParm( cmapGrids[idx], false, debug_ );
        if (ret == Cpptraj::Parm::ERR)
          paramOverwriteWarning("CMAP");
      }
    }
  }
  return 0;
}

/** \param atomTypesOut Output array of atom types and indivudual LJ parameters.
  * \param LJ612out Output array of LJ 6-12 pair parameters.
  * \param LJ14out Output array of LJ 6-12 1-4 pair parameters.
  * \param LJ1012out Output array of LJ 10-12 pair parameters.
  * \param atoms Current array of atoms.
  * \param NB0 Current nonbond parameters.
  */
void GetParams::GetLJAtomTypes(ParmHolder<AtomType>& atomTypesOut,
                               ParmHolder<NonbondType>& LJ612out,
                               ParmHolder<NonbondType>& LJ14out,
                               ParmHolder<double>& LJCout,
                               ParmHolder<HB_ParmType>& LJ1012out,
                               std::vector<Atom> const& atoms,
                               NonbondParmType const& NB0)
const
{
  if (NB0.HasNonbond()) {
    if (debug_ > 0) mprintf("DEBUG: Topology has nonbond parameters.\n");
    bool hasLJ14 = !NB0.LJ14().empty();
    if (debug_ > 0 && hasLJ14)
      mprintf("DEBUG: Topology has 1-4 nonbond parameters.\n");
    bool hasLJC = !NB0.LJC_Array().empty();
    if (debug_ > 0 && hasLJC)
      mprintf("DEBUG: Topology has LJC nonbond paramters.\n");
    // Nonbonded parameters are present.
    for (std::vector<Atom>::const_iterator atm = atoms.begin(); atm != atoms.end(); ++atm)
    {
      if (!atm->HasType()) {
        mprintf("Warning: Topology has nonbond parameters but atom %s has no type.\n", *(atm->Name()));
        continue;
      }
      TypeNameHolder atype( atm->Type() );
      // Check for self parameters to back-calculate LJ depth/radius
      int idx = NB0.GetLJindex( atm->TypeIndex(), atm->TypeIndex() );
      AtomType thisType;
      if (idx > -1) {
        // Has LJ 6-12 parameters
        NonbondType const& LJ = NB0.NBarray( idx );
        thisType = AtomType(LJ.Radius(), LJ.Depth(), atm->Mass(), atm->Polar());
        if (hasLJ14) {
          NonbondType const& lj14 = NB0.LJ14( idx );
          thisType.SetLJ14( LJparmType(lj14.Radius(), lj14.Depth()) );
        }
      } else {
        // Has LJ 10-12 parameters
        thisType = AtomType(atm->Mass(), atm->Polar());
      }
      thisType.SetTypeIdx( atm->TypeIndex() );
      Cpptraj::Parm::RetType ret = atomTypesOut.AddParm( atype, thisType, true );
      if (debug_ > 0 && ret == Cpptraj::Parm::ADDED) {
        mprintf("DEBUG: New atom type: %s R=%g D=%g M=%g P=%g\n", *(atype[0]),
                thisType.LJ().Radius(), thisType.LJ().Depth(), thisType.Mass(), thisType.Polarizability());
        if (hasLJ14)
          mprintf("DEBUG: New LJ14 atom type: %s R=%g D=%g\n", *(atype[0]), thisType.LJ14().Radius(), thisType.LJ14().Depth());
      }
    }
    // Do atom type pairs, check for off-diagonal elements.
    // Explicitly store pairs instead of regenerating to avoid round-off issues.
    //GetLJterms(atomTypesOut, LJ612out, &LJ1012out, NB0, false, debug_);
    //if (hasLJ14)
    //  GetLJterms(LJ14typesOut, LJ14out, 0, NB0, true, debug_);
    unsigned int nModifiedOffDiagonal = 0;
    unsigned int nModified14OffDiagonal = 0;
    for (ParmHolder<AtomType>::const_iterator i1 = atomTypesOut.begin(); i1 != atomTypesOut.end(); ++i1)
    {
      for (ParmHolder<AtomType>::const_iterator i2 = i1; i2 != atomTypesOut.end(); ++i2)
      {
        NameType const& name1 = i1->first[0];
        NameType const& name2 = i2->first[0];
        TypeNameHolder types(2);
        types.AddName( name1 );
        types.AddName( name2 );
        // Extract original nonbonded parameters for this type pair.
        AtomType const& type1 = i1->second;
        AtomType const& type2 = i2->second;
        int idx1 = type1.OriginalIdx();
        int idx2 = type2.OriginalIdx();
        int idx = NB0.GetLJindex( idx1, idx2 );
        if (idx < 0) {
          // This is LJ 10-12.
          //mprinterr("Error: No off-diagonal LJ for  %s %s (%i %i)\n",
          //          *name1, *name2, idx1, idx2);
          //return;
          if (debug_ > 0)
            mprintf("DEBUG: LJ 10-12 parameters detected for %s %s (%i %i)\n",
                    *name1, *name2, idx1, idx2);
          LJ1012out.AddParm( types, NB0.HBarray((-idx)-1), false );
        } else {
          // This is LJ 6-12.
          // Determine what A and B parameters would be.
          NonbondType lj0 = type1.LJ().Combine_LB( type2.LJ() );
        
          NonbondType lj1 = NB0.NBarray( idx );
          // Compare them
          if (lj0 != lj1) {
            nModifiedOffDiagonal++;
            if (debug_ > 0) {
              double deltaA = fabs(lj0.A() - lj1.A());
              double deltaB = fabs(lj0.B() - lj1.B());
              mprintf("DEBUG: Potential off-diagonal LJ: %s %s expect A=%g B=%g, actual A=%g B=%g\n",
                      *name1, *name2, lj0.A(), lj0.B(), lj1.A(), lj1.B());
              mprintf("DEBUG:\tdeltaA= %g    deltaB= %g\n", deltaA, deltaB);
              double pe_a = (fabs(lj0.A() - lj1.A()) / lj0.A());
              double pe_b = (fabs(lj0.B() - lj1.B()) / lj0.B());
              mprintf("DEBUG:\tPEA= %g  PEB= %g\n", pe_a, pe_b);
            }
          }
          LJ612out.AddParm( types, lj1, false );
          if (hasLJ14) {
            // This is LJ 6-12 1-4.
            // Determine what A and B parameters would be.
            lj0 = type1.LJ14().Combine_LB( type2.LJ14() );
        
            lj1 = NB0.LJ14( idx );
            // Compare them
            if (lj0 != lj1) {
              nModified14OffDiagonal++;
              if (debug_ > 0) {
                double deltaA = fabs(lj0.A() - lj1.A());
                double deltaB = fabs(lj0.B() - lj1.B());
                mprintf("DEBUG: Potential off-diagonal LJ 1-4: %s %s expect A=%g B=%g, actual A=%g B=%g\n",
                        *name1, *name2, lj0.A(), lj0.B(), lj1.A(), lj1.B());
                mprintf("DEBUG:\tdeltaA= %g    deltaB= %g\n", deltaA, deltaB);
                double pe_a = (fabs(lj0.A() - lj1.A()) / lj0.A());
                double pe_b = (fabs(lj0.B() - lj1.B()) / lj0.B());
                mprintf("DEBUG:\tPEA= %g  PEB= %g\n", pe_a, pe_b);
              }
            }
            LJ14out.AddParm( types, lj1, false );
          } // END hasLJ14
          if (hasLJC) {
            // This is LJC
            double ljc1 = NB0.LJC_Array( idx );
            LJCout.AddParm( types, ljc1, false );
          }
        }
      } // END inner loop over atom types
    } // END outer loop over atom types
    if (nModifiedOffDiagonal > 0)
      mprintf("Warning: %u modified off-diagonal LJ terms present.\n", nModifiedOffDiagonal);
    if (nModified14OffDiagonal > 0)
      mprintf("Warning: %u modified off-diagonal LJ 1-4 terms present.\n", nModified14OffDiagonal);
  } else {
    //if (!atoms.empty()) mprintf("DEBUG: Topology does not have nonbond parameters.\n");
    // No nonbonded parameters. Just save mass/polarizability.
    for (std::vector<Atom>::const_iterator atm = atoms.begin(); atm != atoms.end(); ++atm)
      if (atm->HasType() > 0)
        atomTypesOut.AddParm( TypeNameHolder(atm->Type()), AtomType(atm->Mass(), atm->Polar()), true );
  }
}

/** \return ParameterSet for this Topology. */
ParameterSet GetParams::GetParameters(Topology const& topIn) const {
  ParameterSet Params;
  // Atom LJ types and other nonbonded parameters.
  GetLJAtomTypes( Params.AT(), Params.NB(), Params.NB14(), Params.LJC(), Params.HB(),
                  topIn.Atoms(), topIn.Nonbond() );
  // Bond parameters.
  GetBondParams( Params.BP(), topIn.Atoms(), topIn.Bonds(), topIn.BondParm() );
  GetBondParams( Params.BP(), topIn.Atoms(), topIn.BondsH(), topIn.BondParm() );
  // Angle parameters.
  GetAngleParams( Params.AP(), topIn.Atoms(), topIn.Angles(), topIn.AngleParm() );
  GetAngleParams( Params.AP(), topIn.Atoms(), topIn.AnglesH(), topIn.AngleParm() );
  // Dihedral parameters.
  GetDihedralParams( Params.DP(), Params.IP(), topIn.Atoms(), topIn.Dihedrals(), topIn.DihedralParm() );
  GetDihedralParams( Params.DP(), Params.IP(), topIn.Atoms(), topIn.DihedralsH(), topIn.DihedralParm() );
  // CHARMM parameters
  if (!topIn.UB().empty()) {
    // UB parameters
    GetBondParams(Params.UB(), topIn.Atoms(), topIn.UB(), topIn.UBparm() );
  }
  if (!topIn.Impropers().empty()) {
    // Impropers
    GetImproperParams( Params.IP(), topIn.Atoms(), topIn.Impropers(), topIn.ImproperParm() );
  }
  // CMAPs
  if (GetCmapParams( Params.CMAP(), topIn.Cmap(), topIn.CmapGrid(), topIn.Atoms(), topIn.Residues() )) {
    mprinterr("Error: Could not get CMAP parameters.\n"); // TODO fatal?
  }

  return Params;
}
/** \return Total number of unique atom types. */
unsigned int GetParams::NuniqueAtomTypes(Topology const& topIn) const {
  ParmHolder<int> currentAtomTypes;
  for (Topology::atom_iterator atm = topIn.begin(); atm != topIn.end(); ++atm)
  {
    if (atm->HasType()) {
      TypeNameHolder atype( atm->Type() );
      // Find in currentAtomTypes.
      bool found;
      currentAtomTypes.FindParam( atype, found );
      if (!found) {
        currentAtomTypes.AddParm( atype, atm->TypeIndex(), false );
      }
    }
  }
  if (debug_ > 0) {
    mprintf("DEBUG: Unique atom types in %s\n", topIn.c_str());
    for (ParmHolder<int>::const_iterator it = currentAtomTypes.begin();
                                         it != currentAtomTypes.end(); ++it)
      mprintf("\t\t%s %i\n", *(it->first[0]), it->second);
  }
  return currentAtomTypes.size();
}
