#include "Exec_Build.h"
#include "AssociatedData_Connect.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h" // For casting DataSet_Parameters to ParameterSet
#include "ParmFile.h"
#include "StringRoutines.h" // integerToString
#include "Trajout_Single.h"
#include "Parm/AssignParams.h"
#include "Parm/LJ1264_Params.h"
#include "Structure/AddIons.h"
#include "Structure/Builder.h"
#include "Structure/Creator.h"
#include "Structure/Disulfide.h"
#include "Structure/HisProt.h"
#include "Structure/PdbCleaner.h"
#include "Structure/ResStatArray.h"
#include "Structure/Solvate.h"
#include "Structure/SugarBuilder.h"
#include "Structure/Sugar.h"
#include "StructureCheck.h"
#include <set> // For warning about missing residue templates
#include <cmath> // fabs

/** CONSTRUCTOR */
Exec_Build::Exec_Build() :
  Exec(GENERAL),
  debug_(0),
  check_box_natom_(5000), // TODO make user specifiable
  check_structure_(true),
  keepMissingSourceAtoms_(false),
  requireAllInputAtoms_(false),
  outCrdPtr_(0)
{}

/** Search in array of atom bonding pairs for given bonding pair. */
bool Exec_Build::hasBondingPair(IParray const& bpairs, Ipair const& bpair) {
  for (IParray::const_iterator it = bpairs.begin(); it != bpairs.end(); ++it)
    if (*it == bpair) return true;
  return false;
}

/** \return True if target residue is in array of residue connections. */
bool Exec_Build::resIsConnected(Iarray const& resConnections, int tgtRes) {
  for (Iarray::const_iterator it = resConnections.begin(); it != resConnections.end(); ++it)
    if (*it == tgtRes) return true;
  return false;
}

/** Use given templates to construct a final molecule. */
int Exec_Build::FillAtomsWithTemplates(Topology& topOut, Frame& frameOut,
                                       Topology const& topIn, Frame const& frameIn,
                                       Cpptraj::Structure::Creator const& creator,
                                       std::vector<BondType> const& topInBonds)
{
  // Array of head/tail connect atoms for each residue
  Iarray resHeadAtoms;
  Iarray resTailAtoms;
  std::vector<Cpptraj::Structure::TerminalType> ResTermTypes;
  resHeadAtoms.reserve( topIn.Nres() );
  resTailAtoms.reserve( topIn.Nres() );
  ResTermTypes.reserve( topIn.Nres() );
  // Array of templates for each residue
  std::vector<DataSet_Coords*> ResTemplates;
  ResTemplates.reserve( topIn.Nres() );

  t_fill_template_.Start();
  // Initial loop to try to match residues to templates
  int newNatom = 0;
  unsigned int n_no_template_found = 0;
  std::set<NameType> missing_templates;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    Residue const& currentRes = topIn.Res(ires);
    if (debug_ > 0)
      mprintf("DEBUG: ---------- Processing Residue %s ---------- \n", topIn.TruncResNameNum(ires).c_str());
    int pres = ires - 1;
    int nres = ires + 1;
    // Determine if this is a terminal residue
    Cpptraj::Structure::TerminalType resTermType;
    if (currentRes.IsTerminal()) {
      resTermType = Cpptraj::Structure::END_TERMINAL;
      if (debug_ > 0)
        mprintf("DEBUG: %s End terminal due to TERMINAL status.\n", topIn.TruncResNameNum(ires).c_str());
    } else if (ires == 0 && topIn.Nres() > 1) {
      resTermType = Cpptraj::Structure::BEG_TERMINAL;
      if (debug_ > 0)
        mprintf("DEBUG: %s Begin terminal due to first residue.\n", topIn.TruncResNameNum(ires).c_str());
    } else if (pres > -1 && topIn.Res(pres).IsTerminal()) {
      resTermType = Cpptraj::Structure::BEG_TERMINAL;
      if (debug_ > 0)
        mprintf("DEBUG: %s Begin terminal due to previous residue TERMINAL status.\n", topIn.TruncResNameNum(ires).c_str());
    } else if (nres < topIn.Nres() && (topIn.Res(nres).ChainID() != currentRes.ChainID())// ||
//                                        (topIn.Res(nres).OriginalResNum() == currentRes.OriginalResNum() &&
//                                         topIn.Res(nres).Icode()          != currentRes.Icode()
//                                        )
//                                      )
              )
    {
      resTermType = Cpptraj::Structure::END_TERMINAL;
      if (debug_ > 0)
        mprintf("DEBUG: %s End terminal due to chain ID.\n", topIn.TruncResNameNum(ires).c_str());
    } else if (nres == topIn.Nres()) {
      resTermType = Cpptraj::Structure::END_TERMINAL;
      if (debug_ > 0)
        mprintf("DEBUG: %s End terminal due to last residue.\n", topIn.TruncResNameNum(ires).c_str());
    } else {
      resTermType = Cpptraj::Structure::NON_TERMINAL;
    }
    if (debug_ > 0)
      mprintf("DEBUG: Residue type: %s terminal (IsTerminal=%i)\n", Cpptraj::Structure::terminalStr(resTermType), (int)currentRes.IsTerminal());
    ResTermTypes.push_back( resTermType );
    // Identify a template based on the residue name.
    DataSet_Coords* resTemplate = creator.IdTemplateFromResname(currentRes.Name(), resTermType);
    if (resTemplate == 0) {
      // Residue has no template.
      n_no_template_found++;
      mprintf("Warning: No template found for residue %s\n", topIn.TruncResNameOnumId(ires).c_str());
      missing_templates.insert( currentRes.Name() );
      newNatom += currentRes.NumAtoms();
      // Head and tail atoms are blank
      resHeadAtoms.push_back( -1 );
      resTailAtoms.push_back( -1 );
    } else {
      // Residue has a template.
      if (debug_ > 0)
        mprintf("\tTemplate %s being used for residue %s\n",
                resTemplate->legend(), topIn.TruncResNameOnumId(ires).c_str());
      int nTgtAtomsMissing = 0;
      if (keepMissingSourceAtoms_)
        nTgtAtomsMissing = creator.CountAtomsMissingFromTemplate( topIn, ires, resTemplate );
      // Save the head and tail atoms
      AssociatedData* ad = resTemplate->GetAssociatedData(AssociatedData::CONNECT);
      if (ad == 0) {
        if (debug_ > 0)
          mprintf("DEBUG: Unit '%s' does not have CONNECT data.\n", resTemplate->legend());
        resHeadAtoms.push_back( -1 );
        resTailAtoms.push_back( -1 );
      } else {
        AssociatedData_Connect const& CONN = static_cast<AssociatedData_Connect const&>( *ad );
        if (CONN.NconnectAtoms() < 2) {
          mprinterr("Error: Not enough connect atoms in unit '%s'\n", resTemplate->legend());
          return 1;
        }
        if (CONN.Connect()[0] > -1)
          resHeadAtoms.push_back( CONN.Connect()[0] + newNatom );
        else
          resHeadAtoms.push_back( -1 );
        if (CONN.Connect()[1] > -1)
          resTailAtoms.push_back( CONN.Connect()[1] + newNatom );
        else
          resTailAtoms.push_back( -1 );
      }
      // Update # of atoms
      newNatom += resTemplate->Top().Natom() + nTgtAtomsMissing;
    }
    ResTemplates.push_back( resTemplate );
  }
  mprintf("\tFinal structure should have %i atoms.\n", newNatom);
  if (n_no_template_found > 0) {
    mprintf("Warning: No template was found for %u residues.\n", n_no_template_found);
    mprintf("Warning: Residue names:");
    for (std::set<NameType>::const_iterator rit = missing_templates.begin();
                                            rit != missing_templates.end(); ++rit)
      mprintf(" %s", *(*rit));
    mprintf("\n");
    mprintf("Warning: This may indicate that you have not loaded a library file\n"
            "Warning:   or have not loaded the correct force field.\n");
  }
  frameOut.SetupFrame( newNatom );
  // Clear frame so that AddXYZ can be used
  frameOut.ClearAtoms();

  // -----------------------------------
  // hasPosition - for each atom in topOut, status on whether atom in frameOut needs building
  Cpptraj::Structure::Builder::Barray hasPosition;
  hasPosition.reserve( newNatom );

  // Hold atom offsets needed when building residues
  Iarray AtomOffsets;
  AtomOffsets.reserve( topIn.Nres() );

  // For existing inter-residue bonding, use residue # and atom name since
  // atom numbering may change if atoms are added from templates.
  // TODO make this a class var so disulfide/sugar prep can use it.
  typedef std::pair<int,NameType> ResAtPair;
  typedef std::vector<ResAtPair> ResAtArray;
  ResAtArray detectedInterResBonds;

  // Hold template atom names corressponding to source atom names.
  typedef std::vector<NameType> NameArray;
  NameArray SourceAtomNames;
  SourceAtomNames.resize( topIn.Natom() );

  // Loop for setting up atoms in the topology from residues or residue templates.
  int nAtomsNotInTemplates = 0;
  int nRefAtomsMissing = 0;
  int nAtomsMissingTypes = 0;
  bool has_bfac = !topIn.Bfactor().empty();
  bool has_occ  = !topIn.Occupancy().empty();
  bool has_pdbn = !topIn.PdbSerialNum().empty();
  if (debug_ > 0)
    mprintf("DEBUG: Input top has_bfac=%i  has_occ=%i  has_pdbn=%i\n",
            (int)has_bfac, (int)has_occ, (int)has_pdbn);
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    if (debug_ > 0)
      mprintf("\tAdding atoms for residue %s\n", topIn.TruncResNameOnumId(ires).c_str());
    int atomOffset = topOut.Natom();
    //mprintf("DEBUG: atom offset is %i\n", atomOffset);
    DataSet_Coords* resTemplate = ResTemplates[ires];
    IParray intraResBonds;
    if (resTemplate == 0) {
      // ----- No template. Just add the atoms. ------------
      Residue const& currentRes = topIn.Res(ires);
      AtomOffsets.push_back( -1 );
      for (int itgt = currentRes.FirstAtom(); itgt != currentRes.LastAtom(); ++itgt)
      {
        // Track intra-residue bonds
        Atom sourceAtom = topIn[itgt];
        if (!sourceAtom.HasType())
          nAtomsMissingTypes++;
        SourceAtomNames[itgt] = sourceAtom.Name();
        int at0 = itgt - currentRes.FirstAtom() + atomOffset;
        for (Atom::bond_iterator bat = sourceAtom.bondbegin(); bat != sourceAtom.bondend(); ++bat) {
          if ( topIn[*bat].ResNum() == ires ) {
            // Intra-residue
            int at1 = *bat - currentRes.FirstAtom() + atomOffset;
            if (at1 > at0) {
              //mprintf("Will add bond between %i and %i (original %i and %i)\n", at0+1, at1+1, itgt+1, *bat + 1);
              intraResBonds.push_back( Ipair(at0, at1) );
            }
          } else {
            // Inter-residue. Only record if bonding to a previous residue.
            if (topIn[*bat].ResNum() < ires) {
              detectedInterResBonds.push_back( ResAtPair(ires, sourceAtom.Name()) );
              detectedInterResBonds.push_back( ResAtPair(topIn[*bat].ResNum(), topIn[*bat].Name()) );
            }
          }
        }
        sourceAtom.ClearBonds(); // FIXME AddTopAtom should clear bonds
        topOut.AddTopAtom( sourceAtom, currentRes );
        frameOut.AddVec3( Vec3(frameIn.XYZ(itgt)) );
        hasPosition.push_back( true );
        if (has_bfac)
          topOut.AddBfactor( topIn.Bfactor()[itgt] );
        if (has_occ)
          topOut.AddOccupancy( topIn.Occupancy()[itgt] );
        if (has_pdbn)
          topOut.AddPdbSerialNum( topIn.PdbSerialNum()[itgt] );
      }
    } else {
      // ----- A template exists for this residue. ---------
      Residue currentRes = topIn.Res(ires);
      // Use template residue name.
      // To match LEaP behavior, if the template name is > 3 characters,
      // truncate to the last 3 characters.
      NameType const& templateResName = resTemplate->Top().Res(0).Name();
      if (templateResName.len() < 4)
        currentRes.SetName( templateResName );
      else
        currentRes.SetName( NameType( (*templateResName) + ( templateResName.len() - 3) ) );
      // Map source atoms to template atoms.
      int nTgtAtomsMissing = 0;
      std::vector<int> map = creator.MapAtomsToTemplate( topIn, ires, resTemplate, SourceAtomNames, nTgtAtomsMissing );
      if (debug_ > 1) {
        mprintf("\t  Atom map:\n");
        // DEBUG - print map
        for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
          mprintf("\t\t%6i %6s =>", iref+1, *(resTemplate->Top()[iref].Name()));
          if (map[iref] == -1)
            mprintf(" No match\n");
          else
            mprintf(" %6i %6s\n", map[iref]+1, *(topIn[map[iref]].Name()));
        }
      }
      // Map template atoms back to source atoms.
      std::vector<int> pdb(currentRes.NumAtoms(), -1);
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        if (map[iref] != -1)
          pdb[map[iref]-currentRes.FirstAtom()] = iref;
      }
      if (debug_ > 1) {
        mprintf("\t  PDB atom map:\n");
        for (int itgt = 0; itgt != currentRes.NumAtoms(); itgt++) {
          mprintf("\t\t%6i %6s =>", itgt+1, *(topIn[itgt+currentRes.FirstAtom()].Name()));
          if (pdb[itgt] == -1)
            mprintf(" Not in template\n");
          else
            mprintf(" %6i %6s\n", pdb[itgt]+1, *(resTemplate->Top()[pdb[itgt]].Name()));
        }
      }
      bool atomsNeedBuilding = false;
      // Loop over template atoms
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        // Track intra-residue bonds from the template.
        Atom templateAtom = resTemplate->Top()[iref];
        int at0 = iref + atomOffset;
        for (Atom::bond_iterator bat = templateAtom.bondbegin(); bat != templateAtom.bondend(); ++bat) {
          int at1 = *bat + atomOffset;
          if (at1 > at0) {
            //mprintf("Will add bond between %i and %i (original %i and %i)\n", at0+1, at1+1, iref+1, *bat + 1);
            intraResBonds.push_back( Ipair(at0, at1) );
          }
        }
        templateAtom.ClearBonds(); // FIXME AddTopAtom should clear bonds
        //mprintf("DEBUG: Adding template %i atom %6s (elt %2s) Res %4s\n", topOut.Natom()+1, templateAtom.c_str(), templateAtom.ElementName(), currentRes.c_str());
        topOut.AddTopAtom( templateAtom, currentRes );
        if (map[iref] == -1) {
          // Template atom not in input structure.
          frameOut.AddVec3( Vec3(0.0) );
          hasPosition.push_back( false );
          nRefAtomsMissing++;
          atomsNeedBuilding = true;
          if (has_bfac) topOut.AddBfactor(0.0);
          if (has_occ)  topOut.AddOccupancy(0.0);
          if (has_pdbn) topOut.AddPdbSerialNum(0);
        } else {
          // Template atom was in input structure.
          int itgt = map[iref];
          frameOut.AddVec3( Vec3(frameIn.XYZ(itgt)) );
          hasPosition.push_back( true );
          if (has_bfac) topOut.AddBfactor( topIn.Bfactor()[itgt] );
          if (has_occ)  topOut.AddOccupancy( topIn.Occupancy()[itgt] );
          if (has_pdbn) topOut.AddPdbSerialNum( topIn.PdbSerialNum()[itgt] );
          //pdb[itgt-currentRes.FirstAtom()] = iref;
          // Check source atoms for inter-residue connections.
/*          Atom const& sourceAtom = topIn[itgt];
          for (Atom::bond_iterator bat = sourceAtom.bondbegin(); bat != sourceAtom.bondend(); ++bat) {
            if ( topIn[*bat].ResNum() < ires ) {
              // Use template atom names. Use saved names in case source name had an alias.
              detectedInterResBonds.push_back( ResAtPair(ires, templateAtom.Name()) );
              detectedInterResBonds.push_back( ResAtPair(topIn[*bat].ResNum(), SourceAtomNames[*bat]) );
              if ( debug_ > 1 ) {
                mprintf("DEBUG: Adding detected interres bond: itgt=%i name=%s res=%i tempName=%s -- bat=%i name=%s res=%i srcName=%s\n",
                        itgt+1, *(sourceAtom.Name()), ires+1, *(templateAtom.Name()),
                        *bat+1, *(topIn[*bat].Name()), topIn[*bat].ResNum()+1, *(SourceAtomNames[*bat]));
              }
            }
          }*/
        }
      } // END loop over template atoms
      if (nTgtAtomsMissing > 0) {
        nAtomsNotInTemplates += nTgtAtomsMissing;
        mprintf("\t%i source atoms not mapped to template.\n", nTgtAtomsMissing);
        if (keepMissingSourceAtoms_) {
          ResAtArray tmpBonds;
          int firstNonTemplateAtom = topOut.Natom();
          for (int itgt = 0; itgt != currentRes.NumAtoms(); itgt++) {
            if (pdb[itgt] == -1) {
              // This PDB atom had no equivalent in the residue template
              int pdbat = itgt + currentRes.FirstAtom();
              Atom pdbAtom = topIn[pdbat];
              if (debug_ > 0)
                mprintf("DEBUG:\t\tInput atom %s missing from template.\n",*(pdbAtom.Name()));
              // Bonds
              for (Atom::bond_iterator bit = pdbAtom.bondbegin(); bit != pdbAtom.bondend(); ++bit)
              {
                Atom bndAt = topIn[*bit];
                NameType bndAtmName;
                if (SourceAtomNames[*bit].len() > 0)
                  bndAtmName = SourceAtomNames[*bit];
                else
                  bndAtmName = bndAt.Name();
                if (debug_ > 0)
                  mprintf("DEBUG:\t\t\tBonded to %s\n", *(bndAtmName));
                tmpBonds.push_back( ResAtPair(ires, pdbAtom.Name()) );
                tmpBonds.push_back( ResAtPair(bndAt.ResNum(), bndAt.Name()) );
              }
              // Add missing atom
              pdbAtom.ClearBonds();
              topOut.AddTopAtom( pdbAtom, currentRes );
              frameOut.AddVec3( Vec3(frameIn.XYZ(pdbat)) );
              hasPosition.push_back( true );
              if (has_bfac)
                topOut.AddBfactor( topIn.Bfactor()[pdbat] );
              if (has_occ)
                topOut.AddOccupancy( topIn.Occupancy()[pdbat] );
              if (has_pdbn)
                topOut.AddPdbSerialNum( topIn.PdbSerialNum()[pdbat] );
              
            }
          } // END loop over input residue atoms
          // Add Bonds
          for (ResAtArray::const_iterator bit = tmpBonds.begin(); bit != tmpBonds.end(); ++bit)
          {
            int ba0 = topOut.FindAtomInResidue( bit->first, bit->second );
            ++bit;
            int ba1 = topOut.FindAtomInResidue( bit->first, bit->second );
            if (ba0 < firstNonTemplateAtom ||
                ba1 < firstNonTemplateAtom ||
                ba0 < ba1)
            {
              if (debug_ > 0)
                mprintf("DEBUG:\t\tAdd bond %i %s -- %i %s\n", ba0+1, *(topOut[ba0].Name()), ba1+1, *(topOut[ba1].Name()));
              topOut.AddBond( ba0, ba1 );
            }
          } // END loop over bonds
        }
      }
      // Save atom offset if atoms need to be built
      if (atomsNeedBuilding)
        AtomOffsets.push_back( atomOffset );
      else
        AtomOffsets.push_back( -1 );
    } // END template exists
    // Add intra-residue bonds
    for (IParray::const_iterator it = intraResBonds.begin(); it != intraResBonds.end(); ++it)
    {
      //mprintf("DEBUG: Intra-res bond: Res %s atom %s (%s) to res %s atom %s (%s)\n",
      //        topOut.TruncResNameOnumId(topOut[it->first].ResNum()).c_str(), *(topOut[it->first].Name()),
      //        topOut.AtomMaskName(it->first).c_str(),
      //        topOut.TruncResNameOnumId(topOut[it->second].ResNum()).c_str(), *(topOut[it->second].Name()),
      //        topOut.AtomMaskName(it->second).c_str());
      topOut.AddBond(it->first, it->second);
    }
  } // END loop over source residues
  t_fill_template_.Stop();
  if (nRefAtomsMissing > 0)
    mprintf("\t%i template atoms missing in source.\n", nRefAtomsMissing);
  if (nAtomsNotInTemplates > 0) {
    if (requireAllInputAtoms_)
      mprinterr("Error: %i input atoms not in templates.\n", nAtomsNotInTemplates);
    else
      mprintf("\t%i input atoms were not in templates and were ignored.\n", nAtomsNotInTemplates);
  }
  if (nAtomsMissingTypes > 0) {
    mprinterr("Error: %i atoms are missing types, either because they did not have\n"
              "Error:  one initially or they could not be matched to a template.\n"
              "Error: This can happen if a parameter file is missing or a force field\n"
              "Error:  file has not been loaded.\n"
              "Error: Build cannot proceed unless all atoms have a type:\n",
              nAtomsMissingTypes);
    std::set<NameType> atoms_missing_types;
    for (int ires = 0; ires != topOut.Nres(); ires++) {
      for (int at = topOut.Res(ires).FirstAtom(); at != topOut.Res(ires).LastAtom(); ++at) {
        if ( !topOut[at].HasType() )
          atoms_missing_types.insert( topOut[at].Name() );
      }
    }
    mprinterr("Error: Atoms missing types:");
    for (std::set<NameType>::const_iterator ait = atoms_missing_types.begin();
                                            ait != atoms_missing_types.end(); ++ait)
      mprinterr(" %s", *(*ait));
    mprinterr("\n");
    if (debug_ > 0) {
      for (int ires = 0; ires != topOut.Nres(); ires++) {
        std::string missingTypes;
        for (int at = topOut.Res(ires).FirstAtom(); at != topOut.Res(ires).LastAtom(); ++at)
          if ( !topOut[at].HasType() )
            missingTypes.append(" " + topOut[at].Name().Truncated() );
        if (!missingTypes.empty())
          mprinterr("Error:\t%s missing types for%s\n", topOut.TruncResNameNum(ires).c_str(), missingTypes.c_str());
      }
    }
    return 1;
  }

  // -----------------------------------
  // DEBUG - Print primary connection atoms
  if (debug_ > 0) {
    for (unsigned int idx = 0; idx != ResTemplates.size(); idx++) {
      if (ResTemplates[idx] != 0) {
        mprintf("DEBUG: Template %s (%s)", ResTemplates[idx]->legend(), Cpptraj::Structure::terminalStr(ResTermTypes[idx]));
        if (resHeadAtoms[idx] > -1) mprintf(" head %s", topOut.AtomMaskName(resHeadAtoms[idx]).c_str());
        if (resTailAtoms[idx] > -1) mprintf(" tail %s", topOut.AtomMaskName(resTailAtoms[idx]).c_str());
        mprintf("\n");
      }
    }
  }

  // Keep track of which residues are connected.
  ResConnectArray ResidueConnections( topOut.Nres() );

  // Try to connect HEAD atoms to previous residue TAIL atoms.
  std::vector<IParray> resBondingAtoms(topOut.Nres());
  for (int ires = 1; ires < topOut.Nres(); ires++) {
    int pres = ires - 1;
    if (resHeadAtoms[ires] != -1) {
      if (ResTermTypes[ires] == Cpptraj::Structure::BEG_TERMINAL) {
        if (debug_ > 0)
          mprintf("DEBUG: Res %s is begin terminal, ignoring head atom.\n",
                  topOut.TruncResNameOnumId(ires).c_str());
      } else {
        if (resTailAtoms[pres] != -1) {
          if (ResTermTypes[pres] == Cpptraj::Structure::END_TERMINAL) {
            if (debug_ > 0)
              mprintf("DEBUG: Res %s is end terminal, ignoring tail atom.\n",
                      topOut.TruncResNameOnumId(pres).c_str());
          } else {
            if (debug_ > 0)
              mprintf("DEBUG: Connecting HEAD atom %s to tail atom %s\n",
                      topOut.AtomMaskName(resHeadAtoms[ires]).c_str(),
                      topOut.AtomMaskName(resTailAtoms[pres]).c_str());
            resBondingAtoms[ires].push_back( Ipair(resHeadAtoms[ires], resTailAtoms[pres]) );
            resBondingAtoms[pres].push_back( Ipair(resTailAtoms[pres], resHeadAtoms[ires]) );
            resHeadAtoms[ires] = -1;
            resTailAtoms[pres] = -1;
            ResidueConnections[ires].push_back( pres );
            ResidueConnections[pres].push_back( ires );
          }
        }
      }
    }
  }

  // Report unused HEAD/TAIL atoms
  for (int ires = 0; ires != topOut.Nres(); ires++) { // TODO should be topIn?
    if (resHeadAtoms[ires] != -1)
      mprintf("Warning: Unused head atom %s\n", topOut.AtomMaskName(resHeadAtoms[ires]).c_str());
    if (resTailAtoms[ires] != -1)
      mprintf("Warning: Unused tail atom %s\n", topOut.AtomMaskName(resTailAtoms[ires]).c_str());
  }

  // Add external bonds
  if (!topInBonds.empty()) {
    for (std::vector<BondType>::const_iterator bnd = topInBonds.begin();
                                               bnd != topInBonds.end(); ++bnd)
    {
      Atom const& At0 = topIn[bnd->A1()];
      Atom const& At1 = topIn[bnd->A2()];
      // Ignore bonds in the same residue
      if ( At0.ResNum() == At1.ResNum()) {
        mprintf("Build Warning: Atoms %s and %s are in the same residue %i. Not adding extra bond.\n",
                *(At0.Name()), *(At1.Name()), At0.ResNum()+1);
        continue;
      }
      int a0 = topOut.FindAtomInResidue( At0.ResNum(), At0.Name() );
      int a1 = topOut.FindAtomInResidue( At1.ResNum(), At1.Name() );
      if (a0 < 0) {
        mprinterr("Error: Atom %s not found in residue %i\n", *(At0.Name()), At0.ResNum()+1);
        return 1;
      }
      if (a1 < 0) {
        mprinterr("Error: Atom %s not found in residue %i\n", *(At1.Name()), At1.ResNum()+1);
        return 1;
      }
      if (resIsConnected(ResidueConnections[At0.ResNum()], At1.ResNum())) {
        mprintf("Warning: Residue %s already connected to residue %s; ignoring\n"
                "Warning:   extra bond %s - %s\n",
                topOut.TruncResNameNum(At0.ResNum()).c_str(),
                topOut.TruncResNameNum(At1.ResNum()).c_str(),
                topOut.AtomMaskName(a0).c_str(),
                topOut.AtomMaskName(a1).c_str());
      } else {
        if (debug_ > 0)
          mprintf("DEBUG: Adding extra bond %s - %s\n",
                  topOut.AtomMaskName(a0).c_str(),
                  topOut.AtomMaskName(a1).c_str());
        resBondingAtoms[At0.ResNum()].push_back( Ipair(a0, a1) );
        resBondingAtoms[At1.ResNum()].push_back( Ipair(a1, a0) );
        ResidueConnections[At0.ResNum()].push_back( At1.ResNum() );
        ResidueConnections[At1.ResNum()].push_back( At0.ResNum() );
      }
    } // END loop over externally passed in bonds
  }

  // Check detected inter-residue bonds
  for (ResAtArray::const_iterator it = detectedInterResBonds.begin();
                                  it != detectedInterResBonds.end(); ++it)
  {
    ResAtPair const& ra0 = *it;
    ++it;
    ResAtPair const& ra1 = *it;
    bool bondInvolvesResWithNoTemplate = ( (ResTemplates[ra0.first] == 0) ||
                                           (ResTemplates[ra1.first] == 0) );
    if (debug_ > 1)
      mprintf("DEBUG: Inter-res bond: Res %i atom %s to res %i atom %s : bondInvolvesResWithNoTemplate=%i\n",
              ra0.first+1, *(ra0.second),
              ra1.first+1, *(ra1.second),
              (int)bondInvolvesResWithNoTemplate);
    int at0 = topOut.FindAtomInResidue(ra0.first, ra0.second);
    if (at0 < 0) {
      mprinterr("Error: Atom %s not found in residue %i\n", *(ra0.second), ra0.first+1);
      return 1;
    }
    int at1 = topOut.FindAtomInResidue(ra1.first, ra1.second);
    if (at1 < 0) {
      mprinterr("Error: Atom %s not found in residue %i\n", *(ra1.second), ra1.first+1);
      return 1;
    }
    // Save detected inter-residue bonding atoms if not already added via
    // template connect atoms. Convention is atom belonging to the current
    // residue is first.
    // NOTE: Only checking at0/at1 here, which should be fine.
    if (!hasBondingPair(resBondingAtoms[ra0.first], Ipair(at0, at1))) {
      // Check if we already have a connection from ra0 to ra1.
      if (resIsConnected(ResidueConnections[ra0.first], ra1.first)) {
        mprintf("Warning: Residue %s already connected to residue %s; ignoring\n"
                "Warning:   potential detected bond %s - %s\n",
                topOut.TruncResNameNum(ra0.first).c_str(),
                topOut.TruncResNameNum(ra1.first).c_str(),
                topOut.AtomMaskName(at0).c_str(),
                topOut.AtomMaskName(at1).c_str());
      } else {
        if (bondInvolvesResWithNoTemplate) {
          mprintf("\tAdding non-template bond %s - %s\n",
                  topOut.AtomMaskName(at0).c_str(),
                  topOut.AtomMaskName(at1).c_str());
          resBondingAtoms[ra0.first].push_back( Ipair(at0, at1) );
          resBondingAtoms[ra1.first].push_back( Ipair(at1, at0) );
          ResidueConnections[ra0.first].push_back( ra1.first );
          ResidueConnections[ra1.first].push_back( ra0.first );
        } else {
          if (debug_ > 0)
            mprintf("DEBUG: Detected non-template bond %s - %s; not adding it.\n",
                    topOut.AtomMaskName(at0).c_str(),
                    topOut.AtomMaskName(at1).c_str());
        }
      }
    }
  }

  // DEBUG print inter-residue bonding atoms
  if (debug_ > 0) {
    for (std::vector<IParray>::const_iterator rit = resBondingAtoms.begin();
                                              rit != resBondingAtoms.end(); ++rit)
    {
      mprintf("\tResidue %s bonding atoms.\n", topOut.TruncResNameNum(rit-resBondingAtoms.begin()).c_str());
      for (IParray::const_iterator it = rit->begin(); it != rit->end(); ++it)
        mprintf("\t\t%s - %s\n", topOut.AtomMaskName(it->first).c_str(), topOut.AtomMaskName(it->second).c_str());
    }
  }

  // -----------------------------------
  // Do some error checking
  if (hasPosition.size() != (unsigned int)newNatom) {
    mprinterr("Internal Error: hasPosition size %zu != newNatom size %i\n", hasPosition.size(), newNatom);
    return 1;
  }
  if (AtomOffsets.size() != (unsigned int)topOut.Nres()) {
    mprinterr("Internal Error: AtomOffsets size %zu != newNres size %i\n", AtomOffsets.size(), topOut.Nres());
    return 1;
  }
  if (SourceAtomNames.size() != (unsigned int)topIn.Natom()) {
    mprinterr("Internal Error: Source atom names length %zu != # input atoms %i\n", SourceAtomNames.size(), topIn.Natom());
    return 1;
  }
  if (has_bfac) {
    if (topOut.Bfactor().size() != (unsigned int)topOut.Natom()) {
      mprinterr("Internal Error: Size of final Bfactor array %zu != # atoms %i\n", topOut.Bfactor().size(), topOut.Natom());
      return 1;
    }
  }
  if (has_occ) {
    if (topOut.Occupancy().size() != (unsigned int)topOut.Natom()) {
      mprinterr("Internal Error: Size of final Occupancy array %zu != # atoms %i\n", topOut.Occupancy().size(), topOut.Natom());
      return 1;
    }
  }
  if (has_pdbn) {
    if (topOut.PdbSerialNum().size() != (unsigned int)topOut.Natom()) {
      mprinterr("Internal Error: Size of final PdbSerialNum array %zu != # atoms %i\n", topOut.PdbSerialNum().size(), topOut.Natom());
      return 1;
    }
  }

  if (debug_ > 0) {
    mprintf("DEBUG: hasPosition:\n");
    for (unsigned int idx = 0; idx != (unsigned int)topOut.Natom(); idx++)
      mprintf("\t%10u %20s : %i\n", idx+1, topOut.AtomMaskName(idx).c_str(), (int)hasPosition[idx]);
  }

  // -----------------------------------
  // Build using internal coords if needed.
  t_fill_build_.Start();
  bool buildFailed = false;
  Cpptraj::Structure::Builder::Barray tmpHasPosition(topOut.Natom(), false);
  int nextTempHasPositionStart = 0;
  for (Iarray::const_iterator it = AtomOffsets.begin(); it != AtomOffsets.end(); ++it)
  {
    long int ires = it-AtomOffsets.begin();
    if (*it > -1) {
      if (debug_ > 0)
        mprintf("DEBUG: ***** BUILD residue %li %s *****\n", ires + 1,
                topOut.TruncResNameOnumId(ires).c_str());
      // Residue has atom offset which indicates it needs something built.
      Cpptraj::Structure::Builder structureBuilder;// = new Cpptraj::Structure::Builder();
      structureBuilder.SetDebug( debug_ );
      if (creator.HasMainParmSet())
        structureBuilder.SetParameters( creator.MainParmSetPtr() );
      // Generate internals from the template, update indices to this topology.
      DataSet_Coords* resTemplate = ResTemplates[ires];
      t_fill_build_internals_.Start();
      Frame templateFrame = resTemplate->AllocateFrame();
      resTemplate->GetFrame( 0, templateFrame );
      if (structureBuilder.GenerateInternals(templateFrame, resTemplate->Top(),
                                             std::vector<bool>(resTemplate->Top().Natom(), true)))
      {
        mprinterr("Error: Generate internals for residue template failed.\n");
        return 1;
      }
      t_fill_build_internals_.Stop();
      structureBuilder.UpdateIndicesWithOffset( *it );
      //mprintf("DEBUG: Residue type: %s terminal\n", Cpptraj::Structure::terminalStr(*termType));
      // Is this residue connected to an earlier residue?
      t_fill_build_link_.Start();
      for (IParray::const_iterator resBonds = resBondingAtoms[ires].begin();
                                   resBonds != resBondingAtoms[ires].end(); ++resBonds)
      {
        if (resBonds->second < resBonds->first) {
          if (debug_ > 0)
            mprintf("\t\tResidue connection: %s - %s\n",
                    topOut.AtomMaskName(resBonds->first).c_str(),
                    topOut.AtomMaskName(resBonds->second).c_str());
          t_fill_build_link_bond_.Start();
          topOut.AddBond(resBonds->first, resBonds->second);
          t_fill_build_link_bond_.Stop();
          // Generate internals around the link
          Residue const& R0 = topIn.Res(topOut[resBonds->first].ResNum());
          for (int at = nextTempHasPositionStart; at < R0.FirstAtom(); at++)
            tmpHasPosition[at] = true;
          nextTempHasPositionStart = R0.FirstAtom();

          if (structureBuilder.GenerateInternalsAroundLink(resBonds->first, resBonds->second,
                                                           frameOut, topOut, hasPosition, Cpptraj::Structure::Builder::BUILD,
                                                           tmpHasPosition))
          {
            mprinterr("Error: Assign torsions around inter-residue link %s - %s failed.\n",
                      topOut.AtomMaskName(resBonds->first).c_str(),
                      topOut.AtomMaskName(resBonds->second).c_str());
            return 1;
          }
        }
      }
      t_fill_build_link_.Stop();
      // Update internal coords from known positions
      if (structureBuilder.UpdateICsFromFrame( frameOut, topOut, hasPosition )) {
        mprinterr("Error: Failed to update internals with values from existing positions.\n");
        return 1;
      }
      t_fill_build_build_.Start();
      if (structureBuilder.BuildFromInternals(frameOut, topOut, hasPosition)) {
        mprinterr("Error: Building residue %s failed.\n",
                  topOut.TruncResNameOnumId(ires).c_str());
        buildFailed = true;
      }
      t_fill_build_build_.Stop();
    } else {
      // All atoms present. Just connect
      // Is this residue connected to an earlier residue?
      for (IParray::const_iterator resBonds = resBondingAtoms[ires].begin();
                                   resBonds != resBondingAtoms[ires].end(); ++resBonds)
      {
        if (resBonds->second < resBonds->first) {
          if (debug_ > 0)
            mprintf("\t\tResidue connection only: %s - %s\n",
                    topOut.AtomMaskName(resBonds->first).c_str(),
                    topOut.AtomMaskName(resBonds->second).c_str());
          topOut.AddBond(resBonds->first, resBonds->second);
        }
      }
    }
  } // END loop over atom offsets
  t_fill_build_.Stop();

  // DEBUG - Print new top/coords
  if (debug_ > 1) {
    for (int iat = 0; iat != topOut.Natom(); iat++)
    {
      Residue const& res = topOut.Res( topOut[iat].ResNum() );
      const double* XYZ = frameOut.XYZ(iat);
      mprintf("%6i %6s %6i %6s (%i) %g %g %g\n",
              iat+1, *(topOut[iat].Name()), res.OriginalResNum(), *(res.Name()),
              (int)hasPosition[iat], XYZ[0], XYZ[1], XYZ[2]);
    }
  }

  if (buildFailed) return 1;
  return 0;
}

/** Given an original topology and bonded atom indices, find those atoms
  * in another topology and ensure they are bonded.
  */
int Exec_Build::transfer_bonds(Topology& topOut, Topology const& topIn,
                               std::vector<BondType> const& bondsIn)
const
{
  for (std::vector<BondType>::const_iterator bnd = bondsIn.begin();
                                             bnd != bondsIn.end(); ++bnd)
  {
    // Get the original atom name and residue number
    Atom const& original_A1 = topIn[bnd->A1()];
    Atom const& original_A2 = topIn[bnd->A2()];
    if (debug_ > 0)
      mprintf("DEBUG: Original bond atoms %i (%s) %i (%s)\n",
              bnd->A1()+1, topIn.AtomMaskName(bnd->A1()).c_str(),
              bnd->A2()+1, topIn.AtomMaskName(bnd->A2()).c_str());
    // Find the atoms in the new topology
    int a1 = topOut.FindAtomInResidue( original_A1.ResNum(), original_A1.Name() );
    if (a1 < 0) {
      mprinterr("Error: Could not find atom %i (%s) in new topology.\n",
                bnd->A1()+1, topIn.AtomMaskName(bnd->A1()).c_str());
      return 1;
    }
    int a2 = topOut.FindAtomInResidue( original_A2.ResNum(), original_A2.Name() );
    if (a2 < 0) {
      mprinterr("Error: Could not find atom %i (%s) in new topology.\n",
                bnd->A2()+1, topIn.AtomMaskName(bnd->A2()).c_str());
      return 1;
    }
    // Add the bond to the new topology
    topOut.AddBond( a1, a2 );
  }
  return 0;
}

// -----------------------------------------------------------------------------
// Exec_Build::Help()
void Exec_Build::Help() const
{
  mprintf("\tname <output COORDS> crdset <COORDS set> [frame <#>]\n"
          "\t[title <title>] [gb <radii>] [verbose <#>] [keepmissingatoms]\n"
          "\t[parmout <topology file>] [crdout <coord file>] [simplecheck]\n"
          "\t[%s]\n"
          "\t[{%s} ...]\n"
          "\t[{%s} ...]\n"
          "\t[{{solvatebox|solvateoct} %s\n"
          "\t          %s |\n"
          "\t  setbox %s}]\n"
          "%s"
          "%s",
          Cpptraj::Structure::Creator::other_keywords_,
          Cpptraj::Structure::Creator::template_keywords_,
          Cpptraj::Structure::Creator::parm_keywords_,
          Cpptraj::Structure::Solvate::SolvateKeywords1(),
          Cpptraj::Structure::Solvate::SolvateKeywords2(),
          Cpptraj::Structure::Solvate::SetboxKeywords(),
          Cpptraj::Structure::HisProt::keywords_,
          Cpptraj::Structure::Disulfide::keywords_
         );
  Cpptraj::Parm::PrintGbRadiiKeywords();
  mprintf("  Build complete topology and parameters from given crdset.\n");
}

// Exec_Build::Execute()
Exec::RetType Exec_Build::Execute(CpptrajState& State, ArgList& argIn)
{
  // Get input coords
  std::string crdset = argIn.GetStringKey("crdset");
  if (crdset.empty()) {
    mprinterr("Error: Must specify input COORDS set with 'crdset'\n");
    return CpptrajState::ERR;
  }
  DataSet* inCrdPtr = State.DSL().FindSetOfGroup( crdset, DataSet::COORDINATES );
  if (inCrdPtr == 0) {
    mprinterr("Error: No COORDS set found matching %s\n", crdset.c_str());
    return CpptrajState::ERR;
  }
  std::string outset = argIn.GetStringKey("name");
  if (outset.empty()) {
    mprinterr("Error: Must specify output COORDS set with 'name'\n");
    return CpptrajState::ERR;
  }

  return BuildStructure(inCrdPtr, outset, State.DSL(), State.Debug(), argIn, Cpptraj::Parm::UNKNOWN_GB);
}

/** Standalone execute. For DataIO_LeapRC. Operate on inCrdPtr */
Exec::RetType Exec_Build::BuildStructure(DataSet* inCrdPtr, 
                                         DataSetList& DSL, int debugIn, ArgList& argIn,
                                         Cpptraj::Parm::GB_RadiiType gbRadIn)
{
  return BuildStructure(inCrdPtr, "", DSL, debugIn, argIn, gbRadIn);
}

/** Standalone execute. For DataIO_LeapRC. */
Exec::RetType Exec_Build::BuildStructure(DataSet* inCrdPtr, std::string const& outset,
                                         DataSetList& DSL, int debugIn, ArgList& argIn,
                                         Cpptraj::Parm::GB_RadiiType gbRadIn)
{
  t_total_.Start();
  if (inCrdPtr == 0) {
    mprinterr("Internal Error: Exec_Build::BuildStructure(): Null input coordinates.\n");
    return CpptrajState::ERR;
  }
  if (inCrdPtr->Group() != DataSet::COORDINATES) {
    mprinterr("Error: Set '%s' is not coordinates, cannot use for building.\n", inCrdPtr->legend());
    return CpptrajState::ERR;
  }
  debug_ = debugIn;
  int verbose = argIn.getKeyInt("verbose", 0);
  std::string title = argIn.GetStringKey("title");
  std::string outputTopologyName = argIn.GetStringKey("parmout");
  std::string outputCoordsName = argIn.GetStringKey("crdout");
  if (argIn.hasKey("simplecheck")) {
    mprintf("\tSimple check: will only check bond lengths.\n");
    check_structure_ = false;
  } else {
    mprintf("\tWill check bond lengths, atomic overlaps, and ring intersections.\n");
    check_structure_ = true;
  }
  if (!outputTopologyName.empty())
    mprintf("\tWill write topology to %s\n", outputTopologyName.c_str());
  if (!outputCoordsName.empty())
    mprintf("\tWill write coords to %s\n", outputCoordsName.c_str());
  // TODO make it so this can be const (cant bc GetFrame)
  DataSet_Coords& coords = static_cast<DataSet_Coords&>( *((DataSet_Coords*)inCrdPtr) );
  // Get frame from input coords
  int tgtframe = argIn.getKeyInt("frame", 1) - 1;
  mprintf("\tUsing frame %i from COORDS set %s\n", tgtframe+1, coords.legend());
  if (tgtframe < 0 || tgtframe >= (int)coords.Size()) {
    mprinterr("Error: Frame is out of range.\n");
    return CpptrajState::ERR;
  }
  Frame frameIn = coords.AllocateFrame();
  coords.GetFrame(tgtframe, frameIn);
  // Get modifiable topology
  //Topology& topIn = *(coords.TopPtr());
  //Topology const& topIn = coords.Top();
  Topology topIn = coords.Top(); // FIXME do not work on the copy, work on the top itself

  std::string solventResName = argIn.GetStringKey("solventresname", "HOH");
  mprintf("\tResidues named '%s' will be recognized as solvent.\n", solventResName.c_str());

  enum SolvateModeType { NO_SOLVATE = 0, SOLVATEBOX, SETBOX };
  SolvateModeType add_solvent = NO_SOLVATE;
  Cpptraj::Structure::Solvate solvator;
  if (argIn.hasKey("solvatebox")) {
    add_solvent = SOLVATEBOX;
    if (solvator.InitSolvate(argIn, false, debug_)) {
      mprinterr("Error: Init solvatebox failed.\n");
      return CpptrajState::ERR;
    }
    solvator.PrintSolvateInfo();
  } else if (argIn.hasKey("solvateoct")) {
    add_solvent = SOLVATEBOX;
    if (solvator.InitSolvate(argIn, true, debug_)) {
      mprinterr("Error: Init solvateoct failed.\n");
      return CpptrajState::ERR;
    }
    solvator.PrintSolvateInfo();
  } else if (argIn.hasKey("setbox")) {
    add_solvent = SETBOX;
    if (solvator.InitSetbox(argIn, debug_)) {
      mprinterr("Error: Init setbox failed.\n");
      return CpptrajState::ERR;
    }
  }

  Cpptraj::Structure::AddIons addIons;
  if (add_solvent == SOLVATEBOX) {
    if (argIn.hasKey("addionsrand")) {
      std::string ion1name = argIn.GetStringKey("ion1");
      int nion1 = argIn.getKeyInt("nion1", -1);
      std::string ion2name = argIn.GetStringKey("ion2");
      int nion2 = argIn.getKeyInt("nion2", -1);
      double minsep = argIn.getKeyDouble("minsep", 0.0);
      int ionseed = argIn.getKeyInt("ionseed", -1);
      if (addIons.InitAddIons(ion1name, nion1, ion2name, nion2, minsep, ionseed, debug_)) {
        mprinterr("Error: Init addions failed.\n");
        return CpptrajState::ERR;
      }
      addIons.PrintAddIonsInfo();
    }
  }

  keepMissingSourceAtoms_ = argIn.hasKey("keepmissingatoms");
  if (keepMissingSourceAtoms_)
    mprintf("\tInput atoms missing from templates will be kept.\n");
  else
    mprintf("\tInput atoms missing from templates will be ignored.\n");
  requireAllInputAtoms_ = argIn.hasKey("requireallinputatoms");
  if (requireAllInputAtoms_)
    mprintf("\tRequire all input atoms to be found in templates.\n");
  else
    mprintf("\tInput atoms not found in templates will be ignored.\n");

  // Do histidine detection before H atoms are removed
  t_hisDetect_.Start();
  if (!argIn.hasKey("nohisdetect")) {
    Cpptraj::Structure::HisProt hisProt;
    if (hisProt.InitHisProt( argIn, debug_ )) {
      mprinterr("Error: Could not initialize histidine detection.\n");
      return CpptrajState::ERR;
    }
    hisProt.HisProtInfo();
    if (hisProt.DetermineHisProt( topIn )) {
      mprinterr("Error: HIS protonation detection failed.\n");
      return CpptrajState::ERR;
    }
  }
  t_hisDetect_.Stop();

  // Clean up structure
  t_clean_.Start();
  Cpptraj::Structure::PdbCleaner pdbCleaner;
  pdbCleaner.SetDebug( debug_ );
  if (pdbCleaner.InitPdbCleaner( argIn, solventResName, std::vector<int>() )) {
    mprinterr("Error: Could not init PDB cleaner.\n");
    return CpptrajState::ERR;
  }
  if (pdbCleaner.SetupPdbCleaner( topIn )) {
    mprinterr("Error: Could not set up PDB cleaner.\n");
    return CpptrajState::ERR;
  }
  pdbCleaner.PdbCleanerInfo();
  if (pdbCleaner.ModifyCoords(topIn, frameIn)) {
    mprinterr("Error: Could not clean PDB.\n");
    return CpptrajState::ERR;
  }
  t_clean_.Stop();

  // Set up Output coords
  if (!outset.empty()) {
    // Separate output COORDS set.
    outCrdPtr_ = DSL.AddSet( DataSet::COORDS, outset );
  } else {
    // In-place output COORDS
    DSL.RemoveSet( inCrdPtr );
    outCrdPtr_ = DSL.AddSet( DataSet::COORDS, inCrdPtr->Meta() );
  }
  if (outCrdPtr_ == 0) {
    mprinterr("Error: Could not allocate output COORDS set with name '%s'\n", outset.c_str());
    return CpptrajState::ERR;
  }
  DataSet_Coords& crdout = static_cast<DataSet_Coords&>( *((DataSet_Coords*)outCrdPtr_) );
  mprintf("\tOutput COORDS set: %s\n", crdout.legend());

  // GB radii set
  Cpptraj::Parm::GB_RadiiType gbradii;
  if (gbRadIn == Cpptraj::Parm::UNKNOWN_GB) {
    // No radii set specified. Check for keyword.
    gbradii = Cpptraj::Parm::MBONDI; // Default
    std::string gbset = argIn.GetStringKey("gb");
    if (!gbset.empty()) {
      gbradii = Cpptraj::Parm::GbTypeFromKey( gbset );
      if (gbradii == Cpptraj::Parm::UNKNOWN_GB) {
        mprinterr("Error: Unknown GB radii set: %s\n", gbset.c_str());
        return CpptrajState::ERR;
      }
    }
  } else {
    // Use passed-in GB radii set
    gbradii = gbRadIn;
  }
  mprintf("\tGB radii set: %s\n", Cpptraj::Parm::GbTypeStr(gbradii).c_str());

  // LJ 12-6-4
  Cpptraj::Parm::LJ1264_Params lj1264;
  if (argIn.hasKey("lj1264")) {
    std::string lj1264mask = argIn.GetStringKey("lj1264mask");
    std::string c4file = argIn.GetStringKey("c4file");
    std::string polfile = argIn.GetStringKey("polfile");
    double tunfactor = argIn.getKeyDouble("tunfactor", 1.0);
    // Try to guess the water model if we are solvating
    Cpptraj::Parm::WaterModelType wm = Cpptraj::Parm::UNKNOWN_WATER_MODEL;
    if (add_solvent == SOLVATEBOX) {
      if (solvator.SolventBoxName() == "TIP3PBOX") wm = Cpptraj::Parm::TIP3P;
      else if (solvator.SolventBoxName() == "TIP4PEWBOX") wm = Cpptraj::Parm::TIP4PEW;
      else if (solvator.SolventBoxName() == "SPCBOX") wm = Cpptraj::Parm::SPCE;
      else if (solvator.SolventBoxName() == "OPC3BOX") wm = Cpptraj::Parm::OPC3;
      else if (solvator.SolventBoxName() == "OPCBOX") wm = Cpptraj::Parm::OPC;
      else if (solvator.SolventBoxName() == "FB3BOX") wm = Cpptraj::Parm::FB3;
      else if (solvator.SolventBoxName() == "FB4BOX") wm = Cpptraj::Parm::FB4;
      else
        mprintf("Warning: Unable to determine solvent model for LJ 12-6-4 from solvent box name %s\n",
                solvator.SolventBoxName().c_str());
    }
    if (wm == Cpptraj::Parm::UNKNOWN_WATER_MODEL) {
      mprintf("Warning: Unable to determine water model for LJ 12-6-4. Using TIP3P.\n");
      wm = Cpptraj::Parm::TIP3P;
    }
    if (lj1264.Init_LJ1264(lj1264mask, c4file, wm, polfile, tunfactor, debug_)) {
      mprinterr("Error: Init of LJ 12-6-4 failed.\n");
      return CpptrajState::ERR;
    }
  }

  // Get templates and parameter sets.
  t_get_templates_.Start();
  Cpptraj::Structure::Creator creator( debug_ );
  if (creator.InitCreator(argIn, DSL, debug_)) {
    return CpptrajState::ERR;
  }
  if (!creator.HasTemplates()) {
    mprintf("Warning: No residue templates loaded.\n");
  }
  if (!creator.HasMainParmSet()) {
    mprinterr("Error: No parameter sets.\n");
    return CpptrajState::ERR;
  }
  t_get_templates_.Stop();
  // FIXME hide behind ifdef?
  creator.TimingInfo(t_get_templates_.Total(), 2);

  // All residues start unknown
  Cpptraj::Structure::ResStatArray resStat( topIn.Nres() );
  std::vector<BondType> DisulfideBonds;
  std::vector<BondType> SugarBonds;

  // Disulfide search
  t_disulfide_.Start();
  if (!argIn.hasKey("nodisulfides")) {
    Cpptraj::Structure::Disulfide disulfide;
    if (disulfide.InitDisulfide( argIn, Cpptraj::Structure::Disulfide::ADD_BONDS, debug_ )) {
      mprinterr("Error: Could not init disulfide search.\n");
      return CpptrajState::ERR;
    }
    if (disulfide.SearchForDisulfides( resStat, topIn, frameIn, DisulfideBonds ))
    {
      mprinterr("Error: Disulfide search failed.\n");
      return CpptrajState::ERR;
    }
  } else {
    mprintf("\tNot searching for disulfides.\n");
  }
  t_disulfide_.Stop();

  // Handle sugars.
  t_sugar_.Start();
  // TODO should be on a residue by residue basis in FillAtomsWithTemplates
  bool prepare_sugars = !argIn.hasKey("nosugars");
  if (!prepare_sugars)
    mprintf("\tNot attempting to prepare sugars.\n");
  else
    mprintf("\tWill attempt to prepare sugars.\n");
  Cpptraj::Structure::SugarBuilder sugarBuilder(debug_);
  if (prepare_sugars) {
    // Init options
    if (sugarBuilder.InitOptions( argIn.hasKey("hasglycam"),
                                  argIn.getKeyDouble("rescut", 8.0),
                                  argIn.getKeyDouble("bondoffset", 0.2),
                                  argIn.GetStringKey("sugarmask"),
                                  argIn.GetStringKey("determinesugarsby", "geometry"),
                                  argIn.GetStringKey("resmapfile") ))
    {
      mprinterr("Error: Sugar options init failed.\n");
      return CpptrajState::ERR;
    }
    bool splitres = !argIn.hasKey("nosplitres");
    if (splitres)
      mprintf("\tWill split off recognized sugar functional groups into separate residues.\n");
    else
      mprintf("\tNot splitting recognized sugar functional groups into separate residues.\n");
    bool c1bondsearch = !argIn.hasKey("noc1search");
    if (c1bondsearch)
      mprintf("\tWill search for missing bonds to sugar anomeric atoms.\n");
    else
      mprintf("\tNot searching for missing bonds to sugar anomeric atoms.\n");
    // May need to modify sugar structure/topology, either by splitting
    // C1 hydroxyls of terminal sugars into ROH residues, and/or by
    // adding missing bonds to C1 atoms.
    // This is done before any identification takes place since we want
    // to identify based on the most up-to-date topology.
    if (sugarBuilder.FixSugarsStructure(topIn, frameIn,
                                        c1bondsearch, splitres, solventResName,
                                        SugarBonds))
    {
      mprinterr("Error: Sugar structure modification failed.\n");
      return CpptrajState::ERR;
    }
    if (sugarBuilder.PrepareSugars(true, resStat, topIn, frameIn, SugarBonds))
    {
      mprinterr("Error: Sugar preparation failed.\n");
      return CpptrajState::ERR;
    }
  }
  t_sugar_.Stop();

  // Fill in atoms with templates
  t_fill_.Start();
  //Topology topOut;
  Topology& topOut = static_cast<Topology&>( *(crdout.TopPtr()) );
  topOut.SetDebug( debug_ );
  if (outset.empty()) {
    // In-place output COORDS. Copy over existing topology metadata
    topOut.CopyTopMetadata( topIn );
  }
  if (!title.empty())
    topOut.SetParmName( title, FileName() );
  else if (topOut.ParmName().empty()) {
    // TODO better default
    title.assign( topIn.c_str() );
    topOut.SetParmName( title, FileName() );
  }
  Frame frameOut;
  if (FillAtomsWithTemplates(topOut, frameOut, topIn, frameIn, creator, SugarBonds)) {
    mprinterr("Error: Could not fill in atoms using templates.\n");
    return CpptrajState::ERR;
  }
  t_fill_.Stop();

  // Add the disulfide/sugar bonds
  int addBondsErr = transfer_bonds( topOut, topIn, DisulfideBonds );
  //addBondsErr += transfer_bonds( topOut, topIn, SugarBonds );
  if (addBondsErr != 0) {
    mprinterr("Error: Adding disulfide/sugar bonds failed.\n");
    return CpptrajState::ERR;
  }

  // Create empty arrays for the TREE, JOIN, and IROTAT arrays
  topOut.AllocTreeChainClassification( );
  topOut.AllocJoinArray();
  topOut.AllocRotateArray();
    // Finalize topology - determine molecules, dont renumber residues
  topOut.CommonSetup(true, false);

  // Solvate/add ions
  std::string checkMaskString("*");
  if (add_solvent == SOLVATEBOX) {
    t_solvate_.Start();
    // Record initial number of atoms and residues
    int initial_natom = topOut.Natom();
    int initial_nres = topOut.Nres();
    // Get solvent unit box
    DataSet_Coords* solventUnitBox = solvator.GetSolventUnit( DSL );
    if (solventUnitBox == 0) {
      mprinterr("Error: Getting solvent unit failed.\n");
      return CpptrajState::ERR;
    }
    if (solvator.SolvateBox( topOut, frameOut, *(creator.MainParmSetPtr()), *solventUnitBox )) {
      mprinterr("Error: Adding solvent failed.\n");
      return CpptrajState::ERR;
    }
    mprintf("\t  Added %i solvent residues.\n", topOut.Nres() - initial_nres);
    // Since by design the SolvateBox routine ensures solvent does not
    // clash with solute, just check the solute.
    checkMaskString.assign( "@1-" + integerToString(initial_natom) );
    // Add ions if needed
    if (addIons.IsSetup()) {
      if (addIons.AddIonsRand( topOut, frameOut, DSL, *(creator.MainParmSetPtr()) )) {
        mprinterr("Error: Adding ions failed.\n");
        return CpptrajState::ERR;
      }
    }
    t_solvate_.Stop();
  } else if (add_solvent == SETBOX) {
    t_solvate_.Start();
    if (solvator.SetVdwBoundingBox( topOut, frameOut, *(creator.MainParmSetPtr()) )) {
      mprinterr("Error: Setting box failed.\n");
      return CpptrajState::ERR;
    }
    mprintf("\tAdding VDW bounding box.\n");
  }

  // Assign parameters. This will create the bond/angle/dihedral/improper
  // arrays as well.
  t_assign_.Start();
  Exec::RetType ret = CpptrajState::OK;
  Cpptraj::Parm::AssignParams AP;
  AP.SetDebug( debug_ );
  AP.SetVerbose( verbose );
  if ( AP.AssignParameters( topOut, *(creator.MainParmSetPtr()) ) ) {
    mprinterr("Error: Could not assign parameters for '%s'.\n", topOut.c_str());
    ret = CpptrajState::ERR;
  }
  // Assign GB parameters
  if (Cpptraj::Parm::Assign_GB_Radii( topOut, gbradii )) {
    mprinterr("Error: Could not assign GB parameters for '%s'\n", topOut.c_str());
    ret = CpptrajState::ERR;
  }
  // Assign LJ 12-6-4 parameters
  if (lj1264.HasC4Params()) {
    if (lj1264.AssignLJ1264(topOut)) {
      mprinterr("Error: Could not assign LJ 12-6-4 parameters for '%s'\n", topOut.c_str());
      ret = CpptrajState::ERR;
    }
  }

  topOut.Summary();
  t_assign_.Stop();

  // Update coords 
  if (crdout.CoordsSetup( topOut, frameOut.CoordsInfo() )) { // FIXME better coordinate info
    mprinterr("Error: Could not set up output COORDS.\n");
    return CpptrajState::ERR;
  }
  crdout.SetCRD(0, frameOut);

  // Structure check
  if (check_structure_) {
    t_check_.Start();
    StructureCheck check;
    if (check.SetOptions( true, // image 
                          true, // check bonds
                          true,  // save problems
                          debug_, // debug
                          checkMaskString, // mask 1
                          "", // mask 2
                          0.8, // nonbond cut. NOTE: leap check cut is 1.5
                          1.15, // bond long offset
                          0.5, // bond short offset
                          -1, // pairlist cut (-1 for heuristic)
                          true, // ring check
                          0, // 0 = default ring check short distance cut
                          0, // 0 = default ring check long distance cut
                          0  // 0 = default ring check angle cut
        ))
    {
      mprinterr("Error: Structure check options failed.\n");
      return CpptrajState::ERR;
    }
    // For larger structures, automatically add a box
    bool box_added = false;
    if (!frameOut.BoxCrd().HasBox() && topOut.Natom() > check_box_natom_) {
      mprintf("\tAdding unit cell for check only.\n");
      box_added = true;
      // Get radii
      std::vector<double> Radii;
      Radii.reserve( topOut.Natom() );
      for (int atnum = 0; atnum != topOut.Natom(); ++atnum) {
        Radii.push_back( topOut.GetVDWradius(atnum) );
        //Radii.push_back( topIn[atnum].ParseRadius() );
        //Radii.push_back( 0.5 );
      }

      if (frameOut.SetOrthoBoundingBox(Radii, 1.0)) {
        mprinterr("Error: Setting orthogonal bounding box failed.\n");
        return CpptrajState::ERR;
      }
      frameOut.BoxCrd().PrintInfo();
    }
    t_check_setup_.Start();
    if (check.Setup( topOut, frameOut.BoxCrd() )) {
      mprinterr("Error: Structure check setup failed.\n");
      return CpptrajState::ERR;
    }
    t_check_setup_.Stop();
    check.PrintTiming(1, t_check_setup_.Total());
    check.Mask1().MaskInfo();
    if (check.ImageOpt().ImagingEnabled())
      mprintf("\tImaging on.\n");
    else
      mprintf("\timaging off.\n");
    // TODO make file a user option
    CpptrajFile check_output;
    check_output.OpenWrite("");
    t_check_overlaps_.Start();
    int Ntotal_problems = check.CheckOverlaps( frameOut );
    t_check_overlaps_.Stop();
    check.WriteProblemsToFile( &check_output, 1, topOut );
    t_check_bonds_.Start();
    Ntotal_problems += check.CheckBonds( frameOut );
    t_check_bonds_.Stop();
    check.WriteProblemsToFile( &check_output, 1, topOut );
    t_check_rings_.Start();
    Ntotal_problems += check.CheckRings( frameOut );
    t_check_rings_.Stop();
    check.WriteProblemsToFile( &check_output, 1, topOut );
    mprintf("\t%i total problems detected.\n", Ntotal_problems);
    // If box was added for check only, remove it
    if (box_added)
      frameOut.ModifyBox().SetNoBox();
    t_check_.Stop();
  } else {
    // Just check bond lengths
    t_check_.Start();
    StructureCheck check;
    if (check.SetOptions( false, // image 
                          true,  // check bonds
                          true,  // save problems
                          debug_, // debug
                          checkMaskString, // mask 1
                          "", // mask 2
                          0.8, // nonbond cut. NOTE: leap check cut is 1.5
                          1.15, // bond long offset
                          0.5, // bond short offset
                          -1, // pairlist cut (-1 for heuristic)
                          false, // ring check
                          0, // 0 = default ring check short distance cut
                          0, // 0 = default ring check long distance cut
                          0  // 0 = default ring check angle cut
        ))
    {
      mprinterr("Error: Structure check options failed.\n");
      return CpptrajState::ERR;
    }
    t_check_setup_.Start();
    if (check.Setup( topOut, frameOut.BoxCrd() )) {
      mprinterr("Error: Structure check setup failed.\n");
      return CpptrajState::ERR;
    }
    t_check_setup_.Stop();
    check.PrintTiming(1, t_check_setup_.Total());
    check.Mask1().MaskInfo();
    // TODO make file a user option
    CpptrajFile check_output;
    check_output.OpenWrite("");
    t_check_bonds_.Start();
    int Ntotal_problems = check.CheckBonds( frameOut );
    t_check_bonds_.Stop();
    check.WriteProblemsToFile( &check_output, 1, topOut );
    mprintf("\t%i total problems detected.\n", Ntotal_problems);
    t_check_.Stop();
  }

  // Total charge check
  double totalSystemCharge = crdout.Top().TotalCharge();
  double absSystemCharge = fabs(totalSystemCharge);
  double q_frac = fabs( absSystemCharge - (double)(int)(absSystemCharge+0.5) );
  // NOTE: These cutoffs are the same as in LEaP
  if ( q_frac > 0.01 )
   mprintf("Warning: The charge of the system is not integral: %f\n", totalSystemCharge);
  if ( absSystemCharge > 0.01 )
    mprintf("Warning: The charge of the system is not zero: %f\n", totalSystemCharge);

  if (!outputTopologyName.empty()) {
    ParmFile pfile;
    if (pfile.WriteTopology(crdout.Top(), outputTopologyName, argIn, ParmFile::UNKNOWN_PARM, debug_)) {
      mprinterr("Error: Could not write topology file %s\n", outputTopologyName.c_str());
      return CpptrajState::ERR;
    }
  }

  if (!outputCoordsName.empty()) {
    Trajout_Single outtraj;
    if (outtraj.PrepareTrajWrite( outputCoordsName, argIn, DSL, crdout.TopPtr(), crdout.CoordsInfo(),
                                  crdout.Size(), TrajectoryFile::UNKNOWN_TRAJ))
    {
      mprinterr("Error: Could not set up output coords file %s\n", outputCoordsName.c_str());
      return CpptrajState::ERR;
    }
    outtraj.PrintInfo(0);
    if ( outtraj.WriteSingle( 0, frameOut ) ) {
      mprinterr("Error: Could not write output coords file %s\n", outputCoordsName.c_str());
      return CpptrajState::ERR;
    }
  }
  t_total_.Stop();

  PrintTiming();
  if (add_solvent)
    t_solvate_.WriteTiming    (2, "Solvate             :", t_total_.Total());
  t_check_.WriteTiming        (2, "Structure check     :", t_total_.Total());
  t_check_overlaps_.WriteTiming(3, "Overlaps :", t_check_.Total());
  t_check_bonds_.WriteTiming   (3, "Bonds    :", t_check_.Total());
  t_check_rings_.WriteTiming   (3, "Rings    :", t_check_.Total());
  t_check_setup_.WriteTiming   (3, "Setup    :", t_check_.Total());
  t_assign_.WriteTiming       (2, "Param./Top. gen.    :", t_total_.Total());
  AP.WriteAssignTiming(3, t_assign_.Total());

  return ret;
}

/** Standalone execute. For DataIO_LeapRC. */
int Exec_Build::CleanAndFillStructure(DataSet* inCrdPtr, int tgtframe, std::string const& outset, std::string const& title,
                                      DataSetList& DSL, int debugIn, Cpptraj::Structure::Creator const& creator)
{
  t_total_.Start();
  if (inCrdPtr == 0) {
    mprinterr("Internal Error: Exec_Build::BuildStructure(): Null input coordinates.\n");
    return 1;
  }
  if (inCrdPtr->Group() != DataSet::COORDINATES) {
    mprinterr("Error: Set '%s' is not coordinates, cannot use for building.\n", inCrdPtr->legend());
    return 1;
  }
  debug_ = debugIn;
//  int verbose = argIn.getKeyInt("verbose", 0);
//  std::string title = argIn.GetStringKey("title");
//  std::string outputTopologyName = argIn.GetStringKey("parmout");
//  std::string outputCoordsName = argIn.GetStringKey("crdout");
//  if (!outputTopologyName.empty())
//    mprintf("\tWill write topology to %s\n", outputTopologyName.c_str());
//  if (!outputCoordsName.empty())
//    mprintf("\tWill write coords to %s\n", outputCoordsName.c_str());
  // TODO make it so this can be const (cant bc GetFrame)
  DataSet_Coords& coords = static_cast<DataSet_Coords&>( *((DataSet_Coords*)inCrdPtr) );
  // Get frame from input coords
  //int tgtframe = argIn.getKeyInt("frame", 1) - 1;
  mprintf("\tUsing frame %i from COORDS set %s\n", tgtframe+1, coords.legend());
  if (tgtframe < 0 || tgtframe >= (int)coords.Size()) {
    mprinterr("Error: Frame is out of range.\n");
    return 1;
  }
  Frame frameIn = coords.AllocateFrame();
  coords.GetFrame(tgtframe, frameIn);
  // Get modifiable topology
  //Topology& topIn = *(coords.TopPtr());
  //Topology const& topIn = coords.Top();
  Topology topIn = coords.Top(); // FIXME do not work on the copy, work on the top itself

  //std::string solventResName = argIn.GetStringKey("solventresname", "HOH");
  //mprintf("\tSolvent residue name: %s\n", solventResName.c_str());

//  enum SolvateModeType { NO_SOLVATE = 0, SOLVATEBOX, SETBOX };
//  SolvateModeType add_solvent = NO_SOLVATE;
//  Cpptraj::Structure::Solvate solvator;
//  if (argIn.hasKey("solvatebox")) {
//    add_solvent = SOLVATEBOX;
//    if (solvator.InitSolvate(argIn, false, debug_)) {
//      mprinterr("Error: Init solvatebox failed.\n");
//      return CpptrajState::ERR;
//    }
//    solvator.PrintSolvateInfo();
//  } else if (argIn.hasKey("solvateoct")) {
//    add_solvent = SOLVATEBOX;
//    if (solvator.InitSolvate(argIn, true, debug_)) {
//      mprinterr("Error: Init solvateoct failed.\n");
//      return CpptrajState::ERR;
//    }
//    solvator.PrintSolvateInfo();
//  } else if (argIn.hasKey("setbox")) {
//    add_solvent = SETBOX;
//    if (solvator.InitSetbox(argIn, debug_)) {
//      mprinterr("Error: Init setbox failed.\n");
//      return CpptrajState::ERR;
//    }
//  }

  keepMissingSourceAtoms_ = false; // FIXME
//  keepMissingSourceAtoms_ = argIn.hasKey("keepmissingatoms");
  if (keepMissingSourceAtoms_)
    mprintf("\tInput atoms missing from templates will be kept.\n");
  else
    mprintf("\tInput atoms missing from templates will be ignored.\n");
  requireAllInputAtoms_ = false; // FIXME
//  requireAllInputAtoms_ = argIn.hasKey("requireallinputatoms");
  if (requireAllInputAtoms_)
    mprintf("\tRequire all input atoms to be found in templates.\n");
  else
    mprintf("\tInput atoms not found in templates will be ignored.\n");

  // Do histidine detection before H atoms are removed
  ArgList argIn; // FIXME
  t_hisDetect_.Start();
//  if (!argIn.hasKey("nohisdetect")) {
    Cpptraj::Structure::HisProt hisProt;
    if (hisProt.InitHisProt( argIn, debug_ )) {
      mprinterr("Error: Could not initialize histidine detection.\n");
      return 1;
    }
    hisProt.HisProtInfo();
    if (hisProt.DetermineHisProt( topIn )) {
      mprinterr("Error: HIS protonation detection failed.\n");
      return 1;
    }
//  }
  t_hisDetect_.Stop();

  // Clean up structure
  t_clean_.Start();
  Cpptraj::Structure::PdbCleaner pdbCleaner;
  pdbCleaner.SetDebug( debug_ );
  std::string solventResName = "HOH"; //FIXME
  if (pdbCleaner.InitPdbCleaner( argIn, solventResName, std::vector<int>() )) {
    mprinterr("Error: Could not init PDB cleaner.\n");
    return 1;
  }
  if (pdbCleaner.SetupPdbCleaner( topIn )) {
    mprinterr("Error: Could not set up PDB cleaner.\n");
    return 1;
  }
  pdbCleaner.PdbCleanerInfo();
  if (pdbCleaner.ModifyCoords(topIn, frameIn)) {
    mprinterr("Error: Could not clean PDB.\n");
    return 1;
  }
  t_clean_.Stop();

  // Set up Output coords
  if (!outset.empty()) {
    // Separate output COORDS set.
    outCrdPtr_ = DSL.AddSet( DataSet::COORDS, outset );
  } else {
    // In-place output COORDS
    DSL.RemoveSet( inCrdPtr );
    outCrdPtr_ = DSL.AddSet( DataSet::COORDS, inCrdPtr->Meta() );
  }
  if (outCrdPtr_ == 0) {
    mprinterr("Error: Could not allocate output COORDS set with name '%s'\n", outset.c_str());
    return 1;
  }
  DataSet_Coords& crdout = static_cast<DataSet_Coords&>( *((DataSet_Coords*)outCrdPtr_) );
  mprintf("\tOutput COORDS set: %s\n", crdout.legend());
/*
  // GB radii set
  Cpptraj::Parm::GB_RadiiType gbradii;
  if (gbRadIn == Cpptraj::Parm::UNKNOWN_GB) {
    // No radii set specified. Check for keyword.
    gbradii = Cpptraj::Parm::MBONDI; // Default
    //std::string gbset = argIn.GetStringKey("gb");
    //if (!gbset.empty()) {
    //  gbradii = Cpptraj::Parm::GbTypeFromKey( gbset );
    //  if (gbradii == Cpptraj::Parm::UNKNOWN_GB) {
    //    mprinterr("Error: Unknown GB radii set: %s\n", gbset.c_str());
    //    return CpptrajState::ERR;
    //  }
    //}
  } else {
    // Use passed-in GB radii set
    gbradii = gbRadIn;
  }
  mprintf("\tGB radii set: %s\n", Cpptraj::Parm::GbTypeStr(gbradii).c_str());
*/
  // Get templates and parameter sets.
/*  t_get_templates_.Start();
  Cpptraj::Structure::Creator creator( debug_ );
  if (creator.InitCreator(argIn, DSL, debug_)) {
    return 1;
  }
  if (!creator.HasTemplates()) {
    mprintf("Warning: No residue templates loaded.\n");
  }
  if (!creator.HasMainParmSet()) {
    mprinterr("Error: No parameter sets.\n");
    return 1;
  }
  t_get_templates_.Stop();
  // FIXME hide behind ifdef?
  creator.TimingInfo(t_get_templates_.Total(), 2);*/

  // All residues start unknown
  Cpptraj::Structure::ResStatArray resStat( topIn.Nres() );
  std::vector<BondType> DisulfideBonds;
  std::vector<BondType> SugarBonds;

  // Disulfide search
  t_disulfide_.Start();
  //if (!argIn.hasKey("nodisulfides")) {
    Cpptraj::Structure::Disulfide disulfide;
    if (disulfide.InitDisulfide( argIn, Cpptraj::Structure::Disulfide::ADD_BONDS, debug_ )) {
      mprinterr("Error: Could not init disulfide search.\n");
      return 1;
    }
    if (disulfide.SearchForDisulfides( resStat, topIn, frameIn, DisulfideBonds ))
    {
      mprinterr("Error: Disulfide search failed.\n");
      return 1;
    }
  //} else {
  //  mprintf("\tNot searching for disulfides.\n");
  //}
  t_disulfide_.Stop();

  // Handle sugars.
  t_sugar_.Start();
  // TODO should be on a residue by residue basis in FillAtomsWithTemplates
  //bool prepare_sugars = !argIn.hasKey("nosugars");
  bool prepare_sugars = true; //FIXME
  if (!prepare_sugars)
    mprintf("\tNot attempting to prepare sugars.\n");
  else
    mprintf("\tWill attempt to prepare sugars.\n");
  Cpptraj::Structure::SugarBuilder sugarBuilder(debug_);
  if (prepare_sugars) {
    // Init options
    // FIXME
    bool hasglycam = false;
    double rescut = 8.0;
    double bondoffset = 0.2;
    std::string sugarmask("");
    std::string determinesugarsby("geometry");
    std::string resmapfile("");
    if (sugarBuilder.InitOptions( hasglycam,
                                  rescut,
                                  bondoffset,
                                  sugarmask,
                                  determinesugarsby,
                                  resmapfile ))
    {
      mprinterr("Error: Sugar options init failed.\n");
      return 1;
    }
    bool splitres = true; // FIXME
    //bool splitres = !argIn.hasKey("nosplitres");
    if (splitres)
      mprintf("\tWill split off recognized sugar functional groups into separate residues.\n");
    else
      mprintf("\tNot splitting recognized sugar functional groups into separate residues.\n");
    //bool c1bondsearch = !argIn.hasKey("noc1search");
    bool c1bondsearch = true; // FIXME
    if (c1bondsearch)
      mprintf("\tWill search for missing bonds to sugar anomeric atoms.\n");
    else
      mprintf("\tNot searching for missing bonds to sugar anomeric atoms.\n");
    // May need to modify sugar structure/topology, either by splitting
    // C1 hydroxyls of terminal sugars into ROH residues, and/or by
    // adding missing bonds to C1 atoms.
    // This is done before any identification takes place since we want
    // to identify based on the most up-to-date topology.
    if (sugarBuilder.FixSugarsStructure(topIn, frameIn,
                                        c1bondsearch, splitres, solventResName,
                                        SugarBonds))
    {
      mprinterr("Error: Sugar structure modification failed.\n");
      return 1;
    }
    if (sugarBuilder.PrepareSugars(true, resStat, topIn, frameIn, SugarBonds))
    {
      mprinterr("Error: Sugar preparation failed.\n");
      return 1;
    }
  }
  t_sugar_.Stop();

  // Fill in atoms with templates
  t_fill_.Start();
  //Topology topOut;
  Topology& topOut = static_cast<Topology&>( *(crdout.TopPtr()) );
  topOut.SetDebug( debug_ );
  if (outset.empty()) {
    // In-place output COORDS. Copy over existing topology metadata
    topOut.CopyTopMetadata( topIn );
  }
  if (!title.empty())
    topOut.SetParmName( title, FileName() );
  else if (topOut.ParmName().empty()) {
    // TODO better default
    topOut.SetParmName( std::string(topIn.c_str()), FileName() );
  }
  Frame frameOut;
  if (FillAtomsWithTemplates(topOut, frameOut, topIn, frameIn, creator, SugarBonds)) {
    mprinterr("Error: Could not fill in atoms using templates.\n");
    return 1;
  }
  t_fill_.Stop();

  // Add the disulfide/sugar bonds
  int addBondsErr = transfer_bonds( topOut, topIn, DisulfideBonds );
  //addBondsErr += transfer_bonds( topOut, topIn, SugarBonds );
  if (addBondsErr != 0) {
    mprinterr("Error: Adding disulfide/sugar bonds failed.\n");
    return 1;
  }

  // Create empty arrays for the TREE, JOIN, and IROTAT arrays
  topOut.AllocTreeChainClassification( );
  topOut.AllocJoinArray();
  topOut.AllocRotateArray();
    // Finalize topology - determine molecules, dont renumber residues
  topOut.CommonSetup(true, false);
  topOut.Summary();
  t_assign_.Stop();

  // Update coords 
  if (crdout.CoordsSetup( topOut, frameOut.CoordsInfo() )) { // FIXME better coordinate info
    mprinterr("Error: Could not set up output COORDS.\n");
    return CpptrajState::ERR;
  }
  crdout.SetCRD(0, frameOut);
  t_total_.Stop();

  PrintTiming();

  return 0;
}

void Exec_Build::PrintTiming() const {
  t_total_.WriteTiming(1, "Build timing:");
  t_hisDetect_.WriteTiming    (2, "Histidine detection :", t_total_.Total());
  t_clean_.WriteTiming        (2, "Structure clean     :", t_total_.Total());
  t_get_templates_.WriteTiming(2, "Get templates/parms :", t_total_.Total());
  
  t_disulfide_.WriteTiming    (2, "Disulfide detection :", t_total_.Total());
  t_sugar_.WriteTiming        (2, "Sugar preparation   :", t_total_.Total());
  t_fill_.WriteTiming         (2, "Fill missing atoms  :", t_total_.Total());
  t_fill_template_.WriteTiming(3, "Get template atoms :", t_fill_.Total());
  t_fill_build_.WriteTiming   (3, "Build atoms        :", t_fill_.Total());
  t_fill_build_internals_.WriteTiming(4, "Internals :", t_fill_build_.Total());
  t_fill_build_build_.WriteTiming    (4, "Build     :", t_fill_build_.Total());
  t_fill_build_link_.WriteTiming     (4, "Link      :", t_fill_build_.Total());
  t_fill_build_link_bond_.WriteTiming(5, "Link Bond :", t_fill_build_link_.Total());
  Cpptraj::Structure::Builder::PrintTiming(5, t_fill_build_link_.Total());
}
