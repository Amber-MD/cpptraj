#include "Traj_PDBfile.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_PDBfile::Traj_PDBfile() :
  radiiMode_(GB),
  terMode_(BY_MOL),
  pdbAtom_(0),
  currentSet_(0),
  ter_num_(0),
  pdbWriteMode_(NONE),
  dumpq_(false),
  pdbres_(false),
  pdbatom_(false),
  write_cryst1_(false),
  include_ep_(false),
  writeConect_(false),
  prependExt_(false),
  pdbTop_(0),
  chainchar_(' ')
{}

//------------------------------------------------------------------------
bool Traj_PDBfile::ID_TrajFormat(CpptrajFile& fileIn) {
  return PDBfile::ID_PDB( fileIn );
}

// Traj_PDBfile::closeTraj()
/** If not writing one PDB per frame write the END record. */
void Traj_PDBfile::closeTraj() {
  if ( (pdbWriteMode_ == SINGLE || pdbWriteMode_ == MODEL) &&
        file_.IsOpen() )
    file_.WriteEND();
  if (pdbWriteMode_ != MULTI)
    file_.CloseFile();
}

// Traj_PDBfile::openTrajin()
int Traj_PDBfile::openTrajin() {
  currentSet_ = 0;
  return file_.OpenFile();
}

// Traj_PDBfile::setupTrajin()
/** Scan PDB file to determine number of frames (models). The first frame will 
  * also be checked to ensure that the atom names match those in the parm file
  * in TrajectoryFile.
  */
int Traj_PDBfile::setupTrajin(FileName const& fname, Topology* trajParm)
{
  int atom;
  pdbWriteMode_ = NONE;
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  if (file_.OpenFile()) return TRAJIN_ERR;
  // Two strats - check for MODEL keywords or see how many times natom ATOM
  // records can be read. Currently employing the latter.
  int Frames = 0;
  int numMismatch = 0;
  bool scanPDB = true;
  Box boxInfo;
  while (scanPDB) {
    atom = 0;
    while ( atom < trajParm->Natom() ) {
      //fprintf(stdout,"DEBUG: PDB Read atom %i\n",atom);
      if ( file_.NextRecord() == PDBfile::END_OF_FILE ) {
        scanPDB = false;
        break;
      }
      //fprintf(stdout,"DEBUG: PDB buffer %i [%s]\n",atom,buffer);
      if (file_.RecType() ==  PDBfile::CRYST1) {
        // Read in box information
        double box_crd[6];
        file_.pdb_Box( box_crd );
        boxInfo.SetBox( box_crd );
      } 
      // Skip non-ATOM records
      if (file_.RecType() != PDBfile::ATOM) continue;
      // If still on first frame, check pdb atom name against the name in the 
      // associated parm file.
      if (Frames==0) {
        Atom pdbAtom = file_.pdb_Atom();
        if ( pdbAtom.Name() != (*trajParm)[atom].Name() ) {
          if (debug_>1) 
            mprintf("Warning: %s: PDB atom %i name [%s] does not match parm atom name [%s]\n",
                    file_.Filename().base(), atom+1, *(pdbAtom.Name()), 
                    *((*trajParm)[atom].Name()));
          ++numMismatch;
        }
      }
      ++atom;
    }
    if (Frames==0) {
      // First frame #atoms 
      pdbAtom_ = atom;
    } else {
      // Check that # atoms read in this frame match the first frame
      if (atom>0 && pdbAtom_!=atom) {
        mprintf("Warning: PDB %s: Reading frame %i, got %i atoms, expected %i.\n",
                file_.Filename().base(), Frames+1, atom, pdbAtom_);
        mprintf("Warning: Only using frames 1-%i\n", Frames);
        scanPDB = false;
        break;
      }
    }  
    if (scanPDB) ++Frames;
  }
  file_.CloseFile(); 
  if (Frames<1) {
    mprinterr("Error: PDB %s: No frames read. atom=%i expected %i.\n",
              file_.Filename().base(), atom, trajParm->Natom());
    return TRAJIN_ERR;
  }
  if (debug_>0) mprintf("Traj_PDBfile: %s has %i atoms, %i frames.\n",
                        file_.Filename().base(), pdbAtom_, Frames);
  // Report mismatches of pdb atom names against parm names
  if (numMismatch > 0)
    mprintf("Warning: In PDB file %s: %i name mismatches with parm %s.\n",
            file_.Filename().base(), numMismatch, trajParm->c_str());
  // Set traj info - no velocity, temperature, time
  SetCoordInfo( CoordinateInfo(boxInfo, false, false, false) );
  return Frames;
}

// Traj_PDBfile::readFrame()
/** Read frame (model) from PDB file. */
int Traj_PDBfile::readFrame(int set, Frame& frameIn)
{
  int atom;
  if (set < currentSet_) {
    file_.Rewind();
    currentSet_ = 0;
  }
  // Position file at group of ATOM keywords for specified set
  while (currentSet_ < set) {
    atom = 0;
    while (atom < pdbAtom_) {
      if ( file_.NextRecord() == PDBfile::END_OF_FILE ) return 1;
      if ( file_.RecType() == PDBfile::ATOM ) ++atom;
    }
    currentSet_++;
  }
  atom = 0;
  double *Xptr = frameIn.xAddress(); 
  while (atom < pdbAtom_) {
    if ( file_.NextRecord() == PDBfile::END_OF_FILE ) return 1;
    // Skip non-ATOM records
    if ( file_.RecType() == PDBfile::CRYST1 )
      file_.pdb_Box( frameIn.bAddress() );
    else if ( file_.RecType() == PDBfile::ATOM ) {
      // Read current PDB record XYZ into Frame
      file_.pdb_XYZ( Xptr );
      ++atom; 
      Xptr += 3;
    }
  }
  currentSet_++;

  return 0;
}

void Traj_PDBfile::WriteHelp() {
  mprintf("\tdumpq      : Write atom charge/GB radius in occupancy/B-factor columns (PQR format).\n"
          "\tparse      : Write atom charge/PARSE radius in occupancy/B-factor columns (PQR format).\n"
          "\tvdw        : Write atom charge/VDW radius in occupancy/B-factor columns (PQR format).\n"
          "\tpdbres     : Use PDB V3 residue names.\n"
          "\tpdbatom    : Use PDB V3 atom names.\n"
          "\tpdbv3      : Use PDB V3 residue/atom names.\n"
          "\tteradvance : Increment record (atom) # for TER records (default no).\n"
          "\tterbyres   : Print TER cards based on residue sequence instead of molecules.\n"
          "\tpdbter     : Print TER cards according to original PDB TER (if available).\n"
          "\tnoter      : Do not write TER cards.\n"
          "\tmodel      : Write to single file separated by MODEL records.\n"
          "\tmulti      : Write each frame to separate files.\n"
          "\tchainid <c>: Write character 'c' in chain ID column.\n"
          "\tsg <group> : Space group for CRYST1 record, only if box coordinates written.\n"
          "\tinclude_ep : Include extra points.\n"
          "\tconect     : Write CONECT records using bond information.\n"
          "\tkeepext    : Keep filename extension; write '<name>.<num>.<ext>' instead (implies 'multi').\n");
}

// Traj_PDBfile::processWriteArgs()
int Traj_PDBfile::processWriteArgs(ArgList& argIn) {
  pdbWriteMode_ = SINGLE;
  if (argIn.hasKey("dumpq")) {
    dumpq_ = true; 
    radiiMode_ = GB;
  } else if (argIn.hasKey("parse")) {
    dumpq_ = true;
    radiiMode_ = PARSE;
  } else if (argIn.hasKey("vdw") || argIn.hasKey("dumpr*")) {
    dumpq_ = true;
    radiiMode_ = VDW;
  }
  if (argIn.hasKey("terbyres"))   terMode_ = BY_RES;
  else if (argIn.hasKey("noter")) terMode_ = NO_TER;
  else if (argIn.hasKey("pdbter"))terMode_ = ORIGINAL_PDB;
  else                            terMode_ = BY_MOL;
  pdbres_ = argIn.hasKey("pdbres");
  pdbatom_ = argIn.hasKey("pdbatom");
  if (argIn.hasKey("pdbv3")) {
    pdbres_ = true;
    pdbatom_ = true;
  }
  if (argIn.hasKey("teradvance"))
    ter_num_ = 1;
  else
    ter_num_ = 0;
  if (argIn.hasKey("model")) pdbWriteMode_ = MODEL;
  if (argIn.hasKey("multi")) pdbWriteMode_ = MULTI;
  include_ep_ = argIn.hasKey("include_ep");
  writeConect_ = argIn.hasKey("conect");
  prependExt_ = argIn.hasKey("keepext"); // Implies MULTI
  if (prependExt_) pdbWriteMode_ = MULTI;
  space_group_ = argIn.GetStringKey("sg");
  std::string temp = argIn.GetStringKey("chainid");
  if (!temp.empty()) chainchar_ = temp[0];
  return 0;
}

// Traj_PDBfile::setupTrajout()
/** Set parm information needed for write, and check write mode against
  * number of frames to be written.
  */ 
int Traj_PDBfile::setupTrajout(FileName const& fname, Topology* trajParm,
                               CoordinateInfo const& cInfoIn,
                               int NframesToWrite, bool append)
{
  if (trajParm==0) return 1;
  SetCoordInfo( cInfoIn );
  pdbTop_ = trajParm;
  pdbAtom_ = pdbTop_->Natom();
  // Set up file
  if (append && pdbWriteMode_ != MULTI) {
    if (file_.SetupAppend( fname, debug_)) return 1;
  } else {
    if (append && pdbWriteMode_ == MULTI)
      mprintf("Warning: 'append' not compatible with 'multi' pdb write.\n");
    if (file_.SetupWrite( fname, debug_ )) return 1;
  }
  // Set a chainID for each residue 
  // TODO: Set different chain ID for solute mols and solvent
  chainID_.clear();
  if (chainchar_ == ' ') {
    chainID_.reserve( trajParm->Nres() );
    for (Topology::res_iterator res = trajParm->ResStart(); res != trajParm->ResEnd(); ++res)
      chainID_.push_back( res->ChainID() );
  } else
    chainID_.resize(trajParm->Nres(), chainchar_);
        
  // Save residue names. If pdbres specified convert to PDBV3 residue names.
  resNames_.clear();
  resNames_.reserve( trajParm->Nres() );
  resIsHet_.clear();
  resIsHet_.reserve( trajParm->Nres() );
  if (pdbres_) {
    for (Topology::res_iterator res = trajParm->ResStart();
                                res != trajParm->ResEnd(); ++res)
    {
      NameType rname = res->Name();
      // First check if this is water.
      if ( res->NameIsSolvent() )
        rname = "HOH ";
      // convert protein residue names back to more like PDBV3 format:
      else if (rname == "HID " || rname == "HIE " ||
               rname == "HIP " || rname == "HIC "   )
        rname = "HIS ";
      else if (rname == "CYX " || rname == "CYM " || rname == "CYZ ")
        rname = "CYS ";
      else if (rname == "MEM ") 
        rname = "MET ";
      else if (rname == "ASH ")
        rname = "ASP ";
      else if (rname == "GLH ")
        rname = "GLU ";
      // also for nucleic acid names:
      else if ( rname == "C3  " )  rname = "  C ";
      else if ( rname == "U3  " )  rname = "  U ";
      else if ( rname == "G3  " )  rname = "  G ";
      else if ( rname == "A3  " )  rname = "  A ";
      else if ( rname == "C5  " )  rname = "  C ";
      else if ( rname == "U5  " )  rname = "  U ";
      else if ( rname == "G5  " )  rname = "  G ";
      else if ( rname == "A5  " )  rname = "  A ";
      else if ( rname == "DC3 " )  rname = " DC ";
      else if ( rname == "DT3 " )  rname = " DT ";
      else if ( rname == "DG3 " )  rname = " DG ";
      else if ( rname == "DA3 " )  rname = " DA ";
      else if ( rname == "DC5 " )  rname = " DC ";
      else if ( rname == "DT5 " )  rname = " DT ";
      else if ( rname == "DG5 " )  rname = " DG ";
      else if ( rname == "DA5 " )  rname = " DA ";
      else if ( rname == "URA " || rname == "URI" )
        rname = "  U ";
      else if ( rname == "THY " )
        rname = " DT ";
      else if ( rname == "GUA" || rname == "ADE" || rname == "CYT" ) {
        // Determine if RNA or DNA via existence of O2'
        bool isRNA = false;
        for (int ratom = res->FirstAtom(); ratom != res->LastAtom(); ++ratom)
          if ( (*trajParm)[ratom].Name() == "O2'" ||
               (*trajParm)[ratom].Name() == "O2*" )
          {
            isRNA = true;
            break;
          }
        if (isRNA) {
          if      (rname[0] == 'G') rname="  G ";
          else if (rname[0] == 'A') rname="  A ";
          else if (rname[0] == 'C') rname="  C ";
        } else {
          if      (rname[0] == 'G') rname=" DG ";
          else if (rname[0] == 'A') rname=" DA ";
          else if (rname[0] == 'C') rname=" DC ";
        }
      }
      resNames_.push_back( rname );
      // Any non-standard residue should get HETATM
      if ( rname == "ALA " ||
           rname == "ARG " ||
           rname == "ASN " ||
           rname == "ASP " ||
           rname == "ASX " ||
           rname == "CYS " ||
           rname == "GLN " ||
           rname == "GLU " ||
           rname == "GLX " ||
           rname == "GLY " ||
           rname == "HIS " ||
           rname == "ILE " ||
           rname == "LEU " ||
           rname == "LYS " ||
           rname == "MET " ||
           rname == "PHE " ||
           rname == "PRO " ||
           rname == "SER " ||
           rname == "THR " ||
           rname == "TRP " ||
           rname == "TYR " ||
           rname == "UNK " ||
           rname == "VAL " ||
           rname == "  C " ||
           rname == "  G " ||
           rname == "  A " ||
           rname == "  U " ||
           rname == "  I " ||
           rname == " DC " ||
           rname == " DG " ||
           rname == " DA " ||
           rname == " DU " ||
           rname == " DT " ||
           rname == " DI "    )
        resIsHet_.push_back( false );
      else
        resIsHet_.push_back( true );
      //mprintf("DEBUG: ResName='%s' IsHet=%i\n", *(resNames_.back()), (int)resIsHet_.back());
    }
  } else {
    for (Topology::res_iterator res = trajParm->ResStart();
                                res != trajParm->ResEnd(); ++res)
      resNames_.push_back( res->Name() );
    resIsHet_.assign( trajParm->Nres(), false );
  }
  // Set up TER cards.
  TER_idxs_.clear();
  if (terMode_ == ORIGINAL_PDB) {
    // Write a TER card only after residues that originally have a TER card
    // or when chain ID changes.
    Topology::res_iterator res = trajParm->ResStart();
    for (; res != trajParm->ResEnd(); ++res) {
      if (res->IsTerminal() || res+1 == trajParm->ResEnd())
        TER_idxs_.push_back( res->LastAtom() - 1 );
      else if ((res+1)->ChainID() != res->ChainID())
        TER_idxs_.push_back( res->LastAtom() - 1 );
    }
  } else if (terMode_ == BY_RES) {
    bool lastResWasSolvent = false;
    // Write a TER card every time residue of atom N+1 is not bonded to any
    // atom of residue of atom N. Do not do this for solvent.
    for (Topology::res_iterator res = trajParm->ResStart();
                                res != trajParm->ResEnd(); ++res)
    {
      bool isIon = false;
      if (trajParm->Nmol() > 0) {
        int molNum = (*trajParm)[ res->FirstAtom() ].MolNum();
        // If this is a one atom molecule assume it is an ion.
        isIon = trajParm->Mol( molNum ).NumAtoms() == 1; 
      }
      if (!res->NameIsSolvent() && !isIon) {
        // If this is the last residue, terminate the chain with final atom.
        // FIXME build this into the loop.
        if ( res+1 == trajParm->ResEnd() )
          TER_idxs_.push_back( res->LastAtom() - 1 );
        else if ( lastResWasSolvent )
          TER_idxs_.push_back( (res-1)->LastAtom() - 1 );
        else {
          int r2_first = (res+1)->FirstAtom();
          int r2_last  = (res+1)->LastAtom();
          bool residues_are_bonded = false;
          for (int r1_at = res->LastAtom()-1; r1_at >= res->FirstAtom(); r1_at--)
          {
            for (Atom::bond_iterator bnd_at = (*trajParm)[r1_at].bondbegin();
                                     bnd_at != (*trajParm)[r1_at].bondend(); ++bnd_at)
            {
              if ( *bnd_at >= r2_first && *bnd_at < r2_last ) {
                residues_are_bonded = true;
                break;
              }
            }
            if (residues_are_bonded) break;
          }
          if (!residues_are_bonded)
            TER_idxs_.push_back( res->LastAtom() - 1 );
        }
        lastResWasSolvent = false;
      } else
        lastResWasSolvent = true;
    }
  } else if (terMode_ == BY_MOL) {
    // Write a TER card at the end of every molecule
    // NOTE: For backwards compat. dont do this for last mol. FIXME Is this ok?
    if ( trajParm->Nmol() > 0 ) {
      Topology::mol_iterator finalMol = trajParm->MolEnd() - 1;
      for (Topology::mol_iterator mol = trajParm->MolStart();
                                  mol != finalMol; ++mol)
        TER_idxs_.push_back( mol->EndAtom() - 1 );
    }
  }
  TER_idxs_.push_back( -1 ); // Indicates that final TER has been written.
  if (debug_ > 0) {
    mprintf("DEBUG: TER indices:");
    for (std::vector<int>::const_iterator idx = TER_idxs_.begin(); idx != TER_idxs_.end(); ++idx)
      mprintf(" %i", *idx + 1);
    mprintf("\n");
  }
  // Allocate space to hold ATOM record #s if writing CONECT records
  if (writeConect_) atrec_.resize( trajParm->Natom() );
  // If number of frames to write > 1 and not doing 1 pdb file per frame,
  // set write mode to MODEL
  if (append || (pdbWriteMode_==SINGLE && NframesToWrite>1)) 
    pdbWriteMode_ = MODEL;
  // TODO: Setup title
  // Open here if writing to single file
  if (pdbWriteMode_ != MULTI) {
    if ( file_.OpenFile() ) return 1;
    if (!Title().empty()) file_.WriteTITLE( Title() );
  }
  write_cryst1_ = (CoordInfo().TrajBox().Type() != Box::NOBOX);
  if (write_cryst1_) {
    if (pdbWriteMode_==MODEL)
      mprintf("Warning: For PDB with MODEL, box coords for first frame only will be written to CRYST1.\n");
    if (space_group_.empty())
      mprintf("Warning: No PDB space group specified.\n");
  }
  // Set up radii
  if (dumpq_) {
    radii_.clear();
    for (int iat = 0; iat != trajParm->Natom(); iat++) {
      switch (radiiMode_) {
        case GB:    radii_.push_back( (*trajParm)[iat].GBRadius() ); break;
        case PARSE: radii_.push_back( (*trajParm)[iat].ParseRadius() ); break;
        case VDW:   radii_.push_back( trajParm->GetVDWradius(iat) ); break;
      }
    }
  }
  return 0;
}

// Traj_PDBfile::writeFrame()
/** Write the frame (model) to PDB file. */
int Traj_PDBfile::writeFrame(int set, Frame const& frameOut) {
  if (pdbWriteMode_==MULTI) {
    // If writing 1 pdb per frame set up output filename and open
    if (file_.OpenWriteNumbered( set + 1, prependExt_ )) return 1;
    if (!Title().empty()) 
      file_.WriteTITLE( Title() );
    if (write_cryst1_)
      file_.WriteCRYST1( frameOut.BoxCrd().boxPtr(), space_group_.c_str() );
  } else {
    // Write box coords, first frame only.
    if (write_cryst1_) {
      file_.WriteCRYST1( frameOut.BoxCrd().boxPtr(), space_group_.c_str() );
      write_cryst1_ = false;
    }
  }
  // If specified, write MODEL keyword
  if (pdbWriteMode_==MODEL)
    file_.WriteMODEL(set + 1); 

  float Occ = 1.0; 
  float B = 0.0;
  char altLoc = ' ';
  int anum = 1; // Actual PDB ATOM record number
  const double *Xptr = frameOut.xAddress();
  std::vector<int>::const_iterator terIdx = TER_idxs_.begin();
  for (int aidx = 0; aidx != pdbTop_->Natom(); aidx++, Xptr += 3) {
    Atom const& atom = (*pdbTop_)[aidx];
    int res = atom.ResNum();
    if (include_ep_ || atom.Element() != Atom::EXTRAPT) {
      PDBfile::PDB_RECTYPE rectype;
      if ( resIsHet_[res] )
        rectype = PDBfile::HETATM;
      else
        rectype = PDBfile::ATOM;
      if (!pdbTop_->Extra().empty()) {
        Occ = pdbTop_->Extra()[aidx].Occupancy();
        B   = pdbTop_->Extra()[aidx].Bfactor();
        altLoc = pdbTop_->Extra()[aidx].AtomAltLoc();
      }
      if (dumpq_) {
        Occ = (float) atom.Charge();
        B = (float) radii_[aidx];
      }
      // If pdbatom change amber atom names to pdb v3
      NameType atomName = atom.Name();
      if (pdbatom_) {
        if      (atomName == "H5'1") atomName = "H5'";
        else if (atomName == "H5'2") atomName = "H5''";
        else if (atomName == "H2'1") atomName = "H2'";
        else if (atomName == "H2'2") atomName = "H2''";
        else if (atomName == "O1P ") atomName = "OP1";
        else if (atomName == "O2P ") atomName = "OP2";
        else if (atomName == "H5T ") atomName = "HO5'";
        else if (atomName == "H3T ") atomName = "HO3'";
        else if (atomName == "HO'2") atomName = "HO2'";
      }
      file_.WriteCoord(rectype, anum, atomName, altLoc, resNames_[res],
                       chainID_[res], pdbTop_->Res(res).OriginalResNum(),
                       pdbTop_->Res(res).Icode(),
                       Xptr[0], Xptr[1], Xptr[2], Occ, B,
                       atom.ElementName(), 0, dumpq_);
      if (writeConect_) atrec_[aidx] = anum; // Store ATOM record #
    }
    anum++;
    // Check and see if a TER card should be written.
    if (aidx == *terIdx) {
      // FIXME: Should anum not be incremented until after? 
      file_.WriteRecordHeader(PDBfile::TER, anum, "", ' ', resNames_[res],
                              chainID_[res], pdbTop_->Res(res).OriginalResNum(),
                              pdbTop_->Res(res).Icode(), atom.ElementName());
      anum += ter_num_;
      ++terIdx;
    }
  }
  // Write CONECT records for each ATOM
  if (writeConect_) {
    for (int aidx = 0; aidx != pdbTop_->Natom(); aidx++)
      file_.WriteCONECT( atrec_[aidx], atrec_, (*pdbTop_)[aidx] );
  }
  if (pdbWriteMode_==MULTI) {
    // If writing 1 pdb per frame, close output file
    file_.WriteEND();
    file_.CloseFile();
  } else if (pdbWriteMode_==MODEL) {
    // If MODEL keyword was written, write corresponding ENDMDL record
    file_.WriteENDMDL();
  }

  return 0;
}

// Traj_PDBfile::Info()
void Traj_PDBfile::Info() {
  mprintf("is a PDB file");
  if (pdbWriteMode_ != NONE) {
    if (pdbWriteMode_==MULTI)
      mprintf(" (1 file per frame)");
    else if (pdbWriteMode_==MODEL)
      mprintf(" (1 MODEL per frame)");
    if (writeConect_) mprintf(" with CONECT records");
    if (dumpq_) {
      mprintf(", writing charges to occupancy column and ");
      switch (radiiMode_) {
        case GB: mprintf("GB radii"); break;
        case PARSE: mprintf("PARSE radii"); break;
        case VDW: mprintf("vdW radii"); break;
      }
      mprintf(" to B-factor column");
    }
    if (pdbres_ && pdbatom_)
      mprintf(", using PDB V3 res/atom names");
    else if (pdbres_)
      mprintf(", using PDB V3 residue names");
    else if (pdbatom_)
      mprintf(", using PDB V3 atom names");
  }
}
#ifdef MPI
/// Not valid for MODEL (checked in setup) so no need to do anything.
int Traj_PDBfile::parallelOpenTrajout(Parallel::Comm const& commIn) { return 0; }

int Traj_PDBfile::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{
  if (pdbWriteMode_ != MULTI) {
    mprinterr("Error: PDB write in parallel requires 'multi' keyword.\n");
    return 1;
  }
  return setupTrajout(fname, trajParm, cInfoIn, NframesToWrite, append);
}

int Traj_PDBfile::parallelWriteFrame(int set, Frame const& frameOut) {
  return ( writeFrame(set, frameOut) );
}
#endif
