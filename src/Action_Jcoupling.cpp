// Action_Jcoupling
#include <cstdlib> //getenv
#include <cstdio> //sscanf
#include <cstring> //strcpy, strlen
#include <cmath> //cos
#include "Action_Jcoupling.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD, RADDEG
#include "TorsionRoutines.h"
#include "BufferedLine.h"

// CONSTRUCTOR
Action_Jcoupling::Action_Jcoupling() :
  debug_(0),
  Nconstants_(0),
  CurrentParm_(0),
  outputfile_(0),
  masterDSL_(0),
  setcount_(1)
{} 

void Action_Jcoupling::Help() const {
  mprintf("\t<mask1> [outfile <filename>] [kfile <param file>] [out <filename>]\n"
          "\t[name <dsname>]\n"
          "  Calculate J-coupling values for all dihedrals found in <mask1>.\n");
}

static inline int fileEOF(std::string const& filename) {
  mprinterr("Error: Unexpected EOF when reading Jcoupling values from file '%s'\n",
            filename.c_str());
  return 1;
}

// Action_Jcoupling::loadKarplus()
/** Load Karplus parameters from input file.
  * Expected format:
  * - {type}<+|-| ><a[4]><+|-| ><b[4]><+|-| ><c[4]><+|-| ><d[4]><A[6]><B[6]><C[6]>{<D[6]>}
  *   <reslabel[4]>* 
  * \return 0 on success, 1 on error
  */
int Action_Jcoupling::loadKarplus(std::string const& filename) {
  BufferedLine KarplusFile;

  if (filename.empty()) {
    mprinterr("Error: jcoupling: Could not find Karplus parameter file.\n");
    return 1;
  }
  if (KarplusFile.OpenFileRead( filename )) {
    mprinterr("Error: jcoupling: Could not read Karplus parameter file %s\n",
              filename.c_str());
    mprinterr("Error: Ensure the file exists and is readable.\n");
    return 1;
  }
  // Read through all lines of the file
  const char* buffer = KarplusFile.Line();
  if (buffer == 0) return fileEOF(filename);
  while (buffer != 0) {
    // Skip empty lines (BufferedLine.Line() removes newlines) and comments
    while (buffer[0]=='\0' || buffer[0]=='#')
      buffer = KarplusFile.Line();
    if (buffer == 0) return fileEOF(filename);
    const char* ptr = buffer;
    // First char is optional type. If optional type is C, then the Karplus 
    // function specified in Perez et al. JACS (2001) 123 will be used, and 
    // A, B, and C will be taken as C0, C1, and C2.
    jcoupleDihedral KC;
    if (ptr[0] == 'C') {
      KC.SetType( PEREZ );
      ptr++;
    }
    // Read 4 atom names with optional preceding character (+, -)
    for (int i = 0; i < 4; i++) {
      // Check for offset char
      if (*ptr == '+')
        KC.SetOffset(i, 1);
      else if (*ptr == '-')
        KC.SetOffset(i, -1);
      ++ptr;
      char aname[5];
      aname[0] = ptr[0];
      aname[1] = ptr[1];
      aname[2] = ptr[2];
      aname[3] = ptr[3];
      aname[4] = '\0';
      KC.SetName(i, aname);
      if (debug_ > 1)
        mprintf("DEBUG:\tAtomName %i [%s] '%s' %c\n", i, *(KC.AtomName(i)), aname, *ptr);
      ptr += 4;
    }
    // Read parameters
    // NOTE: Using sscanf here instead of atof since the 4th parameter is
    //       optional, behavior is undefined for accessing uninitialized
    //       portion of buffer.
    double* cptr = KC.Cptr();
    int nread = sscanf(ptr, "%6lf%6lf%6lf%6lf", cptr, cptr+1, cptr+2, cptr+3);
    if (nread < 3) {
      mprinterr("Error: jcoupling: Expected at least 3 Karplus parameters, got %i\n", nread);
      mprinterr("Error: Line: [%s]\n", buffer);
      return 1;
    }
    // C3 needs to be in radians since the dihedral will be calculated in radians.
    KC.C(3) *= Constants::DEGRAD;
    /// The next line contains all residues
    buffer = KarplusFile.Line();
    if (buffer == 0) {
      mprinterr("Error: EOF encountered in '%s' before residue line could be read.\n",
                filename.c_str());
      return 1;
    }
    // Read the list of residues to which this constant will apply.
    const char* end = buffer + strlen(buffer);
    for (ptr = buffer; ptr < end; ptr += 4) {
      char rname[5];
      rname[0] = ptr[0];
      rname[1] = ptr[1];
      rname[2] = ptr[2];
      rname[3] = ptr[3];
      rname[4] = '\0';
      NameType resName( rname );
      if (debug_ > 1)
        mprintf("DEBUG:\t Residue [%s]\n", *resName);
      JcoupleDihedralMap::iterator it = JcoupleData_.lower_bound( resName );
      if (it == JcoupleData_.end() || it->first != resName ) {
        // List does not exist for residue yet, create it.
        if (debug_ > 1)
          mprintf("DEBUG: Creating new list for residue.\n");
        JcoupleData_.insert( it,
                             JcoupleDihedralPair(resName,
                                                 JcoupleDihedralArray(1, KC)
                                                )
                           );
      } else {
        // Add constants to this residues list.
        if (debug_ > 1)
          mprintf("DEBUG: Adding constant to list for residue.\n");
        it->second.push_back( KC );
      }
      ++Nconstants_;
    } // END loop over residues in residue line
    buffer = KarplusFile.Line();
  } // END Gets over input file
  KarplusFile.CloseFile();
  // DEBUG - Print out all parameters
  if (debug_>0) {
      mprintf("    KARPLUS PARAMETERS:\n");
      for (JcoupleDihedralMap::const_iterator it = JcoupleData_.begin();
                                              it != JcoupleData_.end(); ++it)
      {
        mprintf("\t[%4s]\n", *(it->first));
        for (JcoupleDihedralArray::const_iterator kc = it->second.begin();
                                                  kc != it->second.end(); ++kc)
        {
          mprintf("\t\t%1i ", (int)kc->Type());
          for (int i = 0; i < 4; i++)
            mprintf(" %4s", *(kc->AtomName(i)));
          for (int i = 0; i < 4; i++)
            mprintf(" %1i", kc->Offset(i));
          for (int i = 0; i < 4; i++)
            mprintf(" %6.2f", kc->Constant(i));
          mprintf("\n");
        }
      }
  }
  return 0;
}

// -----------------------------------------------------------------------------
// Action_Jcoupling::Init()
Action::RetType Action_Jcoupling::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  outfile_ = 0;
  // Get Keywords
  outputfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("outfile"), "J-coupling");
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  std::string karpluspath = actionArgs.GetStringKey("kfile");
  setname_ = actionArgs.GetStringKey("name");
  // Get Masks
  if (Mask1_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;

  // If no Karplus params specified check environment vars. 
  if (karpluspath.empty()) {
    // Check if the KARPLUS env var is set.
    const char* env = getenv("KARPLUS");
    if (env != 0) {
      //mprintf("Info: Using parameter file defined by $KARPLUS environment variable.\n");
      karpluspath.assign(env);
      mprintf("Info: Parameter file path from KARPLUS variable: '%s'\n", karpluspath.c_str());
    } 
    // Check if AMBERHOME is set.
    if (karpluspath.empty()) {
      env = getenv("AMBERHOME");
      if (env != 0) {
        karpluspath.assign(env);
        karpluspath += "/dat/Karplus.txt";
      }
      mprintf("Info: Parameter file path from AMBERHOME variable: '%s'\n", karpluspath.c_str());
    }
    // Last, use CPPTRAJHOME
    if (karpluspath.empty()) {
      env = getenv("CPPTRAJHOME");
      if (env != 0) {
        karpluspath.assign(env);
        karpluspath += "/dat/Karplus.txt";
      }
      mprintf("Info: Parameter file path from CPPTRAJHOME variable: '%s'\n", karpluspath.c_str());
    }
    // If no path, bail out
    if (karpluspath.empty()) {
      mprinterr("Error: Either the CPPTRAJHOME or AMBERHOME env. variables must be set, or the\n"
                "Error:   parameter file location must be specified by 'kfile' or the KARPLUS\n"
                "Error:   env. variable.\n");
      return Action::ERR;
    }
  }
  // Load Karplus parameters
  if (loadKarplus(karpluspath)) 
    return Action::ERR;

  mprintf("    J-COUPLING: Searching for dihedrals in mask [%s].\n"
          "\tUsing Karplus parameters in \"%s\"\n"
          "\t%i parameters found for %zu residues.\n",
          Mask1_.MaskString(), karpluspath.c_str(), Nconstants_, JcoupleData_.size());
  if (outfile_ != 0)
    mprintf("\tDataSets will be written to %s\n", outfile_->DataFilename().full());
  if (outputfile_ != 0)
    mprintf("\tWriting fixed-format output to %s\n",outputfile_->Filename().full());
  mprintf("# Citations: Chou et al. JACS (2003) 125 p.8959-8966\n"
          "#            Perez et al. JACS (2001) 123 p.7081-7093\n");
  init.DSL().SetDataSetsPending(true);
  masterDSL_ = init.DslPtr();
  return Action::OK;
}

// Action_Jcoupling::Setup()
/** Set up a j-coupling calculation for dihedrals defined by atoms within
  * the mask.
  */
Action::RetType Action_Jcoupling::Setup(ActionSetup& setup) {
  if ( setup.Top().SetupCharMask(Mask1_) ) return Action::ERR;
  if (Mask1_.None()) {
    mprintf("Warning: Mask specifies no atoms.\n");
    return Action::SKIP;
  }
  // If JcouplingInfo has already been set up, print a warning and reset for
  // new parm.
  if (!JcouplingInfo_.empty()) {
    mprintf("Warning: Jcoupling has been set up for another parm.\n"
            "Warning:   Resetting jcoupling info for new parm %s\n",setup.Top().c_str());
    JcouplingInfo_.clear();
  }

  // For each residue, set up 1 jcoupling calc for each parameter defined in
  // JcoupleData for this residue. Only set up the Jcoupling calc if all
  // atoms involved are present in the mask.
  Range resRange = setup.Top().SoluteResidues();
  for (Range::const_iterator rnum = resRange.begin();
                             rnum != resRange.end(); ++rnum)
  {
    // Skip residue if no atoms within residue are selected.
    if (!Mask1_.AtomsInCharMask(setup.Top().Res(*rnum).FirstAtom(),
                                setup.Top().Res(*rnum).LastAtom())) continue;
    // Try to find list of constants for this residue.
    JcoupleDihedralMap::const_iterator it = JcoupleData_.find( setup.Top().Res(*rnum).Name() );
    // If list does not exist for residue, skip it.
    if (it == JcoupleData_.end() ) {
      mprintf("Warning: J-coupling parameters not found for %s\n",
              setup.Top().TruncResNameNum(*rnum).c_str());
      continue;
    }
    JcoupleDihedralArray const& resConstants = it->second;
    // For each parameter set in the list find the corresponding atoms.
    for (JcoupleDihedralArray::const_iterator kc = resConstants.begin();
                                              kc != resConstants.end(); ++kc)
    {
      // Init jcoupling info. Constants will point inside KarplusConstants.
      jcouplingInfo JC;
      JC.residue = *rnum;
      JC.atom[0] = -1;
      JC.atom[1] = -1;
      JC.atom[2] = -1;
      JC.atom[3] = -1;
      JC.C = kc->Carray();
      JC.type = kc->Type();
      // For each atom in the dihedral specified in this Karplus constant, find
      // corresponding atoms in parm.
      bool allAtomsFound = true;
      for (int idx=0; idx < 4; idx++) {
        JC.atom[idx] = setup.Top().FindAtomInResidue(*rnum + kc->Offset(idx),
                                                      kc->AtomName(idx)       );
        if (JC.atom[idx] == -1) {
          mprintf("Warning: Atom '%s' at position %i not found for residue %i\n",
                    *(kc->AtomName(idx)), idx, *rnum+kc->Offset(idx)+1);
          allAtomsFound = false;
        }
      }
      if (allAtomsFound) { 
        // Check that all the atoms involved in this Jcouple dihedral are
        // in the atom mask. If so, add jcoupling info to the list.
        if (Mask1_.AtomInCharMask(JC.atom[0]) && Mask1_.AtomInCharMask(JC.atom[1]) &&
            Mask1_.AtomInCharMask(JC.atom[2]) && Mask1_.AtomInCharMask(JC.atom[3]))
        {
          // TODO: Look for previously set up matching data set
          if (setname_.empty())
            setname_ = masterDSL_->GenerateDefaultName("JC");
          JC.data_ = masterDSL_->AddSet( DataSet::FLOAT, MetaData(setname_, setcount_++) );
          if ( JC.data_ != 0 ) {
            JC.data_->SetLegend( setup.Top().TruncResNameNum(JC.residue) + "_" +
                                 setup.Top()[JC.atom[0]].Name().Truncated() + "-" +
                                 setup.Top()[JC.atom[1]].Name().Truncated() + "-" +
                                 setup.Top()[JC.atom[2]].Name().Truncated() + "-" +
                                 setup.Top()[JC.atom[3]].Name().Truncated()  );
            if (outfile_ != 0)
              outfile_->AddDataSet( JC.data_ ); 
            JcouplingInfo_.push_back(JC);
          } else {
            mprinterr("Internal Error: Could not set up Jcoupling data set for res %i\n",
                      JC.residue+1);
          }
        }
      }
    } // END loop over karplus parameters for this residue
  } // END loop over all residues

  // Print info for this parm
  mprintf("    J-COUPLING: [%s] Will calculate J-coupling for %zu dihedrals.\n",
          Mask1_.MaskString(), JcouplingInfo_.size());
  if (JcouplingInfo_.empty()) {
    mprintf("Warning: No dihedrals found for J-coupling calculation!\n"
            "Warning:   Check that all atoms of dihedrals are included in mask [%s]\n"
            "Warning:   and/or that dihedrals are defined in Karplus parameter file.\n",
            Mask1_.MaskString());
    return Action::SKIP;
  }
  // DEBUG
  if (debug_>0) {
    int MaxResidues=1;
    for (std::vector<jcouplingInfo>::iterator jc = JcouplingInfo_.begin();
                                              jc != JcouplingInfo_.end(); ++jc) 
    {
      mprintf("%8i [%i:%4s]",MaxResidues,jc->residue, setup.Top().Res(jc->residue).c_str());
      mprintf(" %6i:%-4s",jc->atom[0],setup.Top()[jc->atom[0]].c_str());
      mprintf(" %6i:%-4s",jc->atom[1],setup.Top()[jc->atom[1]].c_str());
      mprintf(" %6i:%-4s",jc->atom[2],setup.Top()[jc->atom[2]].c_str());
      mprintf(" %6i:%-4s",jc->atom[3],setup.Top()[jc->atom[3]].c_str());
      mprintf(" %6.2lf%6.2lf%6.2lf%6.2lf %i\n",jc->C[0],jc->C[1],jc->C[2],jc->C[3],
              jc->type);
      MaxResidues++;
    }
  }
  CurrentParm_ = setup.TopAddress();
  return Action::OK;  
}

// Action_Jcoupling::DoAction()
/** For each dihedral defined in JcouplingInfo, perform the dihedral and
  * Jcoupling calculation.
  */
Action::RetType Action_Jcoupling::DoAction(int frameNum, ActionFrame& frm) {
  double Jval;

  if (outputfile_ != 0)
    outputfile_->Printf("#Frame %i\n",frameNum+1);

  for (std::vector<jcouplingInfo>::iterator jc = JcouplingInfo_.begin();
                                            jc !=JcouplingInfo_.end(); ++jc)
  {
    double phi = Torsion(frm.Frm().XYZ(jc->atom[0]),
                         frm.Frm().XYZ(jc->atom[1]),
                         frm.Frm().XYZ(jc->atom[2]),
                         frm.Frm().XYZ(jc->atom[3]) );
    if (jc->type==1) {
      //phitemp = phi + jc->C[3]; // Only necessary if offsets become used in perez-type calc
      Jval = jc->C[0] + (jc->C[1] * cos(phi)) + (jc->C[2] * cos(phi * 2.0)); 
    } else {
      double phitemp = cos( phi + jc->C[3] );
      Jval = (jc->C[0] * phitemp * phitemp) + (jc->C[1] * phitemp) + jc->C[2];
    }
    float fval = (float)Jval;
    jc->data_->Add(frameNum, &fval);

    int residue = jc->residue;
    // Output
    if (outputfile_ != 0)
      outputfile_->Printf("%5i %-4s%-4s%-4s%-4s%-4s%12f%12f\n",
                         residue+1, CurrentParm_->Res(residue).c_str(),
                         (*CurrentParm_)[jc->atom[0]].c_str(), 
                         (*CurrentParm_)[jc->atom[1]].c_str(),
                         (*CurrentParm_)[jc->atom[2]].c_str(), 
                         (*CurrentParm_)[jc->atom[3]].c_str(),
                         phi*Constants::RADDEG, Jval);
  }

  return Action::OK;
} 
