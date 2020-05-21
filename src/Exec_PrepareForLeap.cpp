#include "Exec_PrepareForLeap.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "CharMask.h"

// Exec_PrepareForLeap::Help()
void Exec_PrepareForLeap::Help() const
{
  mprintf("\tcrdset <coords set> [frame <#>] out <file>\n"
          "\t[cysmask <mask>] [disulfidecut <cut>] [newcysname <name>]\n"
          "\t[sugarmask <mask>]\n"
          "\t[leapunitname <name>]\n"
         );
}

static inline void ChangeResName(Residue& res, NameType const& nameIn) {
  if (res.Name() != nameIn) {
    mprintf("\tChanging residue %s to %s\n", *(res.Name()), *nameIn);
    res.SetName( nameIn );
  }
}

// Exec_PrepareForLeap::Execute()
Exec::RetType Exec_PrepareForLeap::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string crdset = argIn.GetStringKey("crdset");
  if (crdset.empty()) {
    mprinterr("Error: Must specify COORDS set with 'crdset'\n");
    return CpptrajState::ERR;
  }
  DataSet* ds = State.DSL().FindSetOfGroup( crdset, DataSet::COORDINATES );
  if (ds == 0) {
    mprinterr("Error: No COORDS set found matching %s\n", crdset.c_str());
    return CpptrajState::ERR;
  }
  DataSet_Coords& coords = static_cast<DataSet_Coords&>( *((DataSet_Coords*)ds) );
  int tgtframe = argIn.getKeyInt("frame", 1) - 1;
  mprintf("\tUsing frame %i from COORDS set %s\n", tgtframe+1, coords.legend());
  if (tgtframe < 0 || tgtframe >= (int)coords.Size()) {
    mprinterr("Error: Frame is out of range.\n");
    return CpptrajState::ERR;
  }
  Frame frameIn = coords.AllocateFrame();
  coords.GetFrame(tgtframe, frameIn);

  std::string leapunitname = argIn.GetStringKey("leapunitname", "m");
  mprintf("\tUsing leap unit name: %s\n", leapunitname.c_str());

  AtomMask sugarMask;
  std::string sugarmaskstr = argIn.GetStringKey("sugarmask");
  if (!sugarmaskstr.empty()) {
    if (sugarMask.SetMaskString(sugarmaskstr)) {
      mprinterr("Error: Setting sugar mask string.\n");
      return CpptrajState::ERR;
    }
  }

  CpptrajFile* outfile = State.DFL().AddCpptrajFile(argIn.GetStringKey("out"),
                                                    "LEaP Input", DataFileList::TEXT);
  if (outfile == 0) return CpptrajState::ERR;

  // Disulfide search
  double disulfidecut = argIn.getKeyDouble("disulfidecut", 2.1);
  std::string newcysnamestr = argIn.GetStringKey("newcysname", "CYX");
  NameType newcysname(newcysnamestr);
  mprintf("\tCysteine residues involved in disulfide bonds will be changed to: %s\n", *newcysname);
  std::string cysmaskstr = argIn.GetStringKey("cysmask", ":CYS@SG");
  mprintf("\tSearching for disulfide bonds with a cutoff of %g Ang.\n", disulfidecut);
  AtomMask cysmask;
  if (cysmask.SetMaskString( cysmaskstr )) {
    mprinterr("Error: Could not set up CYS mask string %s\n", cysmaskstr.c_str());
    return CpptrajState::ERR;
  }
  if (coords.Top().SetupIntegerMask( cysmask )) return CpptrajState::ERR;
  cysmask.MaskInfo();
  if (cysmask.None())
    mprintf("Warning: No cysteine sulfur atoms selected by %s\n", cysmaskstr.c_str());
  else {
    int nDisulfides = 0;
    double cut2 = disulfidecut * disulfidecut;
    // Try to find potential disulfide sites.
    for (AtomMask::const_iterator at1 = cysmask.begin(); at1 != cysmask.end(); ++at1) {
      for (AtomMask::const_iterator at2 = at1 + 1; at2 != cysmask.end(); ++at2) {
        bool isBonded = false;
        // Check if the bond already exists
        if (coords.Top()[*at1].IsBondedTo(*at2)) {
          mprintf("\tExisting disulfide: %s to %s\n",
                  coords.Top().ResNameNumAtomNameNum(*at1).c_str(),
                  coords.Top().ResNameNumAtomNameNum(*at2).c_str());
          isBonded = true;
        } else {
          // TODO imaging?
          double r2 = DIST2_NoImage(frameIn.XYZ(*at1), frameIn.XYZ(*at2));
          if (r2 < cut2) {
            mprintf("\tPotential disulfide: %s to %s (%g Ang.)\n",
                    coords.Top().ResNameNumAtomNameNum(*at1).c_str(),
                    coords.Top().ResNameNumAtomNameNum(*at2).c_str(), sqrt(r2));
            isBonded = true;
          }
        }
        if (isBonded) {
          nDisulfides++;
          outfile->Printf("bond %s.%i.%s %s.%i.%s\n",
                          leapunitname.c_str(), coords.Top()[*at1].ResNum()+1, *(coords.Top()[*at1].Name()),
                          leapunitname.c_str(), coords.Top()[*at2].ResNum()+1, *(coords.Top()[*at2].Name()));
          ChangeResName(coords.TopPtr()->SetRes(coords.Top()[*at1].ResNum()), newcysname);
          ChangeResName(coords.TopPtr()->SetRes(coords.Top()[*at2].ResNum()), newcysname);
          // TODO add the bond?
        }
      }
    }
    mprintf("\tDetected %i disulfide bonds.\n", nDisulfides);
  }

  // Prepare sugars
  if (sugarMask.MaskStringSet()) {
    mprintf("\tPreparing sugars selected by '%s'\n", sugarMask.MaskString());
    if (coords.Top().SetupIntegerMask( sugarMask )) return CpptrajState::ERR;
    sugarMask.MaskInfo();
    if (sugarMask.None())
      mprintf("Warning: No sugar atoms selected by %s\n", sugarMask.MaskString());
    else {
      CharMask cmask( sugarMask.ConvertToCharMask(), sugarMask.Nselected() );
      // Get sugar residue numbers
      std::vector<int> sugarResNums = coords.Top().ResnumsSelectedBy( sugarMask );
      // For each sugar residue, see if it is bonded to a non-sugar residue.
      // If it is, remove that bond but record it.
      for (std::vector<int>::const_iterator rnum = sugarResNums.begin();
                                            rnum != sugarResNums.end(); ++rnum)
      {
        Residue const& Res = coords.Top().Res(*rnum);
        std::vector<int> bondedAtoms;
        for (int at = Res.FirstAtom(); at != Res.LastAtom(); at++) {
          for (Atom::bond_iterator bat = coords.Top()[at].bondbegin();
                                   bat != coords.Top()[at].bondend(); ++bat)
          {
            if (coords.Top()[*bat].ResNum() != *rnum) {
              if (!cmask.AtomInCharMask(*bat)) {
                mprintf("\tSugar %s bonded to non-sugar %s\n",
                        coords.Top().ResNameNumAtomNameNum(at).c_str(),
                        coords.Top().ResNameNumAtomNameNum(*bat).c_str());
              } else {
                 mprintf("\tSugar %s bonded to sugar %s\n",
                        coords.Top().ResNameNumAtomNameNum(at).c_str(),
                        coords.Top().ResNameNumAtomNameNum(*bat).c_str());
              }
            }
          } // END loop over bonded atoms
        } // END loop over residue atoms
      } // END loop over sugar residues
    }
  }
  
  return CpptrajState::OK;
}
