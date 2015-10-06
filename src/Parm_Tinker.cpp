#include "Parm_Tinker.h"
#include "TinkerFile.h"
#include "CpptrajStdio.h"

// Parm_Tinker::ID_ParmFormat() 
bool Parm_Tinker::ID_ParmFormat(CpptrajFile& fileIn) {
  return TinkerFile::ID_Tinker(fileIn);
}
    
// Parm_Tinker::ReadParm()
/** Read file as a Tinker file. */
int Parm_Tinker::ReadParm(FileName const& fname, Topology &parmOut) {
  TinkerFile infile;
  infile.SetTinkerName( fname );
  if (infile.OpenTinker()) return 1;
  mprintf("\tReading Tinker file %s as topology file.\n",infile.Filename().base());
  // Allocate memory for coordinates.
  Frame Coords;
  std::vector<int> Bonds;
  std::vector<Atom> Atoms = infile.ReadTinkerAtoms(Coords, Bonds);
  if (Atoms.empty()) return 1;
  // Use up to first 3 chars of title as residue name.
  std::string resname;
  for (std::string::const_iterator c = infile.TinkerTitle().begin();
                                   c != infile.TinkerTitle().end(); ++c)
    resname += *c;
  if (resname.size() > 3) resname.resize(3);
  Residue tinker_res( resname, 0, ' ', ' ' );
  // Put atoms into topology
  for (std::vector<Atom>::const_iterator atom = Atoms.begin();
                                         atom != Atoms.end();
                                       ++atom)
    parmOut.AddTopAtom( *atom, tinker_res );
  // Add bond information
  for (std::vector<int>::const_iterator bond = Bonds.begin();
                                        bond != Bonds.end(); bond += 2)
    parmOut.AddBond( *bond, *(bond+1) );
  // Try to set up residue info based on bonds.
  if (parmOut.Setup_NoResInfo()) return 1;
  // Set topology box info.
  parmOut.SetParmBox( infile.TinkerBox() );
  parmOut.SetParmName( infile.TinkerTitle(), infile.Filename() );
  mprintf("\tTinker file contains %i atoms, %i residues,\n", parmOut.Natom(),parmOut.Nres());
  //mprintf("    %i bonds to H, %i other bonds.\n", parmOut.NbondsWithH,parmOut.NbondsWithoutH);

  infile.CloseFile();

  return 0;
}
