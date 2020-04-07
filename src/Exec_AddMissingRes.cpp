#include "Exec_AddMissingRes.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include <cstdlib> // atoi
#include <cstring> //strncmp

// Exec_AddMissingRes::Help()
void Exec_AddMissingRes::Help() const
{

}

typedef std::vector<std::string> Sarray;

/// Record a gap in PDB
class Gap {
  public:
    Gap() {}
    /// Start gap with res name, res num, and chain
    Gap(std::string const& startNameIn, int startResIn, std::string const& startChainIn) :
      resNames_(1, startNameIn), startRes_(startResIn), stopRes_(-1), chainId_(startChainIn[0])
      {}
    /// Start gap with res num and chain
    Gap(int startResIn, std::string const& startChainIn) :
      startRes_(startResIn), stopRes_(-1), chainId_(startChainIn[0])
      {}

    std::string const& FirstName() const { return resNames_.front(); }
    std::string const& LastName()  const { return resNames_.back(); }
    int StartRes()                 const { return startRes_; }
    int StopRes()                  const { return stopRes_; }
    char Chain()                   const { return chainId_; }
    unsigned int Nres()            const { return resNames_.size(); }
    typedef Sarray::const_iterator name_iterator;
    name_iterator nameBegin()      const { return resNames_.begin(); }
    name_iterator nameEnd()        const { return resNames_.end(); }

    void SetStopRes(int s) { stopRes_ = s; }
    void AddGapRes(std::string const& r) { resNames_.push_back( r ); }
  private:
    Sarray resNames_; ///< Residue names in the Gap
    int startRes_;    ///< pdb start residue number for the gap 
    int stopRes_;     ///< pdb stop residue number for the gap
    char chainId_;    ///< chain ID of the gap
};

// Exec_AddMissingRes::Execute()
Exec::RetType Exec_AddMissingRes::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string pdbname = argIn.GetStringKey("pdbname");
  if (pdbname.empty()) {
    mprinterr("Error: provide PDB name.\n");
    return CpptrajState::ERR;
  }
  std::string outname = argIn.GetStringKey("out");

  mprintf("\tPDB name: %s\n", pdbname.c_str());
  CpptrajFile outfile; // TODO use DFL
  outfile.OpenWrite(outname);
  mprintf("\tOutput file: %s\n", outfile.Filename().full());

  // Get gap info from PDB
  typedef std::vector<Gap> Garray;
  Garray Gaps;
  BufferedLine infile;
  if (infile.OpenFileRead( pdbname )) return CpptrajState::ERR;
  const char* linePtr = infile.Line();
  int inMissing = 0;
  int nmissing = 0;
  bool firstGap = true;
  bool atTheEnd = false;
  std::string lastname, lastchain;
  int lastres = 0;
  std::string currentname, currentchain;
  int currentres;
  while (linePtr != 0) {
    if (strncmp(linePtr, "REMARK", 6) == 0)
    {
      ArgList line(linePtr);
      if (line.Nargs() > 2) {
        if (inMissing == 0) {
          // MISSING section not yet encountered.
          if (line[0] == "REMARK" && line[2] == "MISSING") {
            inMissing = 1;
          }
        } else if (inMissing == 1) {
          // In MISSING, looking for start of missing residues
          if (line[0] == "REMARK" && line[2] == "M") {
            inMissing = 2;
            nmissing = 0;
            firstGap = true;
            atTheEnd = false;
          }
        } else if (inMissing == 2) {
          // Reading MISSING residues
          if (line[1] != "465") {
            atTheEnd = true;
          } else
            nmissing++;
          if (firstGap) {
            // The very first gap TODO check Gaps.empty
            Gaps.push_back( Gap(line[2], atoi(line[4].c_str()), line[3]) );
            lastname = Gaps.back().LastName();
            lastres = Gaps.back().StartRes();
            lastchain = Gaps.back().Chain();
            firstGap = false;
          } else {
            currentname = line[2];
            currentres = atoi(line[4].c_str());
            currentchain = line[3];
            if (atTheEnd || currentres - lastres > 1 || currentchain != lastchain) {
              // New sequence starting or end. Finish current.
              Gaps.back().SetStopRes(lastres);
              /*
              mprintf("  Gap %c %4s %6i to %4s %6i %6u\n",
                      Gaps.back().Chain(),
                      Gaps.back().FirstName().c_str(), Gaps.back().StartRes(),
                      Gaps.back().LastName().c_str(), Gaps.back().StopRes(),
                      Gaps.back().Nres());*/
              if (atTheEnd) {
                break;
              }
              Gaps.push_back( Gap(currentres, currentchain) );
            }
            // Continue the current sequence
            Gaps.back().AddGapRes(currentname);
            lastname = currentname;
            lastres = currentres;
            lastchain = currentchain;
          }
        } // END inMissing == 2
      } // END nargs > 2
    } // END REMARK
    linePtr = infile.Line();
  } // END while linePtr != 0

  // Printout
  for (Garray::const_iterator it = Gaps.begin(); it != Gaps.end(); ++it) {
    outfile.Printf("  Gap %c %4s %6i to %4s %6i %6u\n",
                   it->Chain(),
                   it->FirstName().c_str(), it->StartRes(),
                   it->LastName().c_str(), it->StopRes(),
                   it->Nres());
    // Print residues
    unsigned int col = 1;
    for (Gap::name_iterator name = it->nameBegin(); name != it->nameEnd(); ++name) {
      outfile.Printf("%c", Residue::ConvertResName(*name));
      col++;
      if (col > 80) {
        outfile.Printf("\n");
        col = 1;
      }
    }
    if (col > 1)
      outfile.Printf("\n");
  }
  outfile.Printf("%i missing residues.\n", nmissing);
  if (Gaps.empty()) {
    mprintf("Warning: No gaps found.\n");
    return CpptrajState::OK;
  }

  return CpptrajState::OK;
}
