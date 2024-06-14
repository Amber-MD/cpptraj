#include <cstring>
#include "ArgList.h"
#include "CIFfile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToDouble

static inline int LineError(const char* msg, int num, const char* ptr) {
  mprinterr("Error: CIF line %i: %s\n", num, msg);
  if (ptr != 0) mprinterr("Error: '%s'\n", ptr);
  return 1;
}

/** Separate data header with format '_HEADER.ID' into _HEADER and ID. */
int CIFfile::DataBlock::ParseData( std::string const& sIn, std::string& Header, 
                                   std::string& ID)
{
  size_t found = sIn.find_first_of(".");
  if (found == std::string::npos) {
    mprinterr("Error: No '.' in data record: %s\n", sIn.c_str());
    return 1;
  }
  ID = sIn.substr(found+1);
  Header = sIn.substr(0, found);
  //mprintf("DEBUG:\t\t\tHeader=%s  ID=%s", Header.c_str(), ID.c_str());
  return 0;
}

/** If this is the first entry in a data block, make this the data blocks
  * header. If not, make sure this matches the previously set header.
  */
int CIFfile::DataBlock::AddHeader(std::string const& Header) {
  if (dataHeader_.empty()) { // First entry in this data block
    dataHeader_ = Header;
  } else if (dataHeader_ != Header) {
    mprinterr("Error: Data header in CIF file changes from %s to %s\n",
              dataHeader_.c_str(), Header.c_str());
    return 1;
  }
  return 0;
}

static inline bool IsQuoteChar(char qc) {
  return (qc == '\'' || qc == '"' || qc == ';');
}

/// \return true if string has end quote. Skip any terminal whitespace.
static inline bool HasEndQuote(std::string const& strIn) {
  std::string::const_reverse_iterator it = strIn.rbegin();
  while (it != strIn.rend() && isspace(*it)) --it;
  if ( IsQuoteChar(*it) ) return true;
  return false;
}

static inline std::string RemoveEndQuote(std::string const& strIn) {
  std::string tmps = strIn.substr(0, strIn.size()-1);
  return tmps;
}

/** Split given line into a certain number of tokens. Data might be
  * split across multiple lines.
  * \param NexpectedCols Number of expected data cols.
  * \param infile File being read.
  * \param isSerial true if inside a serial data block.
  */
int CIFfile::DataBlock::GetColumnData(int NexpectedCols, BufferedLine& infile, bool isSerial)
{
  const char* SEP = " \t";
  // Allocate for a line of data
  //columnData_.push_back( Sarray() );
  // Tokenize the initial line
  int nReadCols = 0;
  int Ncols = infile.TokenizeLine(SEP);
  int idx = 0;
  bool insideQuote = false;
  bool insideSemi = false;
  while (nReadCols < NexpectedCols) {
    // Load up the next line if needed
    if (idx == Ncols) {
      if (infile.Line() == 0) break;
      Ncols = infile.TokenizeLine(SEP);
      idx = 0;
    }
    const char *tkn = infile.NextToken();
    // Skip blanks
    if (tkn == 0) continue;
    idx++;
    //mprintf("DEBUG: Token %i '%s'\n", idx, tkn);
    if (isSerial && nReadCols == 0) {
      // First column for serial data is header.id
      std::string ID, Header;
      if (ParseData( std::string(tkn), Header, ID )) return 1;
      //mprintf("  Ndata=%i  Data=%s\n", serialData.Nargs(), serialData[1].c_str());
      if (AddHeader( Header )) return 1;
      columnHeaders_.push_back( ID );
      nReadCols++;
    } else if (insideQuote) {
      // Append this to the current data column.
      columnData_.back().back().append( " " + std::string(tkn) );
      // Check for an end quote.
      if (HasEndQuote( columnData_.back().back() )) {
        // Remove that end quote.
        columnData_.back().back() = RemoveEndQuote( columnData_.back().back() );
        insideQuote = false;
        nReadCols++;
      }
    } else if (insideSemi) {
      // End if line begins with semicolon, otherwise append.
      if (tkn[0] == ';') {
        insideSemi = false;
        nReadCols++;
      } else
        columnData_.back().back().append( std::string(tkn) );
    } else {
      // Add new data column
      if (idx == 1 && tkn[0] == ';') {
        // Semicolon indicates more lines to be read.
        columnData_.back().push_back( std::string(tkn+1) );
        insideSemi = true;
      } else {
        columnData_.back().push_back( std::string(tkn) );
        // Check if column began and did not end with a quote.
        if (IsQuoteChar(columnData_.back().back()[0])) {
          // Remove leading quote.
          std::string tmps = columnData_.back().back().substr(1);
          columnData_.back().back() = tmps;
          if ( !HasEndQuote((columnData_.back().back())) ) {
            // Still need to look for the end quote.
            insideQuote = true;
          } else {
            // We have the end quote. Remove it.
            columnData_.back().back() = RemoveEndQuote( columnData_.back().back() );
          }
        }
      }
      if (!insideSemi && !insideQuote) nReadCols++;
    }
  }
  if (nReadCols != NexpectedCols) {
    mprinterr("Error: Line %i: '%s': Read %i columns, expected %i\n",
              infile.LineNumber(), dataHeader_.c_str(), nReadCols, NexpectedCols);
    return 1;
  }
  return 0;
}
      
/** Add entries to a serial data block. */
int CIFfile::DataBlock::AddSerialDataRecord( const char* ptr, BufferedLine& infile ) {
  if (ptr == 0) return 1;
  if (GetColumnData(2, infile, true)) return 1;
  return 0;
}

/** Start a serial data block. */
void CIFfile::DataBlock::StartSerialDataBlock() {
  // Allocate for a line of data
  columnData_.push_back( Sarray() );
}

/** Add column label from loop section. */
int CIFfile::DataBlock::AddLoopColumn( const char* ptr, BufferedLine& infile ) {
  if (ptr == 0) return 1;
  // Expect header.id
  int Ncols = infile.TokenizeLine(" \t");
  if ( Ncols > 1 ) {
    mprinterr("Error: Data record expected to have ID only.\n"
              "Error: '%s'\n", ptr);
    return 1;
  }
  std::string ID, Header;
  if (ParseData( std::string(infile.NextToken()), Header, ID )) return 1;
  //mprintf("\n"); // DEBUG
  if (AddHeader( Header )) return 1;
  columnHeaders_.push_back( ID );

  return 0;
}  

/** Add loop data. */
int CIFfile::DataBlock::AddLoopData( const char* ptr, BufferedLine& infile ) {
  // Allocate for a line of data
  columnData_.push_back( Sarray() );
  // Should be as much data as there are column headers
  if (GetColumnData( columnHeaders_.size(), infile, false )) return 1;
  return 0;
} 

/** List data currently stored in the block. */
void CIFfile::DataBlock::ListData() const {
  mprintf("DataBlock: %s\n", dataHeader_.c_str());
  for (Sarray::const_iterator colname = columnHeaders_.begin();
                              colname != columnHeaders_.end(); ++colname)
    mprintf("  Col %li name: %s\n", colname - columnHeaders_.begin(), colname->c_str());
  for (std::vector<Sarray>::const_iterator rec = columnData_.begin();
                                           rec != columnData_.end(); ++rec)
  {
    mprintf("    [%li] ", rec - columnData_.begin());
    for (Sarray::const_iterator col = rec->begin();
                                col != rec->end(); ++col)
      mprintf(" {%s}", col->c_str());
    mprintf("\n");
  }
}

/** Append given DataBlock to this one. */
void CIFfile::DataBlock::Append(DataBlock const& rhs) {
  //mprintf("\tAppending into block '%s'\n", dataHeader_.c_str());
  for (Sarray::const_iterator colname = rhs.columnHeaders_.begin();
                              colname != rhs.columnHeaders_.end(); ++colname)
    columnHeaders_.push_back( *colname );
  for (std::vector<Sarray>::const_iterator rec = rhs.columnData_.begin();
                                           rec != rhs.columnData_.end(); ++rec)
  {
    columnData_.push_back(Sarray());
    for (Sarray::const_iterator col = rec->begin();
                                col != rec->end(); ++col)
      columnData_.back().push_back( *col );
  }
  //mprintf("DEBUG: Post append:\n");
  //ListData(); // DEBUG
}

/** \return the index of the specified column, -1 if not present. */
int CIFfile::DataBlock::ColumnIndex(std::string const& headerIn) const {
  for (Sarray::const_iterator col = columnHeaders_.begin();
                              col != columnHeaders_.end(); ++col)
    if (headerIn == *col)
      return (int)(col - columnHeaders_.begin());
  return -1;
}

/** Unlike loop data, each column for serial data in columnData is sequential. */
std::string CIFfile::DataBlock::Data(std::string const& idIn) const {
  if (columnHeaders_.empty() || columnData_.empty()) return std::string("");
  int colnum = ColumnIndex( idIn );
  if (colnum == -1) return std::string("");
  //return columnData_[colnum].front();
  return columnData_[0][colnum];
}

// -----------------------------------------------------------------------------

/** CONSTRUCTOR */
CIFfile::CIFdata::CIFdata(std::string const& nameIn) : dataName_(nameIn) {}

// CIFfile::GetDataBlock()
CIFfile::DataBlock const& CIFfile::CIFdata::GetDataBlock(std::string const& header) const {
  CIF_DataType::const_iterator it = cifdata_.find( header );
  if (it == cifdata_.end()) {
    //mprinterr("Error: CIF data block '%s' not found.\n", header.c_str());
    return emptyblock;
  }
  return it->second;
}

// CIFfile::AddDataBlock()
int CIFfile::CIFdata::AddDataBlock( DataBlock const& block ) {
  if (block.Header().empty()) {
    mprinterr("Internal Error: Attempting to add empty CIF data block.\n");
    return 1;
  }
  CIF_DataType::iterator it = cifdata_.find( block.Header() );
  if (it != cifdata_.end()) {
    //mprinterr("Error: Duplicate CIF block found: '%s'\n", block.Header().c_str());
    //return 1;
    it->second.Append(block);
  } else
    cifdata_.insert( std::pair<std::string, DataBlock>(block.Header(), block) );
  return 0;
}

/** Print CIF data blocks */
void CIFfile::CIFdata::PrintDataBlocks() const {
  for (CIF_DataType::const_iterator it = cifdata_.begin(); it != cifdata_.end(); ++it)
  {
    mprintf("\tData block: %s\n", it->first.c_str());
    it->second.ListData();
  }
}

// -----------------------------------------------------------------------------
/// Used to return empty block for GetDataBlock
const CIFfile::DataBlock CIFfile::emptyblock = DataBlock();

/** Determine if fileIn is a CIF file. Look for entries beginning with 
  * an underscore (indicating data block), and a 'loop_' keyword or
  * '_entry.id' block.
  */
bool CIFfile::ID_CIF(CpptrajFile& fileIn) {
  // NOTE: ASSUME FILE SET UP FOR READ
  if (fileIn.OpenFile()) return false;
  int ndata = 0; // Number of '_XXX' entries seen
  bool foundLoop = false;
  bool foundEntryID = false;
  bool foundDataBlock = false;
  for (int i = 0; i < 10; i++) {
    std::string lineIn = fileIn.GetLine();
    if (lineIn[0] == '_') ndata++;
    if (lineIn.compare(0,5,"data_")==0) foundDataBlock = true;
    if (lineIn.compare(0,5,"loop_")==0) foundLoop = true;
    if (lineIn.compare(0,9,"_entry.id")==0) foundEntryID = true;
  }
  fileIn.CloseFile();
  if (foundDataBlock || foundLoop || foundEntryID) return true;
  return ( ndata > 2 && (foundLoop || foundEntryID || foundDataBlock) );
}

// CIFfile::Read()
int CIFfile::Read(FileName const& fnameIn, int debugIn) {
  if (file_.OpenFileRead( fnameIn )) return 1;
  const char* ptr = file_.Line();
  mode currentMode = UNKNOWN;
  while (ptr != 0) {
    /// There are 4 places we can be; a data_<name> statement, a data block
    /// or a looped data block corresponding to a previous data_<name>
    /// statement, or unknown.
    if (strncmp(ptr, "data_", 5) == 0) {
         ArgList dataLine(ptr, " \r\n");
         ArgList dataName(dataLine[0], "_");
         if (dataName.Nargs() < 2) {
           mprinterr("Error: malformed 'data_' name.\n");
           mprinterr("%s\n", ptr);
           return 1;
         }
         if (debugIn > 0) mprintf("\tGathering data for '%s'\n", dataName[1].c_str());
         data_.push_back(CIFdata(dataName[1]));
         ptr = file_.Line();
    } else if ( currentMode == UNKNOWN ) {
      // See if we are at a data block yet
      if (ptr[0] == '_')
        currentMode = SERIAL;
      else if ( strncmp(ptr, "loop_", 5) == 0 )
        currentMode = LOOP;
      else
        ptr = file_.Line();

    } else if ( currentMode == SERIAL ) {
      // SERIAL data block
      DataBlock serial;
      serial.StartSerialDataBlock();
      while ( ptr != 0 && ptr[0] == '_' ) {
        serial.AddSerialDataRecord(ptr, file_);
        ptr = file_.Line();
      }
      if (debugIn > 1) serial.ListData();
      currentMode = UNKNOWN;
      if (data_.back().AddDataBlock( serial )) return 1;
    } else if ( currentMode == LOOP ) {
      DataBlock loop;
      ptr = file_.Line();
      if (ptr == 0 || ptr[0] != '_')
        return LineError("In CIF file, malformed loop.", file_.LineNumber(), ptr);
      while (ptr != 0 && ptr[0] == '_') {
        loop.AddLoopColumn(ptr, file_);
        ptr = file_.Line();
      }
      // Should now be positioned at loop data
      if (ptr == 0)
        return LineError("In CIF file, no loop data.", file_.LineNumber(), ptr);
      while (ptr != 0 && ptr[0] != '_' && ptr[0] != '#') {
        loop.AddLoopData(ptr, file_);
        ptr = file_.Line();
      }
      if (debugIn > 1) loop.ListData();
      currentMode = UNKNOWN;
      if (data_.back().AddDataBlock( loop )) return 1;
    }
  }
  if (debugIn > 0)    
    mprintf("\tCIF file '%s', %i lines.\n", file_.Filename().full(), file_.LineNumber());
  if (data_.size() > 1)
    mprintf("\tCIF file contains %zu data entries.\n", data_.size());
  return 0;
}

/** Vector with DataBlocks corresponding to given header and value. */
/*CIFfile::DataBlock const& CIFfile::GetBlockWithColValue(
                                                   std::string const& header,
                                                   std::string const& col,
                                                   std::string const& value)
const
{
  std::vector<std::vector<CIFdata>::const_iterator> ret;
  for (std::vector<CIFdata>::const_iterator it = data_.begin();
                                            it != data_.end(); ++it)
  {
    DataBlock const& tempBlock = it->GetDataBlock( header );
    if (!tempBlock.empty()) {
      std::string data = tempBlock.Data( col );
      if (data == value)
        return tempBlock;
    }
  }
  return emptyblock;
}*/

/** List all data currently in the CIFfile. */
void CIFfile::ListAllData() const {
  for (std::vector<CIFdata>::const_iterator it = data_.begin();
                                            it != data_.end(); ++it)
  {
    mprintf("CIF data: %s\n", it->DataName().c_str());
    it->PrintDataBlocks();
  }
}

/** Get box info from _cell block.
  * \return 1 if box seems invalid, -1 if not box, 0 otherwise.
  */
int CIFfile::cif_Box_verbose(double* cif_box) const {
  if (cif_box == 0) {
    mprinterr("Internal Error: CIFfile::cif_Box_verbose: Null box passed in.\n");
    return 1;
  }
  int box_stat = 0;
  DataBlock const& cellblock = GetDataBlock("_cell");
  if (cellblock.empty()) {
    cif_box[0] = 0;
    cif_box[1] = 0;
    cif_box[2] = 0;
    cif_box[3] = 0;
    cif_box[4] = 0;
    cif_box[5] = 0;
    box_stat = -1;
  } else {
    cif_box[0] = convertToDouble( cellblock.Data("length_a") );
    cif_box[1] = convertToDouble( cellblock.Data("length_b") );
    cif_box[2] = convertToDouble( cellblock.Data("length_c") );
    if (cif_box[0] == 1.0 && cif_box[1] == 1.0 && cif_box[2] == 1.0) {
      mprintf("Warning: CIF cell lengths are all 1.0 Ang.;"
              " this usually indicates an invalid box.\n");
      box_stat = 1;
    }
    cif_box[3] = convertToDouble( cellblock.Data("angle_alpha") );
    cif_box[4] = convertToDouble( cellblock.Data("angle_beta" ) );
    cif_box[5] = convertToDouble( cellblock.Data("angle_gamma") );
    mprintf("\tRead cell info from CIF: a=%g b=%g c=%g alpha=%g beta=%g gamma=%g\n",
              cif_box[0], cif_box[1], cif_box[2], cif_box[3], cif_box[4], cif_box[5]);
  }
  return box_stat;
}

