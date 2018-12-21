#include <cstring>
#include "CIFfile.h"
#include "CpptrajStdio.h"

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
  columnData_.push_back( Sarray() );
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
  // Should be as much data as there are column headers
  if (GetColumnData( columnHeaders_.size(), infile, false )) return 1;
  return 0;
} 

/** List data currently stored in the block. */
void CIFfile::DataBlock::ListData() const {
  mprintf("DataBlock: %s\n", dataHeader_.c_str());
  for (Sarray::const_iterator colname = columnHeaders_.begin();
                              colname != columnHeaders_.end(); ++colname)
    mprintf("  Col %u name: %s\n", colname - columnHeaders_.begin(), colname->c_str());
  for (std::vector<Sarray>::const_iterator rec = columnData_.begin();
                                           rec != columnData_.end(); ++rec)
  {
    mprintf("    [%u] ", rec - columnData_.begin());
    for (Sarray::const_iterator col = rec->begin();
                                col != rec->end(); ++col)
      mprintf(" {%s}", col->c_str());
    mprintf("\n");
  }
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
  return columnData_[colnum].front();
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
  for (int i = 0; i < 10; i++) {
    std::string lineIn = fileIn.GetLine();
    if (lineIn[0] == '_') ndata++;
    if (lineIn.compare(0,5,"loop_")==0) foundLoop = true;
    if (lineIn.compare(0,9,"_entry.id")==0) foundEntryID = true;
  }
  fileIn.CloseFile();
  return ( ndata > 2 && (foundLoop || foundEntryID) );
}

// CIFfile::Read()
int CIFfile::Read(FileName const& fnameIn, int debugIn) {
  if (file_.OpenFileRead( fnameIn )) return 1;
  const char* ptr = file_.Line();
  mode currentMode = UNKNOWN;
  while (ptr != 0) {
    /// There are 3 places we can be; a data block, a looped data block,
    /// or unknown.
    if ( currentMode == UNKNOWN ) {
      // Are we at a data block yet?
      if (ptr[0] == '_')
        currentMode = SERIAL;
      else if ( strncmp(ptr, "loop_", 5) == 0 )
        currentMode = LOOP;
      else
        ptr = file_.Line();

    } else if ( currentMode == SERIAL ) {
      // SERIAL data block
      DataBlock serial;
      while ( ptr != 0 && ptr[0] == '_' ) {
        serial.AddSerialDataRecord(ptr, file_);
        ptr = file_.Line();
      }
      if (debugIn > 1) serial.ListData();
      currentMode = UNKNOWN;
      if (AddDataBlock( serial )) return 1;
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
      if (AddDataBlock( loop )) return 1;
    }
  }
  if (debugIn > 0)    
    mprintf("\tCIF file '%s', %i lines.\n", file_.Filename().full(), file_.LineNumber());
  return 0;
}

// CIFfile::GetDataBlock()
CIFfile::DataBlock const& CIFfile::GetDataBlock(std::string const& header) const {
  CIF_DataType::const_iterator it = cifdata_.find( header );
  if (it == cifdata_.end()) {
    //mprinterr("Error: CIF data block '%s' not found.\n", header.c_str());
    return emptyblock;
  }
  return it->second;
}

// CIFfile::AddDataBlock()
int CIFfile::AddDataBlock( DataBlock const& block ) {
  if (block.Header().empty()) {
    mprinterr("Internal Error: Attempting to add empty CIF data block.\n");
    return 1;
  }
  CIF_DataType::const_iterator it = cifdata_.find( block.Header() );
  if (it != cifdata_.end()) {
    mprinterr("Error: Duplicate CIF block found: '%s'\n", block.Header().c_str());
    return 1;
  }
  cifdata_.insert( std::pair<std::string, DataBlock>(block.Header(), block) );
  return 0;
}
