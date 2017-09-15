#include <cmath> // DEBUG
#include <locale>
#include <stack> // Tokenize()
#include "MaskToken.h"
#include "ArgList.h" // Tokenize()
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToInteger, convertToDouble
#include "DistRoutines.h" // selection by distance

MaskToken::MaskToken() :
  type_(OP_NONE),
  res1_(-1),
  res2_(-1),
  name_(""),
  onStack_(false),
  d_within_(false),
  d_atom_(false),
  distance_(0.0)
{ }

const char* MaskToken::MaskTypeString[] = {
  "OP_NONE", "ResNum", "ResName", "AtomNum", "AtomName", "AtomType", "AtomElement", "SelectAll",
  "OP_AND", "OP_OR", "OP_NEG", "OP_DIST"
};

void MaskToken::Print() const {
  mprintf("TOKEN: [%s]",MaskTypeString[type_]);
  switch (type_) {
    case ResName:
    case AtomName: mprintf(" Name=[%s]",*name_); break;
    case ResNum:
    case AtomNum: mprintf(" First=%i  Second=%i",res1_,res2_); break;
    case OP_DIST: 
      mprintf(" within=%i  d_atom=%i  distance^2=%lf",
              (int)d_within_, (int)d_atom_, distance_);
      break;
    default: mprintf(" ");
  }
  mprintf(" OnStack=%i\n",(int)onStack_);
/*
  mprintf("TOKEN: [%s] Res1=%i  Res2=%i  Name=[%s]  OnStack=%i\n",
          MaskTypeString[type_], res1_, res2_, *name_, (int)onStack_);
  mprintf("            within=%i  d_atom=%i  distance^2=%lf\n",
          (int)d_within_, (int)d_atom_, distance_);*/
} 

const char *MaskToken::TypeName() const {
  return MaskTypeString[type_];
}

void MaskToken::MakeNameType() {
  if (type_ == ResNum)
    type_ = ResName;
  else if (type_ == AtomNum)
    type_ = AtomName;
}

/** Basic : or @ operand. */
int MaskToken::SetToken( MaskTokenType typeIn, std::string const& tokenString ) {
  std::locale loc;
  if (tokenString.empty()) return 1;
  // Set initial token type
  type_ = typeIn;
  onStack_ = false;
  // Does this token argument have an asterisk? If its at the first position
  // make this an ALL token, otherwise make this a name token.
  size_t asteriskPosition = tokenString.find_first_of("*");
  if ( asteriskPosition != std::string::npos ) {
    if (asteriskPosition == 0) {
      type_ = SelectAll;
      return 0;
    } else {
      MakeNameType();
    }
  }
  // Check that all chars are digits or - for number range 
  if (type_ == ResNum || type_ == AtomNum) {
    for (std::string::const_iterator p = tokenString.begin(); p != tokenString.end(); ++p) {
      if (*p != '-' && isalpha(*p, loc)) {
        //mprintf("DEBUG: making name type because of %c\n",*p);
        MakeNameType();
        break;
      } 
    }
  }
  if (type_ == ResNum || type_ == AtomNum) {
    // Does this token argument have a dash? Only valid for number ranges.
    size_t dashPosition = tokenString.find_first_of("-");
    if (dashPosition != std::string::npos) {
      // Get first and second args. If first arg is blank negative number specified.
      std::string arg1(tokenString.begin(), tokenString.begin()+dashPosition);
      if (arg1.empty()) {
        mprinterr("Error: Mask expressions cannot contain negative numbers (%s)\n",
                  tokenString.c_str());
        return 1;
      }
      std::string arg2(tokenString.begin()+dashPosition+1, tokenString.end());
      if (arg2.empty()) {
        mprinterr("Error: Incomplete number range given (%s).\n", tokenString.c_str());
        return 1;
      }
      res1_ = convertToInteger( arg1 );
      res2_ = convertToInteger( arg2 );
    } else {
      // Get the number arg
      res1_ = convertToInteger( tokenString );
      res2_ = res1_;
    }
    // Ensure that res1 and res2 are valid
    if (res2_ < res1_) {
      mprinterr("Error: Mask range, second num (%i) less than first (%i).\n",res2_,res1_);
      return 1;
    }
    // It is expected that number args will start from 1
    if (res1_ < 1 || res2_ < 1) {
      mprinterr("Error: One or both numbers of mask arg (%s) < 1 (%i, %i)\n",
                tokenString.c_str(), res1_,res2_);
      return 1;
    }
  } else {
    // This is a string arg.
    // Use AssignNoFormat so that * not converte to '
    //name_.AssignNoFormat( tokenString.c_str() ); // TODO: Convert directly from string
    name_ = tokenString;
  }
  return 0;
}

/** [<|>][@|:]<dist> */
int MaskToken::SetDistance(std::string &distop) {
  if (distop.empty()) return 1;
  type_ = OP_DIST;
  onStack_ = false;
  // Min size is 3 chars
  if (distop.size() < 3) {
    mprinterr("Error: Malformed distance operator [%s]\n",distop.c_str());
    return 1;
  }
  // 1st char indicates within (<) or without (>)
  if (distop[0]=='<')
    d_within_ = true;
  else if (distop[0]=='>')
    d_within_ = false;
  else {
    mprinterr("Error: Malformed distance operator: expected '<' or '>' (%c)\n",distop[0]);
    return 1;
  }
  // 2nd char indidcates atoms (@) or residues (:)
  if (distop[1]=='@')
    d_atom_ = true;
  else if (distop[1]==':')
    d_atom_ = false;
  else {
    mprinterr("Error: Malformed distance operator: expected ':' or '@' (%c)\n",distop[1]);
    return 1;
  }
  // 3rd char onwards is the distance argument
  std::string distarg(distop.begin()+2, distop.end());
  distance_ = convertToDouble(distarg);
  // Pre-square the distance
  distance_ *= distance_;
  return 0;
}

// =============================================================================
// Class: MaskTokenArray
MaskTokenArray::MaskTokenArray() : debug_(10) {}

bool MaskTokenArray::IsOperator(char op) {
  if (op=='!') return true;
  if (op=='&') return true;
  if (op=='|') return true;
  if (op=='<') return true;
  if (op=='>') return true;
  return false;
}

bool MaskTokenArray::IsOperand(char op) {
  std::locale loc;
  if (op=='*')  return true;
  if (op=='/')  return true;
  if (op=='\\')  return true;
  if (op=='%')  return true;
  if (op=='-')  return true;
  if (op=='?')  return true;
  if (op==',')  return true;
  if (op=='\'') return true;
  if (op=='.')  return true;
  if (op=='=')  return true;
  if (op=='+')  return true;
  if (isalnum(op, loc)) return true;
  return false;
}

int MaskTokenArray::OperatorPriority(char op) {
  if (op == '>') return(6);
  if (op == '<') return(6);
  if (op == '!') return(5);
  if (op == '&') return(4);
  if (op == '|') return(3);
  if (op == '(') return(2);
  if (op == '_') return(1);

  mprinterr("OperatorPriority(): unknown operator ==%c== on stack when processing atom mask",op);
  return(0);
}

/** STEP1: preprocess the input string:
  * - remove spaces 
  * - isolate 'operands' into brackets [...]
  * - split expressions of the type :1-10@CA,CB into two parts;
  *   the two parts are joined with '&' operator and (for the sake
  *   of preserving precedence of other operators) enclosed into (..)
  *   :1-10@CA,CB    is split into  (:1-10 & @CA,CB)
  * - do basic error checking
  *
  * STEP2: convert to RPN
  * - operands (part enclosed in [..]) are copied directly to 'postfix'
  * - left parentheses are always pushed onto the stack
  * - when a right parenthesis is encountered the symbol at the top
  *   of the stack is popped off the stack and copied to 'postfix'
  *   until the symbol at the top of the stack is a left parenthesis.
  *   When that occurs both parentheses are discarded
  * - if the symbol scanned from 'infix' has a higher precedence 
  *   then the symbol at the top of the stack, the symbol being 
  *   scanned is pushed onto the stack
  * -  if the precedence of the symbol being scanned is lower than 
  *   or equal to the precedence of the symbol at the top of the 
  *   stack, the stack is popped to 'postfix' until the condition
  *   holds
  * - when the terminating symbol '_' is reached on the input scan
  *   the stack is popped to 'postfix' until the terminating symbol
  *   is also reached on the stack. Then the algorithm terminates.
  * - if the top of the stack is '(' and the terminating symbol '_'
  *   is scanned, or ')' is scanned when '_' is at the top of the
  *   stack, the parentheses of the atom expression were unbalanced
  *   and an unrecoverable error has occurred.
  *
  * \authors Daniel R. Roe, based on the tokenize() routine from PTRAJ
  *          by Viktor Hornak (2003).
  */
int MaskTokenArray::Tokenize() {
  std::string infix;
  std::string buffer;
  std::locale loc;
  std::string postfix;
  std::stack<char> Stack;

  // 0 means new operand or operand was just completed, and terminated with ']', 
  // 1 means operand with ":" read,
  // 2 means operand with "@" read
  // 3 means '<' or '>' read, waiting for numbers.
  int flag = 0;

  for (std::string::iterator p = maskExpression_.begin(); p != maskExpression_.end(); p++)
  {
    // Skip spaces and newlines
    if ( isspace(*p, loc) ) continue;
    
    if ( IsOperator(*p) || *p == '(' || *p == ')' ) {
      if (flag > 0) {
        buffer += "])";
        flag = 0;
        infix += buffer;
      }

      infix += *p;

      if ( *p == '>' || *p == '<' ) {
        buffer.assign("([");
        buffer += *p;
        ++p;
        buffer += *p;
        flag = 3;
        if ( *p != ':' && *p != '@' ) {
          --p;
          mprinterr("Error: Tokenize: Wrong syntax for distance mask [%c]\n",*p);
          return 1;
        }
      }
    } else if ( IsOperand(*p) || isalnum(*p, loc) ) {
      if (flag==0) {
        buffer.assign("([");
        flag = 1;
        if ( *p != '*') {
          mprinterr("Error: Tokenize: Wrong syntax [%c]\n",*p);
          return 1;
        }
      }
      if (*p == '=') { // The AMBER9 definition of wildcard '=' is equivalent to '*'.
        if (flag > 0)
          *p = '*';
        else {
          mprinterr("Error: Tokenize: '=' not in name list syntax\n");
          return 1;
        }
      }
      buffer += *p;
    } else if ( *p == ':' ) {
      if (flag == 0) {
        buffer.assign("([:");
        flag = 1;
      } else {
        buffer += "])|([:";
        flag = 1;
      }
    } else if ( *p == '@' ) {
      if (flag == 0) {
        buffer.assign("([@");
        flag = 2;
      } else if (flag == 1) {
        buffer += "]&[@";
        flag = 2;
      } else if (flag == 2) {
        buffer += "])|([@";
        flag = 2;
      }
    } else {
      mprinterr("Error: Tokenize: Unknown symbol (%c) expression when parsing atom mask [%s]\n",
                *p, maskExpression_.c_str());
      return 1;
    }
  } // END for loop over maskExpression_
  // Terminate buffer if necessary
  if (flag > 0) {
    buffer += "])";
    infix += buffer;
  }

  // TODO: Check for malformed tokens?
  // Add terminal symbol '_', needed for RPN conversion
  infix += "_";

  if (debug_ > 0)
    mprintf("DEBUG: NEW_INFIX ==%s==\n",infix.c_str());

  // -----------------------------------
  // Convert to RPN
  //postfix.clear(); 
  // push terminal symbol '_' to stack
  Stack.push('_');

  // 1 when start with "[", 0 when finished.
  flag = 0;
  char pp = ' ';
  for (std::string::const_iterator p = infix.begin(); p != infix.end(); p++) {
    if (*p == '[') {
      postfix += *p;
      flag = 1;
    } else if (*p == ']') {
      postfix += *p;
      flag = 0;
    } else if (flag == 1) {
      postfix += *p;
    } else if (*p == '(') {
      Stack.push(*p);
    } else if (*p == ')') {
      while ( (pp = Stack.top()) != '(') {
        Stack.pop();
        if (pp == '_') {
          mprinterr("Error: Mask::ToRPN: unbalanced parentheses in expression\n");
          return 1;
        }
        postfix += pp;
      }
      Stack.pop(); // Discard '('
      // At this point both parentheses are discarded
    } else if (*p == '_') {
      while ( (pp = Stack.top()) != '_') {
        Stack.pop();
        if (pp == '(') {
          mprinterr("Error: Mask::ToRPN: unbalanced parentheses in expression\n");
          return 1;
        }
        postfix += pp;
      }
      Stack.pop(); // Discard '_'
    } else if ( IsOperator(*p) ) {
      int P1 = OperatorPriority( *p );
      int P2 = OperatorPriority( Stack.top() );
      if ( P1==0 || P2==0 ) return 1; // 0 indicates error in op
      if (P1 > P2) {
        Stack.push( *p );
      } else {
        while ( P1 <= P2 ) {
          pp = Stack.top();
          Stack.pop();
          postfix += pp;
          P1 = OperatorPriority( *p );
          P2 = OperatorPriority( Stack.top() );
          if ( P1==0 || P2==0 ) return 1; // 0 indicates error in op
        }
        Stack.push( *p );
      }
    } else {
      mprinterr("Error: ToRPN: Unknown symbol in atom mask (%c)\n", *p);
      return 1;
    } 
  } // END for loop over infix
  if (debug_ > 0)
  mprintf("DEBUG: NEW_POSTFIX ==%s==\n",postfix.c_str());

  // Convert to MaskTokens in same order. The postfix expression is composed
  // of operands enclosed within brackets, and single character operators.
  // The exception is the distance operator, which is also an operand. An
  // operand can have multiple entries separated by commas (e.g. [:1,2-7,5] 
  // has 3). Once an operand is complete the OnStack bit of the token is set
  // to indicate the mask should go on the stack for processing by operators.
  // Operators store the result of their operation on the mask on top of
  // the stack so they dont need to be pushed.
  std::string tokenString;
  MaskToken token;
  maskTokens_.clear();
  for (std::string::const_iterator p = postfix.begin(); p != postfix.end(); p++) 
  {  // Operand begins here
    if (*p == '[')
      buffer.clear();
    // Operand is completed
    else if (*p == ']') {
      //mprintf("PROCESSING MASK OPERAND [%s]\n",buffer.c_str());
      if (buffer[0]=='<' || buffer[0]=='>') {
        // Distance criterion operand/operator
        // Since operator, resulting mask doesnt need to go on stack.
        token.SetDistance( buffer );
        maskTokens_.push_back( token );
      } else if (buffer[0]=='*') {
        // Select all operand - result goes on stack
        token.SetOperator( MaskToken::SelectAll );
        maskTokens_.push_back( token );
        maskTokens_.back().SetOnStack();
      } else {
        // Basic Operand list. After last entry in list processed result goes
        // on stack.
        // Determine type from first char. Default to Num; MaskToken::SetToken
        // will convert to Name if appropriate.
        MaskToken::MaskTokenType tokenType = MaskToken::OP_NONE;
        if (buffer[0]==':') // Residue
          tokenType = MaskToken::ResNum; 
        else if (buffer[0]=='@') { // Atom
          tokenType = MaskToken::AtomNum; 
          if      (buffer[1]=='%') tokenType = MaskToken::AtomType;
          else if (buffer[1]=='/') tokenType = MaskToken::AtomElement;
        }
        if (tokenType==MaskToken::OP_NONE) {
          mprinterr("Error: Unrecognized token type.\n");
          maskTokens_.clear();
          return 1;
        }
        // Create new string without type character(s)
        if (tokenType==MaskToken::ResNum || tokenType==MaskToken::AtomNum)
          tokenString.assign( buffer.begin()+1, buffer.end() );
        else
          tokenString.assign( buffer.begin()+2, buffer.end() );
        if (tokenString.empty()) {
          mprinterr("Error: empty token for '%c'\n",buffer[0]);
          return 1;
        }
        // DEBUG
        //mprintf("DEBUG: buffer=[%s]  tokenString=[%s]\n",buffer.c_str(),tokenString.c_str());
        // Split operand by comma
        ArgList commaList(tokenString, ",");
        //commaList.PrintList();
        // Assign each comma-separated arg to a new token
        for (ArgList::const_iterator arg = commaList.begin(); arg != commaList.end(); ++arg) {
          if (token.SetToken( tokenType, *arg ))
            return 1;
          maskTokens_.push_back( token );
        }
        // Indicate that after last token is processed the resulting mask should 
        // go on the stack.
        maskTokens_.back().SetOnStack();
      }
    // operand is a part inside [...]
    } else if ( IsOperand( *p ) || *p == ':' || *p == '@' || *p == '<' || *p == '>' ) {
      buffer += *p;
    // Operators
    } else if (*p == '!' ) {
      token.SetOperator( MaskToken::OP_NEG );
      maskTokens_.push_back( token );
    } else if (*p == '&' ) {
      token.SetOperator( MaskToken::OP_AND );
      maskTokens_.push_back( token );
    } else if (*p == '|' ) {
      token.SetOperator( MaskToken::OP_OR );
      maskTokens_.push_back( token );
    // Distance operator; No longer used, operand is the operator
    } else if (*p == '<' || *p == '>') {
      continue;
    } else {
      mprinterr("Error: Unknown symbol while evaluating mask (%c)\n",*p);
      maskTokens_.clear();
      return 1;
    }
  } // END loop over postfix
  // Test that operators will work correctly.
  if (!maskTokens_.empty()) {
    std::stack<char> tempStack;
    bool validMask = true;
    std::vector<MaskToken>::const_iterator T = maskTokens_.begin();
    for ( ; T != maskTokens_.end(); ++T)
    {
      if (T->Type() == MaskToken::OP_AND ||
          T->Type() == MaskToken::OP_OR) // Requires 2 operands
      {
        if (tempStack.empty()) { validMask = false; break; }
        tempStack.pop();
        if (tempStack.empty()) { validMask = false; break; }
      } else if (T->Type() == MaskToken::OP_NEG  ||
                 T->Type() == MaskToken::OP_DIST) // Requires 1 operand
      {
        if (tempStack.empty()) { validMask = false; break; }
      }
      if ( T->OnStack() )
        tempStack.push(SelectedChar_);
    } 
    if ( !validMask ) {
      mprinterr("Error: Misplaced operator %s.\n", T->TypeName());
      maskTokens_.clear();
      return 1;
    } 
  }
  if (debug_ > 0)
    for (std::vector<MaskToken>::const_iterator T = maskTokens_.begin(); 
                                                T != maskTokens_.end(); T++)
      T->Print();

  return 0;
}

/** Take the given mask expression and preprocess it for subsequent use
  * with the mask parser. Convert to infix, then postfix notation.
  */
int MaskTokenArray::SetMaskString(const char* maskStringIn) {
  if (maskStringIn!=0)
    maskExpression_.assign( maskStringIn );
  else
    maskExpression_.assign( "*" );

  if (debug_ > 0) mprintf("expression: ==%s==\n", maskExpression_.c_str());

  // Convert mask expression to maskTokens 
  if (Tokenize()) return 1;

  return 0;
}

/** If the input string is not empty, set the AtomMask expression to 
  * input string.
  * \return 0 if string was set and successfully tokenized.
  * \return 1 if tokenization failed.
  */
int MaskTokenArray::SetMaskString( std::string const& maskStringIn ) {
  if (!maskStringIn.empty())
    return ( SetMaskString( maskStringIn.c_str() ) );
  else
    return ( SetMaskString( 0 ) );
}

void MaskTokenArray::MaskInfo() const {
  mprintf("\tMask [%s] corresponds to %i atoms.\n", maskExpression_.c_str(), Nselected());
}

void MaskTokenArray::BriefMaskInfo() const {
  mprintf(" [%s](%i)", maskExpression_.c_str(), Nselected());
}

// =============================================================================
char MaskTokenArray::SelectedChar_ = 'T';
char MaskTokenArray::UnselectedChar_ = 'F';

/** \return Array of char, same size as atoms_, with T for selected atoms and F otherwise.
  */
char* MaskTokenArray::ParseMask(std::vector<Atom> const& atoms_,
                                std::vector<Residue> const& residues_,
                                const double* XYZ) const
{
  std::stack<char*> Stack;
  char *pMask = 0; 
  char *pMask2 = 0;
  int err = 0;

  for (token_iterator token = maskTokens_.begin(); token != maskTokens_.end(); ++token)
  {
    if (pMask==0) {
      // Create new blank mask
      pMask = new char[ atoms_.size() ];
      std::fill(pMask, pMask + atoms_.size(), UnselectedChar_);
    }
    switch ( token->Type() ) {
      case MaskToken::ResNum : 
        MaskSelectResidues( residues_, token->Res1(), token->Res2(), pMask );
        break;
      case MaskToken::ResName :
        MaskSelectResidues( residues_, token->Name(), pMask );
        break;
      case MaskToken::AtomNum :
        MaskSelectAtoms( atoms_, token->Res1(), token->Res2(), pMask );
        break;
      case MaskToken::AtomName :
        MaskSelectAtoms( atoms_, token->Name(), pMask );
        break;
      case MaskToken::AtomType :
        MaskSelectTypes( atoms_, token->Name(), pMask );
        break;
      case MaskToken::AtomElement :
        MaskSelectElements( atoms_, token->Name(), pMask );
        break;
      case MaskToken::SelectAll :
        std::fill(pMask, pMask + atoms_.size(), SelectedChar_);
        break;
      case MaskToken::OP_AND :
        pMask2 = Stack.top();
        Stack.pop();
        Mask_AND( Stack.top(), pMask2, atoms_.size() );
        delete[] pMask2;
        break;
      case MaskToken::OP_OR :
        pMask2 = Stack.top();
        Stack.pop();
        Mask_OR( Stack.top(), pMask2, atoms_.size() );
        delete[] pMask2;
        break;
      case MaskToken::OP_NEG :
        Mask_NEG( Stack.top(), atoms_.size() );
        break;
      case MaskToken::OP_DIST :
        err = Mask_SelectDistance( XYZ, Stack.top(), *token, atoms_, residues_); 
        break;
      default:
        mprinterr("Error: Invalid mask token (Mask [%s], type [%s]).\n",
                  MaskString(), token->TypeName() );
    }
    // If an error occurred, exit the loop.
    if (err != 0 ) break;
    // Check if this mask should now go on the stack
    if ( token->OnStack() ) {
      mprintf("DEBUG: Pushing Mask on stack, last Token [%s]\n",token->TypeName());
      Stack.push( pMask );
      pMask = 0;
    }
  }
  // If pMask is not null it is probably a blank leftover
  if (pMask!=0) delete[] pMask;

  // If stack is empty here there was an error.
  if (Stack.empty()) {
    mprinterr("Error: Could not parse mask [%s].\n", MaskString());
    return 0;
  }

  // Top of the stack should point to the final mask
  pMask = Stack.top();
  Stack.pop();
  // Stack should be empty now
  if (!Stack.empty()) {
    mprinterr("Error: Mask stack is not empty.\n");
    while (!Stack.empty()) {
      delete[] Stack.top();
      Stack.pop();
    }
    delete[] pMask;
    return 0;
  }
  return pMask;
}

// Topology::Mask_SelectDistance()
int MaskTokenArray::Mask_SelectDistance(const double* REF, char *mask,
                                        MaskToken const& token,
                                        AtomArrayT const& atoms_,
                                        ResArrayT const& residues_) const
{
  int endatom, resi;
  bool selectresidue;
  int atomi, idx, atomj;
  double d2;
  double dcut2 = token.Distance();
  const double* i_crd;

  if (REF == 0) {
    mprinterr("Error: No reference set, cannot select by distance.\n");
    return 1;
  }
  // Distance has been pre-squared.
  // Create temporary array of atom #s currently selected in mask. Also
  // reset mask, it will be the output mask.
  std::vector<unsigned int> selected;
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    if (mask[i]==SelectedChar_) {
      selected.push_back( i );
      mask[i] = UnselectedChar_;
    }
  }
  if (selected.empty()) {
    mprinterr("Error: Mask_SelectDistance: No atoms in prior selection.\n");
    return 1;
  }
  if (debug_ > 1) {
    mprintf("\t\tDistance Op: Within=%i  byAtom=%i  distance^2=%lf\n",
            (int)token.Within(), (int)token.ByAtom(), token.Distance());
    mprintf("\t\tSearch Mask=[");
    for (std::vector<unsigned int>::const_iterator at = selected.begin();
                                                   at != selected.end(); ++at)
      mprintf(" %u",*at + 1);
    mprintf(" ]\n");
  }

  if (token.ByAtom()) { // Select by atom
    // Loop over all atoms
    int n_of_atoms = (int)atoms_.size();
    if (token.Within()) {
      // Select all atoms within dcut2 of any selected atom.
#     ifdef _OPENMP
#     pragma omp parallel private(atomi, idx, atomj, d2, i_crd)
      {
#     pragma omp for
#     endif
      for (atomi = 0; atomi < n_of_atoms; atomi++) {
        // Starts out unselected
        mask[atomi] = UnselectedChar_;
        i_crd = REF + (atomi * 3);
        // Loop over initially selected atoms
        for (idx = 0; idx < (int)selected.size(); idx++) {
          atomj = selected[idx];
          d2 = DIST2_NoImage(i_crd, REF + (atomj*3));
          if (d2 < dcut2) {
            mask[atomi] = SelectedChar_;
            break;
          }
        } // END loop over initially selected atoms
      } // END loop over all atoms
#     ifdef _OPENMP
      } // END pragma omp parallel
#     endif
    } else {
      // Select all atoms outside dcut2 of all selected atoms
#     ifdef _OPENMP
#     pragma omp parallel private(atomi, idx, atomj, d2, i_crd)
      {
#     pragma omp for
#     endif
      for (atomi = 0; atomi < n_of_atoms; atomi++) {
        // Starts out selected
        mask[atomi] = SelectedChar_;
        i_crd = REF + (atomi * 3);
        // Loop over initially selected atoms
        for (idx = 0; idx < (int)selected.size(); idx++) {
          atomj = selected[idx];
          d2 = DIST2_NoImage(i_crd, REF + (atomj*3));
          if (d2 < dcut2) {
            mask[atomi] = UnselectedChar_;
            break;
          }
        } // END loop over initially selected atoms
      } // END loop over all atoms
#     ifdef _OPENMP
      } // END pragma omp parallel
#     endif
    }
  } else { // Select by residue
    int n_of_res = (int)residues_.size();
#ifdef _OPENMP
#pragma omp parallel private(atomi, idx, atomj, d2, resi, selectresidue, endatom, i_crd)
{
#pragma omp for
#endif
    for (resi = 0; resi < n_of_res; resi++) {
      selectresidue = false;
      // Determine end atom for this residue
      endatom = residues_[resi].LastAtom();
      // Loop over mask atoms
      for (idx = 0; idx < (int)selected.size(); idx++) {
        atomj = selected[idx];
        i_crd = REF + (atomj * 3);
        // Loop over residue atoms
        for (atomi = residues_[resi].FirstAtom(); atomi < endatom; atomi++) {
          d2 = DIST2_NoImage(REF + (atomi * 3), i_crd);
          if (token.Within()) {
            if (d2 < dcut2) selectresidue = true;
          } else {
            if (d2 > dcut2) selectresidue = true;
          }
          if (selectresidue) break; 
        }
        if (selectresidue) break;
      }
      if (selectresidue) {
        for (atomi = residues_[resi].FirstAtom(); atomi < endatom; atomi++)
          mask[atomi] = SelectedChar_;
        continue;
      }
    }
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  }
  return 0;
}

// MaskTokenArray::Mask_AND()
void MaskTokenArray::Mask_AND(char *mask1, char *mask2, unsigned int N) const {
  //mprintf("\t\t\tPerforming AND on masks.\n");
  for (unsigned int i = 0; i != N; i++) {
    //mprintf(" [%c|%c]",mask1[i],mask2[i]);
    if (mask1[i]==UnselectedChar_ || mask2[i]==UnselectedChar_)
      mask1[i] = UnselectedChar_;
    // Otherwise mask1 should already be T
  }
  //mprintf("\n");
}

// MaskTokenArray::Mask_OR()
void MaskTokenArray::Mask_OR(char *mask1, char *mask2, unsigned int N) const {
  //mprintf("\t\t\tPerforming OR on masks.\n");
  for (unsigned int i = 0; i != N; i++) {
    if (mask1[i]==SelectedChar_ || mask2[i]==SelectedChar_)
      mask1[i] = SelectedChar_;
    else
      mask1[i] = UnselectedChar_;
  }
}

// MaskTokenArray::Mask_NEG()
void MaskTokenArray::Mask_NEG(char *mask1, unsigned int N) const {
  //mprintf("\t\t\tNegating mask.\n");
  for (unsigned int i = 0; i != N; i++) {
    if (mask1[i]==SelectedChar_)
      mask1[i] = UnselectedChar_;
    else
      mask1[i] = SelectedChar_;
  }
}

// MaskTokenArray::MaskSelectResidues()
void MaskTokenArray::MaskSelectResidues(ResArrayT const& residues_, NameType const& name,
                                        char *mask) const
{
  //mprintf("\t\t\tSelecting residues named [%s]\n",*name);
  for (std::vector<Residue>::const_iterator res = residues_.begin();
                                            res != residues_.end(); res++)
  {
    if ( res->Name().Match( name ) ) {
      std::fill(mask + res->FirstAtom(), mask + res->LastAtom(), SelectedChar_);
    }
  }
}

// MaskTokenArray::MaskSelectResidues()
// Mask args expected to start from 1
void MaskTokenArray::MaskSelectResidues(ResArrayT const& residues_, int res1, int res2,
                                        char *mask) const
{
  int endatom;
  int nres = (int) residues_.size();
  //mprintf("\t\t\tSelecting residues %i to %i\n",res1,res2);
  // Check start atom. res1 and res2 are checked by MaskToken
  if (res1 > nres) {
    if (debug_>0)
      mprintf("Warning: Select residues: res 1 out of range (%i)\n",res1);
    return;
  }
  // If last res > nres, make it nres
  if ( res2 >= nres )
    endatom = residues_.back().LastAtom();
  else
    endatom = residues_[res2-1].LastAtom();
  // Select atoms
  std::fill(mask + residues_[res1-1].FirstAtom(), mask + endatom, SelectedChar_);
}

// MaskTokenArray::MaskSelectElements()
void MaskTokenArray::MaskSelectElements(AtomArrayT const& atoms_, NameType const& element,
                                        char* mask ) const
{
  unsigned int m = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
                                         atom != atoms_.end(); ++atom, ++m)
  {
    NameType atom_element( atom->ElementName() );
    if ( atom_element.Match( element ) )
      mask[m] = SelectedChar_;
  } 
}

// MaskTokenArray::MaskSelectTypes()
void MaskTokenArray::MaskSelectTypes(AtomArrayT const& atoms_, NameType const& type,
                                     char* mask ) const
{
  unsigned int m = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
                                         atom != atoms_.end(); ++atom, ++m)
  {
    if ( atom->Type().Match( type ) )
      mask[m] = SelectedChar_;
  } 
}

// MaskTokenArray::MaskSelectAtoms()
void MaskTokenArray::MaskSelectAtoms(AtomArrayT const& atoms_, NameType const& name,
                                     char *mask) const
{
  //mprintf("\t\t\tSelecting atoms named [%s]\n",*name);
  unsigned int m = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
                                         atom != atoms_.end(); atom++, ++m)
  {
    //mprintf("\t\t\t%u PARM[%s]  NAME[%s]",m,(*atom).c_str(),*name);
    if ( atom->Name().Match( name ) )
      mask[m] = SelectedChar_;
    //mprintf(" %c\n",mask[m]);
  }
}

// Mask args expected to start from 1
void MaskTokenArray::MaskSelectAtoms(AtomArrayT const& atoms_, int atom1, int atom2,
                                     char *mask) const
{
  int startatom, endatom;
  //mprintf("\t\t\tSelecting atoms %i to %i\n",atom1,atom2);
  if (atom1 > (int)atoms_.size()) {
    if (debug_>0) 
      mprintf("Warning: Select atoms: atom 1 out of range (%i)\n",atom1);
    return;
  }
  startatom = atom1 - 1;
  if (atom2 > (int)atoms_.size()) 
    //mprinterr("Error: Select atoms: atom 2 out of range (%i)\n",atom2)
    endatom = atoms_.size();
  else
    endatom = atom2;
  // Select atoms
  std::fill(mask + startatom, mask + endatom, SelectedChar_);
}
