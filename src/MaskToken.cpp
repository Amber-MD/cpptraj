#include <locale>
#include <stack> // Tokenize()
#include "MaskToken.h"
#include "ArgList.h" // Tokenize()
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToInteger, convertToDouble
#include "DistRoutines.h" // selection by distance

MaskToken::MaskToken() :
  distance2_(0.0),
  name_(""),
  type_(OP_NONE),
  distOp_(BY_ATOM),
  idx1_(-1),
  idx2_(-1),
  onStack_(false),
  d_within_(false)
{ }

const char* MaskToken::MaskTypeString[] = {
  "OP_NONE",
  "ResNum", "ResName", "ResChain", "OriginalResNum",
  "AtomNum", "AtomName", "AtomType", "AtomElement",
  "MolNum",
  "SelectAll", "OP_AND", "OP_OR", "OP_NEG", "OP_DIST"
};

void MaskToken::Print() const {
  mprintf("TOKEN: [%s]",MaskTypeString[type_]);
  switch (type_) {
    case ResName:
    case ResChain:
    case AtomType:
    case AtomElement:
    case AtomName: mprintf(" Name=[%s]",*name_); break;
    case MolNum:
    case ResNum:
    case OresNum:
    case AtomNum: mprintf(" First=%i  Second=%i",idx1_,idx2_); break;
    case OP_DIST: 
      mprintf(" within=%i  distOp=%i  distance^2=%f",
              (int)d_within_, (int)distOp_, distance2_);
      break;
    default: mprintf(" ");
  }
  mprintf(" OnStack=%i\n",(int)onStack_);
/*
  mprintf("TOKEN: [%s] Res1=%i  Res2=%i  Name=[%s]  OnStack=%i\n",
          MaskTypeString[type_], idx1_, idx2_, *name_, (int)onStack_);
  mprintf("            within=%i  distOp=%i  distance^2=%f\n",
          (int)d_within_, (int)distOp_, distance2_);*/
} 

const char *MaskToken::TypeName() const { return MaskTypeString[type_]; }

/** Convert number type to name type if possible. */
int MaskToken::MakeNameType() {
  if (type_ == ResNum)
    type_ = ResName;
  else if (type_ == AtomNum)
    type_ = AtomName;
  else if (type_ == MolNum) {
    mprinterr("Internal Error: Molecule name not yet supported.\n");
    return 1;
  } else if (type_ == OresNum) {
    mprinterr("Internal Error: Only digits supported for original residue number.\n");
    return 1;
  }
  return 0;
}

/** Basic : @ or ^ operand. */
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
      if (MakeNameType()) return 1;
    }
  }
  // Check that all chars are digits or - for number range 
  if (type_ == ResNum || type_ == AtomNum || type_ == MolNum || type_ == OresNum) {
    for (std::string::const_iterator p = tokenString.begin(); p != tokenString.end(); ++p) {
      if (*p != '-' && isalpha(*p, loc)) {
        //mprintf("DEBUG: making name type because of %c\n",*p);
        if (MakeNameType()) return 1;
        break;
      } 
    }
  }
  if (type_ == ResNum || type_ == AtomNum || type_ == MolNum || type_ == OresNum) {
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
      idx1_ = convertToInteger( arg1 );
      idx2_ = convertToInteger( arg2 );
    } else {
      // Get the number arg
      idx1_ = convertToInteger( tokenString );
      idx2_ = idx1_;
    }
    // Ensure that res1 and res2 are valid
    if (idx2_ < idx1_) {
      mprinterr("Error: Mask range, second num (%i) less than first (%i).\n",idx2_,idx1_);
      return 1;
    }
    // It is expected that number args will start from 1
    if (idx1_ < 1 || idx2_ < 1) {
      mprinterr("Error: One or both numbers of mask arg (%s) < 1 (%i, %i)\n",
                tokenString.c_str(), idx1_,idx2_);
      return 1;
    }
  } else {
    // This is a string arg.
    // Use AssignNoFormat so that * not convert to '
    //name_.AssignNoFormat( tokenString.c_str() ); // TODO: Convert directly from string
    name_ = tokenString;
  }
  return 0;
}

/** Distance by distance. [<|>][@|:|^]<dist> */
int MaskToken::SetDistance(std::string const& distop) {
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
  // 2nd char indidcates atoms (@), residues (:), or molecules (^)
  if (distop[1]=='@')
    distOp_ = BY_ATOM;
  else if (distop[1]==':')
    distOp_ = BY_RES;
  else if (distop[1]=='^')
    distOp_ = BY_MOL;
  else {
    mprinterr("Error: Malformed distance operator: expected '^', ':', or '@' (%c)\n",distop[1]);
    return 1;
  }
  // 3rd char onwards is the distance argument
  std::string distarg(distop.begin()+2, distop.end());
  distance2_ = convertToDouble(distarg);
  // Pre-square the distance
  distance2_ *= distance2_;
  return 0;
}

// =============================================================================
// Class: MaskTokenArray
MaskTokenArray::MaskTokenArray() : debug_(10) {}

// TODO include parentheses?
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
  if (op=='*' ) return true; // Wildcard
  if (op=='/' ) return true; // Atom element, res chain ID
  if (op=='\\') return true;
  if (op=='%' ) return true; // Atom type
  if (op==';' ) return true; // Original res number
  if (op=='-' ) return true;
  if (op=='?' ) return true; // Wildcard single char
  if (op==',' ) return true;
  if (op=='\'') return true;
  if (op=='.' ) return true;
  if (op=='=' ) return true; // Wildcard
  if (op=='+' ) return true;
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

  mprinterr("OperatorPriority(): unknown operator '%c' on stack when processing atom mask",op);
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

  // 0 means new operand or operand was just completed, and terminated with ']'.
  // 1 means operand with ":" read.
  // 2 means operand with "@" read.
  // 3 means '<' or '>' read, waiting for numbers.
  // 4 means operand with "^" read.
  // 5 means operand with ":" just read; check for additional ":".
  int flag = 0;

  for (std::string::iterator p = maskExpression_.begin(); p != maskExpression_.end(); p++)
  {
    // Skip spaces and newlines
    if ( isspace(*p, loc) ) continue;
    if ( flag == 5 && *p != ':' ) flag = 1;
    if ( IsOperator(*p) || *p == '(' || *p == ')' ) {
      //mprintf("DEBUG: Operator or parentheses: %c\n", *p);
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
        if ( *p != ':' && *p != '@' && *p != '^' ) {
          --p;
          mprinterr("Error: Tokenize: Wrong syntax for distance mask [%c]\n",*p);
          return 1;
        }
      }
    } else if ( IsOperand(*p) ) {
      //mprintf("DEBUG: Operand: %c\n", *p);
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
      // Residue character
      if (flag == 0) {
        buffer.assign("([:");
        flag = 5;
      } else if (flag == 4) {
        // Molecule AND residue
        buffer += ("]&[:");
        flag = 5;
      } else if (flag == 5) {
        // Second of two ':', just append.
        buffer += *p;
        flag = 1;
      } else {
        buffer += "])|([:";
        flag = 5;
      }
    } else if ( *p == '@' ) {
      // Atom character
      if (flag == 0) {
        buffer.assign("([@");
        flag = 2;
      } else if (flag == 1) {
        // Residue AND atom
        buffer += "]&[@";
        flag = 2;
      } else if (flag == 4) {
        // Molecule AND atom
        buffer += "]&[@";
        flag = 2;
      } else if (flag == 2) {
        // Atom OR Atom
        buffer += "])|([@";
        flag = 2;
      }
    } else if ( *p == '^' ) {
      // Molecule character
      if (flag == 0) {
        buffer.assign("([^");
        flag = 4;
      } else {
        buffer += "])|([^";
        flag = 4;
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
    mprintf("DEBUG: NEW_INFIX: %s\n",infix.c_str());

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
    mprintf("DEBUG: NEW_POSTFIX: %s\n",postfix.c_str());

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
        if (buffer[0]==':') {
          // Residue
          tokenType = MaskToken::ResNum;
          if      (buffer[1] ==':') tokenType = MaskToken::ResChain;
          else if (buffer[1] == ';') tokenType = MaskToken::OresNum;
        } else if (buffer[0]=='@') {
          // Atom
          tokenType = MaskToken::AtomNum; 
          if      (buffer[1]=='%') tokenType = MaskToken::AtomType;
          else if (buffer[1]=='/') tokenType = MaskToken::AtomElement;
        } else if (buffer[0]=='^') {
          // Molecule
          tokenType = MaskToken::MolNum;
        }
        if (tokenType==MaskToken::OP_NONE) {
          mprinterr("Error: Unrecognized token type.\n");
          maskTokens_.clear();
          return 1;
        }
        // Create new string without type character(s)
        if (tokenType==MaskToken::ResNum  ||
            tokenType==MaskToken::AtomNum ||
            tokenType==MaskToken::MolNum)
          tokenString.assign( buffer.begin()+1, buffer.end() );
        else
          tokenString.assign( buffer.begin()+2, buffer.end() );
        if (tokenString.empty()) {
          mprinterr("Error: empty token for '%c'\n",buffer[0]);
          return 1;
        }
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
    } else if ( IsOperand( *p ) || *p == ':' || *p == '@' || *p == '^' || *p == '<' || *p == '>' ) {
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
    MTarray::const_iterator T = maskTokens_.begin();
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
    for (MTarray::const_iterator T = maskTokens_.begin(); T != maskTokens_.end(); ++T)
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

  if (debug_ > 0) mprintf("DEBUG: expression: %s\n", maskExpression_.c_str());

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

/** \return Array of char, same size as atoms, with T for selected atoms and F otherwise.
  */
char* MaskTokenArray::ParseMask(AtomArrayT const& atoms,
                                ResArrayT const& residues,
                                MolArrayT const& molecules,
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
      pMask = new char[ atoms.size() ];
      std::fill(pMask, pMask + atoms.size(), UnselectedChar_);
    }
    switch ( token->Type() ) {
      case MaskToken::ResNum : 
        SelectResNum( residues, token->Idx1(), token->Idx2(), pMask );
        break;
      case MaskToken::OresNum :
        SelectOriginalResNum( residues, token->Idx1(), token->Idx2(), pMask );
        break;
      case MaskToken::ResName :
        SelectResName( residues, token->Name(), pMask );
        break;
      case MaskToken::ResChain :
        SelectChainID( residues, token->Name(), pMask );
        break;
      case MaskToken::AtomNum :
        SelectAtomNum( atoms, token->Idx1(), token->Idx2(), pMask );
        break;
      case MaskToken::AtomName :
        SelectAtomName( atoms, token->Name(), pMask );
        break;
      case MaskToken::AtomType :
        SelectAtomType( atoms, token->Name(), pMask );
        break;
      case MaskToken::AtomElement :
        SelectElement( atoms, token->Name(), pMask );
        break;
      case MaskToken::MolNum :
        SelectMolNum( molecules, token->Idx1(), token->Idx2(), pMask );
        break;
      case MaskToken::SelectAll :
        std::fill(pMask, pMask + atoms.size(), SelectedChar_);
        break;
      case MaskToken::OP_AND :
        pMask2 = Stack.top();
        Stack.pop();
        Mask_AND( Stack.top(), pMask2, atoms.size() );
        delete[] pMask2;
        break;
      case MaskToken::OP_OR :
        pMask2 = Stack.top();
        Stack.pop();
        Mask_OR( Stack.top(), pMask2, atoms.size() );
        delete[] pMask2;
        break;
      case MaskToken::OP_NEG :
        Mask_NEG( Stack.top(), atoms.size() );
        break;
      case MaskToken::OP_DIST :
        err = SelectDistance( XYZ, Stack.top(), *token, atoms, residues, molecules);
        break;
      default:
        mprinterr("Error: Invalid mask token (Mask [%s], type [%s]).\n",
                  MaskString(), token->TypeName() );
    }
    // If an error occurred, exit the loop.
    if (err != 0 ) break;
    // Check if this mask should now go on the stack
    if ( token->OnStack() ) {
      //mprintf("DEBUG: Pushing Mask on stack, last Token [%s]\n",token->TypeName());
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

/** \param REF reference coordinates.
  * \param mask Initial atom selection; will be set with final output mask.
  * \param token Describe how atoms are to be selected.
  * \param atoms Atom array.
  * \param residues Residue array.
  * \return 0 if successful, 1 if an error occurs.
  */
int MaskTokenArray::SelectDistance(const double* REF, char *mask,
                                   MaskToken const& token,
                                   AtomArrayT const& atoms,
                                   ResArrayT const& residues,
                                   MolArrayT const& molecules) const
{
  if (REF == 0) {
    mprinterr("Error: No reference set, cannot select by distance.\n");
    return 1;
  }
  // Distance cutoff has been pre-squared.
  double dcut2 = token.Distance2();
  // Create temporary array of atom #s currently selected in mask.
  // These are the atoms the search is based on.
  typedef std::vector<unsigned int> Uarray;
  Uarray Idx;
  for (unsigned int i = 0; i < atoms.size(); i++) {
    if (mask[i] == SelectedChar_)
      Idx.push_back( i*3 );
  }
  if (Idx.empty()) {
    mprinterr("Error: Mask_SelectDistance: No atoms in prior selection.\n");
    return 1;
  }
  //if (debug_ > 1) {
  //  mprintf("\t\tDistance Op: Within=%i  DistOp=%i  distance^2=%f\n",
  //          (int)token.Within(), (int)token.DistOp(), token.Distance2());
  //  mprintf("\t\tSearch Mask=[");
  //  for (Uarray::const_iterator at = Idx.begin(); at != Idx.end(); ++at)
  //    mprintf(" %u",*at/3 + 1);
  //  mprintf(" ]\n");
  //}

  char char0, char1;
  if (token.Within()) {
    // Select all atoms within dcut2 of any selected atom.
    // Select all residues with atoms within dcut2 of any selected atom.
    char0 = UnselectedChar_;
    char1 = SelectedChar_;
  } else {
    // Select all atoms outside dcut2 of all selected atoms.
    // Select all residues with atoms outside dcut2 of all selected atoms.
    char0 = SelectedChar_;
    char1 = UnselectedChar_;
  }

  if (token.DistOp() == MaskToken::BY_ATOM) {
    // Select by atom
    int n_of_atoms = (int)atoms.size();
    int atomi;
#   ifdef _OPENMP
#   pragma omp parallel private(atomi)
    {
#   pragma omp for
#   endif
    for (atomi = 0; atomi < n_of_atoms; atomi++) {
      // Initial state
      mask[atomi] = char0;
      const double* i_crd = REF + (atomi * 3);
      // Loop over initially selected atoms
      for (Uarray::const_iterator idx = Idx.begin(); idx != Idx.end(); ++idx) {
        double d2 = DIST2_NoImage(i_crd, REF + *idx);
        if (d2 < dcut2) {
          // State changes
          mask[atomi] = char1;
          break;
        }
      } // END loop over initially selected atoms
    } // END loop over all atoms
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
  } else if (token.DistOp() == MaskToken::BY_RES) {
    // Select by residue
    int n_of_res = (int)residues.size();
    int resi;
    // Loop over all residues
#   ifdef _OPENMP
#   pragma omp parallel private(resi)
    {
#   pragma omp for
#   endif
    for (resi = 0; resi < n_of_res; resi++) {
      // Initial state
      char schar = char0;
      int atomi = residues[resi].FirstAtom();
      const double* i_crd = REF + (atomi * 3);
      // Loop over residue atoms
      for (; atomi != residues[resi].LastAtom(); atomi++, i_crd += 3) {
        // Loop over initially selected atoms
        for (Uarray::const_iterator idx = Idx.begin(); idx != Idx.end(); ++idx) {
          double d2 = DIST2_NoImage(i_crd, REF + *idx);
          if (d2 < dcut2) {
            // State changes
            schar = char1;
            break;
          }
        } // END loop over initially selected atoms
        if (schar == char1) break;
      } // END loop over residue atoms
      // Set residue selection status
      for (atomi = residues[resi].FirstAtom();
           atomi != residues[resi].LastAtom(); atomi++)
        mask[atomi] = schar;
    } // END loop over all residues
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
  } else {
    // Select by molecule
    if (molecules.empty()) {
      mprinterr("Error: No molecule info. Cannot select molecules by distance.\n");
      return 1;
    }
    int n_of_mol = (int)molecules.size();
    int moli;
    // Loop over all molecules
#   ifdef _OPENMP
#   pragma omp parallel private(moli)
    {
#   pragma omp for
#   endif
    for (moli = 0; moli < n_of_mol; moli++) {
      // Initial state
      char schar = char0;
      int atomi = molecules[moli].BeginAtom();
      const double* i_crd = REF + (atomi * 3);
      // Loop over molecule atoms
      for (; atomi != molecules[moli].EndAtom(); atomi++, i_crd += 3) {
        // Loop over initially selected atoms
        for (Uarray::const_iterator idx = Idx.begin(); idx != Idx.end(); ++idx) {
          double d2 = DIST2_NoImage(i_crd, REF + *idx);
          if (d2 < dcut2) {
            // State changes
            schar = char1;
            break;
          }
        } // END loop over initially selected atoms
        if (schar == char1) break;
      } // END loop over molecule atoms
      // Set molecule selection status
      for (atomi = molecules[moli].BeginAtom();
           atomi != molecules[moli].EndAtom(); atomi++)
        mask[atomi] = schar;
    } // END loop over all molecules
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif

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

/** Select by residue name. */
void MaskTokenArray::SelectResName(ResArrayT const& residues, NameType const& name,
                                   char *mask) const
{
  //mprintf("\t\t\tSelecting residues named [%s]\n",*name);
  for (ResArrayT::const_iterator res = residues.begin(); res != residues.end(); ++res)
  {
    if ( res->Name().Match( name ) ) {
      std::fill(mask + res->FirstAtom(), mask + res->LastAtom(), SelectedChar_);
    }
  }
}

/** Select by residue number. */
void MaskTokenArray::SelectResNum(ResArrayT const& residues, int res1, int res2,
                                  char *mask) const
{
  int endatom;
  int nres = (int) residues.size();
  //mprintf("\t\t\tSelecting residues %i to %i\n",res1,res2);
  // Check start atom. res1 and res2 are checked by MaskToken
  // Mask args expected to start from 1
  if (res1 > nres) {
    mprintf("Warning: Select residues: res 1 out of range (%i > %i)\n",res1, nres);
    return;
  }
  // If last res > nres, make it nres
  if ( res2 >= nres )
    endatom = residues.back().LastAtom();
  else
    endatom = residues[res2-1].LastAtom();
  // Select atoms
  std::fill(mask + residues[res1-1].FirstAtom(), mask + endatom, SelectedChar_);
}

/** Select by original residue number. */
void MaskTokenArray::SelectOriginalResNum(ResArrayT const& residues, int res1, int res2,
                                          char* mask) const
{
  for (ResArrayT::const_iterator res = residues.begin(); res != residues.end(); ++res)
    if ( res->OriginalResNum() >= res1 && res->OriginalResNum() <= res2 )
      std::fill(mask + res->FirstAtom(), mask + res->LastAtom(), SelectedChar_);
}

/** Select residues by chain ID. */
void MaskTokenArray::SelectChainID(ResArrayT const& residues, NameType const& name,
                                   char* mask) const
{
  for (ResArrayT::const_iterator res = residues.begin();
                                 res != residues.end(); ++res)
    if ( res->ChainID() == name[0] )
      std::fill(mask + res->FirstAtom(), mask + res->LastAtom(), SelectedChar_);
}

/** Select by molecule number. */
void MaskTokenArray::SelectMolNum(MolArrayT const& molecules, int mol1, int mol2,
                                  char* mask) const
{
  if (molecules.empty()) {
    mprintf("Warning: No molecule information, cannot select by molecule.\n");
    return;
  }
  int endatom;
  int nmol = (int)molecules.size();
  // Mask args expected to start from 1
  if (mol1 > nmol) {
    mprintf("Warning: Select molecules: mol 1 out of range (%i > %i)\n", mol1, nmol);
    return;
  }
  // If last mol > nmol, make it nmol
  if ( mol2 >= nmol )
    endatom = molecules.back().EndAtom();
  else
    endatom = molecules[mol2-1].EndAtom();
  // Select atoms
  std::fill(mask + molecules[mol1-1].BeginAtom(), mask + endatom, SelectedChar_);
}

/** Select by atomic element. */
void MaskTokenArray::SelectElement(AtomArrayT const& atoms, NameType const& element,
                                   char* mask ) const
{
  unsigned int m = 0;
  for (AtomArrayT::const_iterator atom = atoms.begin(); atom != atoms.end(); ++atom, ++m)
  {
    NameType atom_element( atom->ElementName() );
    if ( atom_element.Match( element ) )
      mask[m] = SelectedChar_;
  } 
}

/** Select by atom type. */
void MaskTokenArray::SelectAtomType(AtomArrayT const& atoms, NameType const& type,
                                    char* mask ) const
{
  unsigned int m = 0;
  for (AtomArrayT::const_iterator atom = atoms.begin(); atom != atoms.end(); ++atom, ++m)
  {
    if ( atom->Type().Match( type ) )
      mask[m] = SelectedChar_;
  } 
}

/** Select by atom name. */
void MaskTokenArray::SelectAtomName(AtomArrayT const& atoms, NameType const& name,
                                    char *mask) const
{
  //mprintf("\t\t\tSelecting atoms named [%s]\n",*name);
  unsigned int m = 0;
  for (AtomArrayT::const_iterator atom = atoms.begin(); atom != atoms.end(); ++atom, ++m)
  {
    //mprintf("\t\t\t%u PARM[%s]  NAME[%s]",m,(*atom).c_str(),*name);
    if ( atom->Name().Match( name ) )
      mask[m] = SelectedChar_;
    //mprintf(" %c\n",mask[m]);
  }
}

/** Select by atom number. */
void MaskTokenArray::SelectAtomNum(AtomArrayT const& atoms, int atom1, int atom2,
                                   char *mask) const
{
  int startatom, endatom;
  //mprintf("\t\t\tSelecting atoms %i to %i\n",atom1,atom2);
  // Mask args expected to start from 1
  if (atom1 > (int)atoms.size()) {
    mprintf("Warning: Select atoms: atom 1 out of range (%i > %zu)\n",atom1, atoms.size());
    return;
  }
  startatom = atom1 - 1;
  if (atom2 > (int)atoms.size())
    //mprinterr("Error: Select atoms: atom 2 out of range (%i)\n",atom2)
    endatom = atoms.size();
  else
    endatom = atom2;
  // Select atoms
  std::fill(mask + startatom, mask + endatom, SelectedChar_);
}
