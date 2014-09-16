#include <cstdio> // DEBUG
#include <cstdarg> // DEBUG
#include <cmath>
#include <sstream>
#include <locale>
#include <stack>
#include "RPNcalc.h"

// DEBUG FIXME Replace with CpptrajStdio functions
static void mprintf(const char *format, ...) {
  va_list args;
  va_start(args,format);
  vfprintf(stdout,format,args);
  va_end(args);
}
static void mprinterr(const char *format, ...) {
  va_list args;
  va_start(args,format);
  vfprintf(stderr,format,args);
  va_end(args);
}

// -----------------------------------------------------------------------------
// CONSTRUCTOR
RPNcalc::RPNcalc() {}

/** Convert infix expression to RPN in tokens_ array. This uses a
  * shunting-yard algorithm which has been slightly modified to
  * recognize unary right-associative operators.
  */
int RPNcalc::ProcessExpression(std::string const& expression) {
  std::locale loc;
  if (expression.empty()) return 1;
  if (debug_ > 0) mprintf("Parsing expression: '%s'\n", expression.c_str());
  tokens_.clear();
  std::stack<Token> op_stack;
  std::string::const_iterator ptr = expression.begin();
  bool lastTokenWasOperator = true;
  while ( ptr != expression.end() )
  {
    //mprintf("DEBUG: Start of loop, char is '%c'\n", *ptr);
    // Skip whitespace
    if (isspace(*ptr, loc)) { ++ptr; continue; }

    // NUMBER ------------------------------------
    if (*ptr == '.' || isdigit(*ptr, loc))
    { // Start of a number
      std::string number;
      bool decimal_point = (*ptr == '.');
      if (decimal_point)
        number.push_back( *(ptr++) );
      while ( ptr != expression.end() && isdigit(*ptr, loc) )
      {
        number.push_back( *(ptr++) );
        // Check the next character
        if (*ptr == '.')
        {
          if (decimal_point)
          {
            mprinterr("Error: Two decimal points encountered in number: %s\n", number.c_str());
            return 1;
          }
          decimal_point = true;
          number.push_back( *(ptr++) );
        } 
      }
      if (debug_ > 0) mprintf("Number detected: %s\n", number.c_str());
      std::istringstream iss(number);
      double val;
      if (!(iss >> val)) {
        mprinterr("Error: Invalid number: %s\n", number.c_str());
        return 1;
      }
      tokens_.push_back( Token( val ) );
      lastTokenWasOperator = false;
    }
    // ALPHA (FUNCTION/VAR) ----------------------
    else if ( isalpha(*ptr, loc) )
    { // Look for Function name
      size_t pos = ptr - expression.begin();
      if (expression.compare(pos, 4, "sqrt")==0)
      {
        op_stack.push( Token(FN_SQRT) );
        ptr += 4;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 3, "exp")==0)
      {
        op_stack.push( Token(FN_EXP) );
        ptr += 3;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 2, "ln")==0)
      {
        op_stack.push( Token(FN_LN) );
        ptr += 2;
        lastTokenWasOperator = true;
      }
      // -----------------------------------------
      else
      { // Assume variable name.
        std::string varname;
        bool has_colon = false; // For index
        enum BracketState { NONE, OPEN, CLOSED };
        BracketState bracket = NONE;
        while ( ptr != expression.end() && (isalpha(*ptr,loc) || isdigit(*ptr,loc)) )
        {
          //mprintf("DEBUG: Var '%c'\n", *ptr);
          varname.push_back( *(ptr++) );
          // Check for brackets (Aspect)
          if (*ptr == '[' || *ptr == ']')
          {
            if (bracket == CLOSED)
            {
              mprinterr("Error: Multiple bracket sets encountered in set name: %s\n",
                        varname.c_str());
              return 1;
            }
            else if (*ptr == ']')
            {
              if (bracket == NONE)
              {
                mprinterr("Error: ']' encountered before '[' in set name: %s\n", varname.c_str());
                return 1;
              }
              bracket = CLOSED;
            }
            else if (*ptr == '[')
            {
              if (bracket == OPEN)
              {
                mprinterr("Error: Multiple '[' encountered in set name: %s\n", varname.c_str());
                return 1;
              }
              bracket = OPEN;
            }
            varname.push_back( *(ptr++) );
          }
          // Colon (index) can be after bracket or after name
          if (*ptr == ':')
          {
            if (has_colon)
            {
              mprinterr("Error: Two colons encountered in set name: %s\n", varname.c_str());
              return 1;
            }
            has_colon = true;
            varname.push_back( *(ptr++) );
          }
          //mprintf("DEBUG: Post Var '%c'\n", *ptr);
        }
        if (debug_>0) mprintf("Variable name detected: %s\n", varname.c_str());
        tokens_.push_back( Token(varname) );
        lastTokenWasOperator = false;
      }
    }
    // -----------------------------------------
    else if (*ptr == '(')
    {
      op_stack.push( Token(LPAR) );
      ++ptr;
    }
    else if (*ptr == ')')
    {
      if (op_stack.empty()) {
        mprinterr("Error: Mismatched parentheses.\n");
        return 1;
      }
      while (!op_stack.empty() && op_stack.top().Type() != LPAR)
      {
        tokens_.push_back( op_stack.top() );
        op_stack.pop();
      }
      if (op_stack.empty() || op_stack.top().Type() != LPAR) {
        mprinterr("Error: Mismatched parentheses.\n");
        return 1;
      }
      // Pop left parentheses off the stack.
      op_stack.pop();
      // Supposed to put function on top of stack in output
      if (!op_stack.empty() && op_stack.top().IsFunction()) {
        tokens_.push_back( op_stack.top() );
        op_stack.pop();
      }
      ++ptr;
    }
    // -----------------------------------------
    else
    { // Operator
      Token O1;
      if (*ptr == '+')
        O1.SetType(OP_PLUS);
      else if (*ptr == '-') {
        if (lastTokenWasOperator)
          O1.SetType(OP_NEG);
        else
          O1.SetType(OP_MINUS);
      } else if (*ptr == '*')
        O1.SetType(OP_MULT);
      else if (*ptr == '/')
        O1.SetType(OP_DIV);
      else if (*ptr == '^')
        O1.SetType(OP_POW);
      else {
        mprinterr("Error: Unrecognized character in expression: %c\n", *ptr);
        return 1;
      }
      ++ptr;
      while (!op_stack.empty() && op_stack.top().IsOperator())
      {
        if ((O1.IsLeftAssociative() && O1.Priority() <= op_stack.top().Priority()) ||
            (O1.Priority() < op_stack.top().Priority()))
        {
          tokens_.push_back( op_stack.top() );
          op_stack.pop();
        } else
          break;
      }
      if (debug_>0) mprintf("Operator detected: %s\n", O1.Description());
      lastTokenWasOperator = true;
      op_stack.push( O1 );
    }
  } // END loop over expression
  while (!op_stack.empty()) {
    if (op_stack.top().IsOperator()) {
      tokens_.push_back( op_stack.top() );
      op_stack.pop();
    } else {
      mprinterr("Error: Missing or mismatched parentheses.\n");
      return 1;
    }
  }
  if (debug_ > 0) {
    mprintf("Final RPN expression:\n");
    for (unsigned int i = 0; i != tokens_.size(); i++) {
      mprintf("\t%u: %s", i, tokens_[i].Description());
      if (tokens_[i].IsValue()) {
        if (tokens_[i].Type() == NUMBER)
          mprintf(" %lf", tokens_[i].Value());
        else
          mprintf(" %s", tokens_[i].name());
      }
      mprintf("\n");
    }
  }
  if (tokens_.empty()) {
    mprinterr("Error: No valid tokens detected.\n");
    return 1;
  }
  return 0;
}

inline int GetOperands(unsigned int nOps, std::stack<double>& Stack, double* Dval)
{
  if (Stack.size() < nOps) {
    mprinterr("Error: Not enough operands.\n");
    return 1;
  }
  for (unsigned int i = 0; i != nOps; i++) {
    Dval[i] = Stack.top();
    Stack.pop();
  }
  return 0;
}

int RPNcalc::Evaluate() const {
  if (tokens_.empty()) {
    mprinterr("Error: Expression was not set.\n");
    return 1;
  }
  //TODO: This may need to enapsulate data sets as well eventually.
  std::stack<double> Stack;
  double Dval[2];
  for (Tarray::const_iterator T = tokens_.begin(); T != tokens_.end(); ++T)
  {
    if ( T->Type() == NUMBER )
      Stack.push( T->Value() );
    else if ( T->Type() == VARIABLE ) {
      mprinterr("Error: Not able to handle variables yet.\n");
      return 1;
    } else {
      // Operand or function
      if (GetOperands(T->numOperands(), Stack, Dval)) return 1;
      switch (T->Type()) {
        case OP_MINUS: Stack.push( Dval[1] - Dval[0] ); break;
        case OP_PLUS: Stack.push( Dval[1] + Dval[0] ); break;
        case OP_DIV: Stack.push( Dval[1] / Dval[0] ); break;
        case OP_MULT: Stack.push( Dval[1] * Dval[0] ); break;
        case OP_POW: Stack.push( pow(Dval[1], Dval[0]) ); break;
        case OP_NEG: Stack.push( -Dval[0] ); break;
        // ---------------------
        case FN_SQRT: Stack.push( sqrt(Dval[0]) ); break;
        case FN_EXP: Stack.push( exp(Dval[0]) ); break;
        case FN_LN: Stack.push( log(Dval[0]) ); break;
        default:
          mprinterr("Error: Invalid token type: '%s'\n", T->Description());
          return 1;
      }
    }
  }
  if (Stack.size() != 1) {
    mprinterr("Error: Unbalanced expression.\n");
    return 1;
  }
  mprintf("Result: %f\n", Stack.top());
  return 0;
}

// -----------------------------------------------------------------------------
/// Priority, #operands, associativity, class, description.
const RPNcalc::OpType RPNcalc::Token::OpArray_[] = {
  { 0, 0, NO_A,  NO_C,  "None"        }, // NONE
  { 0, 0, NO_A,  VALUE, "Number"      }, // NUMBER
  { 0, 0, NO_A,  VALUE, "Variable"    }, // VARIABLE
  { 1, 2, LEFT,  OP,    "Minus"       }, // OP_MINUS
  { 1, 2, LEFT,  OP,    "Plus"        }, // OP_PLUS
  { 2, 2, LEFT,  OP,    "Divide"      }, // OP_DIV
  { 2, 2, LEFT,  OP,    "Multiply"    }, // OP_MULT
  { 3, 2, LEFT,  OP,    "Power"       }, // OP_POW
  { 4, 1, RIGHT, OP,    "Unary minus" }, // OP_NEG
  { 0, 1, NO_A,  FN,    "Square root" }, // FN_SQRT
  { 0, 1, NO_A,  FN,    "Exponential" }, // FN_EXP
  { 0, 1, NO_A,  FN,    "Natural log" }, // FN_LN
  { 0, 0, NO_A,  NO_C,  "Left Par"    }, // LPAR
  { 0, 0, NO_A,  NO_C,  "Right Par"   }, // RPAR
};
