#include <cmath>
#include <sstream>
#include <locale>
#include <stack>
#include "RPNcalc.h"
#include "DataSet_1D.h"
#include "CpptrajStdio.h"

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
      else if (*ptr == '=')
        O1.SetType(OP_ASSIGN);
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

double RPNcalc::DoOperation(double d1, double d2, TokenType op_type) {
  switch (op_type) {
    case OP_MINUS: return d2 - d1;
    case OP_PLUS: return d2 + d1;
    case OP_DIV: return d2 / d1;
    case OP_MULT: return d2 * d1;
    case OP_POW: return pow(d2, d1);
    case OP_NEG: return -d1;
    // ---------------------
    case FN_SQRT: return sqrt(d1);
    case FN_EXP: return exp(d1);
    case FN_LN: return log(d1);
    default:
      mprinterr("Error: Invalid token type.\n");
  }
  return 0.0;
}

int RPNcalc::Evaluate(DataSetList& DSL) const {
  if (tokens_.empty()) {
    mprinterr("Error: Expression was not set.\n");
    return 1;
  }
  std::stack<ValType> Stack;
  ValType Dval[2];
  DataSetList LocalList;
  // Are we going to be assigning this?
  DataSet* output = 0;
  if (tokens_.front().IsValue() && tokens_.back().Type() == OP_ASSIGN) {
    if (tokens_.size() < 3) {
      mprinterr("Error: Cannot assign nothing.\n");
      return 1;
    }
    if (tokens_.front().Type() != VARIABLE) {
      mprinterr("Error: Must assign to a data set.\n");
      return 1;
    }
    output = DSL.AddSet(DataSet::DOUBLE, tokens_.front().Name(), "CALC");
    if (output == 0) return 1;
  }
    
  for (Tarray::const_iterator T = tokens_.begin(); T != tokens_.end(); ++T)
  {
    if ( T->Type() == NUMBER )
      Stack.push( ValType(T->Value()) );
    else if ( T->Type() == VARIABLE ) {
      DataSet* ds = 0;
      if (output != 0 && T == tokens_.begin())
        ds = output;
      else
        ds= DSL.GetDataSet( T->Name() );
      if (ds == 0) {
        mprinterr("Error: Data set with name '%s' not found.\n", T->name());
        return 1;
      }
      Stack.push( ValType( ds ) );
    } else {
      Dval[0].Reset();
      Dval[1].Reset();
      // Operand or function. Get operand(s)
      unsigned int nOps = (unsigned int)T->numOperands();
      if (Stack.size() < nOps) {
        mprinterr("Error: Not enough operands for '%s'.\n", T->Description());
        return 1;
      }
      for (unsigned int i = 0; i != nOps; i++) {
        Dval[i] = Stack.top();
        Stack.pop();
      }
      if (T->Type() == OP_ASSIGN) {
        // Assignment. This should be the last operation.
        if (!Dval[1].IsDataSet()) {
          mprinterr("Error: Attempting to assign to something that is not a data set.\n");
          return 1;
        }
        if (Dval[1].DS() != output) {
          mprinterr("Internal Error: Assigning to wrong data set!\n");
          return 1;
        }
        if (Dval[0].IsDataSet()) {
          if (debug_>0)
            mprintf("DEBUG: Assigning '%s' to '%s'\n", Dval[0].DS()->Legend().c_str(),
                    Dval[1].DS()->Legend().c_str());
          // Should be 1D by definition, allocated below in LocalList
          DataSet_1D const& D1 = static_cast<DataSet_1D const&>( *Dval[0].DS() );
          for (unsigned int n = 0; n != D1.Size(); n++) {
            double dval = D1.Dval(n); // TODO: Direct copy
            output->Add(n, &dval);
          }
        } else {
          if (debug_>0)
            mprintf("DEBUG: Assigning %f to '%s'\n", Dval[0].Value(),
                    Dval[1].DS()->Legend().c_str());
          double dval = Dval[0].Value();
          output->Add(0, &dval); 
        }
        Stack.push(ValType(output));
      } else if (!Dval[0].IsDataSet() && !Dval[1].IsDataSet()) {
        // Neither operand is a data set
        if (debug_>0)
          mprintf("DEBUG: '%f' [%s] '%f'\n", Dval[1].Value(), T->Description(), Dval[0].Value());
        Stack.push(ValType(DoOperation(Dval[0].Value(), Dval[1].Value(), T->Type())));
      } else {
        // One or both operands is a DataSet
        // Check that final output data set has been allocated.
        if (output == 0) {
          mprinterr("Error: DataSet math must be assigned to a variable.\n");
          return 1;
        }
        // Set up temporary data set
        DataSet* tempDS = LocalList.AddSetIdx(DataSet::DOUBLE, "TEMP", T-tokens_.begin());
        if (tempDS == 0) return 1;
        // Handle 2 or 1 operand
        if (T->numOperands() == 2) {
          if (Dval[0].IsDataSet() && Dval[1].IsDataSet()) {
            // Both are DataSets. Must have same size.
            DataSet* ds1 = Dval[0].DS();
            DataSet* ds2 = Dval[1].DS();
            if (debug_>0)
              mprintf("DEBUG: '%s' [%s] '%s' => '%s'\n", ds2->Legend().c_str(), T->Description(),
                      ds1->Legend().c_str(), tempDS->Legend().c_str());
            if (ds1->Size() != ds2->Size()) {
              mprinterr("Error: Sets '%s' and '%s' do not have same size, required for %s\n",
                        ds1->Legend().c_str(), ds2->Legend().c_str(), T->name());
              return 1;
            }
            if (ds1->Ndim() != 1 && ds2->Ndim() != 1) {
              mprinterr("Error: Data set math currently restricted to 1D data sets.\n");
              return 1;
            }
            DataSet_1D const& D1 = static_cast<DataSet_1D const&>( *ds1 );
            DataSet_1D const& D2 = static_cast<DataSet_1D const&>( *ds2 );
            for (unsigned int n = 0; n != D1.Size(); n++) {
              double dval = DoOperation(D1.Dval(n), D2.Dval(n), T->Type());
              tempDS->Add(n, &dval);
            }
          } else {
            // DataSet OP Value or Value OP DataSet
            if (Dval[0].IsDataSet()) {
              // DataSet OP Value
              DataSet* ds1 = Dval[0].DS();
              if (debug_ > 0)
                mprintf("DEBUG: %f [%s] '%s' => '%s'\n", Dval[1].Value(), T->Description(),
                        ds1->Legend().c_str(), tempDS->Legend().c_str());
              if (ds1->Ndim() != 1) {
                mprinterr("Error: Data set math currently restricted to 1D data sets.\n");
                return 1;
              }
              DataSet_1D const& D1 = static_cast<DataSet_1D const&>( *ds1 );
              double d2 = Dval[1].Value();
              for (unsigned int n = 0; n != D1.Size(); n++) {
                double dval = DoOperation(D1.Dval(n), d2, T->Type());
                tempDS->Add(n, &dval);
              }
            } else {
              // Value OP DataSet
              DataSet* ds2 = Dval[1].DS();
              if (debug_ > 0)
                mprintf("DEBUG: '%s' [%s] '%f' => '%s'\n", ds2->Legend().c_str(), T->Description(),
                        Dval[0].Value(), tempDS->Legend().c_str());
              if (ds2->Ndim() != 1) {
                mprinterr("Error: Data set math currently restricted to 1D data sets.\n");
                return 1;
              }
              DataSet_1D const& D2 = static_cast<DataSet_1D const&>( *ds2 );
              double d1 = Dval[0].Value();
              for (unsigned int n = 0; n != D2.Size(); n++) {
                double dval = DoOperation(d1, D2.Dval(n), T->Type());
                tempDS->Add(n, &dval);
              }
            }
          }
        } else {
          // Only 1 operand and it is a DataSet
          DataSet* ds1 = Dval[0].DS();
          if (debug_ > 0)
            mprintf("DEBUG: [%s] '%s' => '%s'\n", T->Description(),
                    ds1->Legend().c_str(), tempDS->Legend().c_str());
          if (ds1->Ndim() != 1) {
            mprinterr("Error: Data set math currently restricted to 1D data sets.\n");
            return 1;
          }
          DataSet_1D const& D1 = static_cast<DataSet_1D const&>( *ds1 );
          for (unsigned int n = 0; n != D1.Size(); n++) {
            double dval = DoOperation(D1.Dval(n), 0.0, T->Type());
            tempDS->Add(n, &dval);
          }
        }
        Stack.push(ValType(tempDS));
      }
    }
  }
  if (Stack.size() != 1) {
    mprinterr("Error: Unbalanced expression.\n");
    return 1;
  }
  if (output == 0) mprintf("Result: %f\n", Stack.top().Value());
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
  { 0, 2, RIGHT, OP,    "Assignment"  }, // OP_ASSIGN
  { 0, 1, NO_A,  FN,    "Square root" }, // FN_SQRT
  { 0, 1, NO_A,  FN,    "Exponential" }, // FN_EXP
  { 0, 1, NO_A,  FN,    "Natural log" }, // FN_LN
  { 0, 0, NO_A,  NO_C,  "Left Par"    }, // LPAR
  { 0, 0, NO_A,  NO_C,  "Right Par"   }, // RPAR
};
