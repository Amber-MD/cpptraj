#ifndef INC_RPNCALC_H
#define INC_RPNCALC_H
#include <vector>
#include <string>
/// Reverse Polish notation calculator
class RPNcalc {
  public:
    RPNcalc();
    void SetDebug(int d) { debug_ = d; }
    int ProcessExpression(std::string const&);
    int Evaluate() const;
  private:
    class Token;
    enum TokenType { NONE = 0, NUMBER, VARIABLE,
                     // Left-associative operators 
                     OP_MINUS, OP_PLUS, OP_DIV, OP_MULT, OP_POW,
                     // Right-associative operators
                     OP_NEG,
                     // Functions
                     FN_SQRT, FN_EXP, FN_LN,
                     // Parentheses (for infix conversion only)
                     LPAR, RPAR };
    enum Associativity { NO_A = 0, LEFT, RIGHT };
    enum OpClass { NO_C = 0, VALUE, OP, FN };
    /// Hold information about ops
    struct OpType {
      int priority_;
      int nOperands_;
      Associativity assoc_;
      OpClass opclass_;
      const char* description_;
    };

    typedef std::vector<Token> Tarray;
    Tarray tokens_;
    int debug_;
};
/// Hold values/operators for RPN calculator.
class RPNcalc::Token {
  public:
    Token() : type_(NONE), value_(0.0) {}
    /// CONSTRUCTOR - Numerical values
    Token(double val) : type_(NUMBER), value_(val) {}
    /// CONSTRUCTOR - Variables
    Token(std::string const& s) : type_(VARIABLE), value_(0.0), name_(s) {}
    /// CONSTRUCTOR - Operators
    Token(TokenType t) : type_(t), value_(0.0) {}
    /// Set token type. Intended for use with operator.
    void SetType(TokenType t) { type_ = t; }
    /// \return Token type.
    TokenType Type() const { return type_; }
    /// \return Token value.
    double Value() const { return value_; }
    /// \return Token variable name.
    const char* name() const { return name_.c_str(); }
    /// \return true if token is a variable or number
    inline bool IsValue() const { return (OpArray_[type_].opclass_ == VALUE); }
    /// \return true if token is an operator
    inline bool IsOperator() const { return (OpArray_[type_].opclass_ == OP); }
    /// \return true if token is a function.
    inline bool IsFunction() const { return (OpArray_[type_].opclass_ == FN); }
    /// \return true if OP is left associative.
    inline bool IsLeftAssociative() const { return (OpArray_[type_].assoc_ == LEFT); }
    /// \return string indicating token type.
    const char* Description() const { return OpArray_[type_].description_; }
    /// \return operator priority
    int Priority() const { return OpArray_[type_].priority_; }
    /// \return expected number of operands
    int numOperands() const { return OpArray_[type_].nOperands_; }
  private:
    static const OpType OpArray_[];
 
    TokenType type_;
    double value_; ///< Numerical value.
    std::string name_; ///< Variable name
};
#endif
