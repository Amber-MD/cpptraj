#ifndef INC_RANGE_H
#define INC_RANGE_H
/*
 * Range
 */
#include <list>
class Range : public std::list<int> {
    char *rangeArg;
  public:
    Range();
    ~Range();

    int SetRange(char *);
    int SetRange(int,int);

    char *RangeArg();
    void PrintRange(const char*);
};
#endif
