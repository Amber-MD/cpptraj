#ifndef INC_RANGE_H
#define INC_RANGE_H
#include <list>
#include <string>
// Class: Range
/// The Range class is used to hold an ordered list of numbers. 
/** Range can either be set using a range expression, e.g. 
  * X-Y,C-D (X to Y and C to D) etc, or can be set using a beginning 
  * and end number.
*/
class Range {
    std::string rangeArg;
    std::list<int> rangeList;
    std::list<int>::iterator rangeNum;
  public:
    Range();
    ~Range();
    Range(const Range&);
    Range & operator=(const Range&);

    int SetRange(char *);
    int SetRange(int,int);
    void SetRange(Range*);

    char *RangeArg();
    void PrintRange(const char*,int);

    void ShiftBy(int);
    void AddToRange(int);
    void Begin();
    bool NextInRange(int *);
    bool End();
    void Next();
    int Current();
    void RemoveFromRange(int);

    bool Empty() { return rangeList.empty();      }
    int Size()   { return (int) rangeList.size(); }
};
#endif
