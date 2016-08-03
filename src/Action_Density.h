// -*- mode: c++; -*-

#ifndef INC_ACTION_DENSITY_H
#define INC_ACTION_DENSITY_H

#include <vector>
#include <map>

#include "Action.h"
#include "OnlineVarT.h"

#define ROUTINE_VERSION_STRING "1.0.2"

/** Calculate density along a coordinate.
  * \author Hannes H. Loeffler.
  */
class Action_Density : public Action {
public:
  Action_Density();

    DispatchObject* Alloc() const { return (DispatchObject*) new Action_Density(); }
    
    void Help() const;

private:
  Action::RetType Init(ArgList&, ActionInit&, int);
  Action::RetType Setup(ActionSetup&);
  Action::RetType DoAction(int, ActionFrame&);
  void Print();

  typedef StatsMap<long,double> statmap;
  void Output(long, long, std::vector<statmap>&);

  static const std::string emptystring;
  static const char* PropertyStr_[];
  static const char AxisStr_[];

  enum DirectionType {DX = 0, DY, DZ};
  enum PropertyType {NUMBER = 0, MASS, CHARGE, ELECTRON};

  DirectionType axis_;
  DirectionType area_coord_[2];
  PropertyType property_;

  double delta_;
  Stats<double> area_;

  std::vector<AtomMask> masks_;

  typedef std::vector<DataSet*> DSarray;
  DSarray AvSets_; ///< Hold average data sets for each mask
  DSarray SdSets_; ///< Hold SD data sets for each mask

  // std::unordered_map may be better but it is C++11, some STL's may have
  // hash_map but not same number of params,
  // two separate maps are used to ensure we store negative and positive zeros,
  // i.e. to properly bin ]-1,0[ and [0,1[
  std::vector<statmap> minus_histograms_, plus_histograms_;

  std::vector<std::vector<double> > properties_;
};
#endif    
