#ifndef INC_MAPATOM_H
#define INC_MAPATOM_H
#include "Atom.h"
class MapAtom : public Atom {
  public :
    MapAtom();
    MapAtom(const MapAtom&);
    MapAtom(const Atom&);
    MapAtom& operator=(const MapAtom&);

    bool IsChiral() { return isChiral_; }
    bool IsMapped() { return isMapped_; }
    bool Complete() { return complete_; }
    bool IsUnique() { return (Nduplicated_ == 0); }
    const std::string& AtomID() { return atomID_; }
    const std::string& Unique() { return unique_;}
    int Nduplicated() { return Nduplicated_; }
    char CharName() { return name_; }

    void IsDuplicated() { ++Nduplicated_; }

    void SetMapped() { isMapped_ = true; } 
    void SetComplete() { complete_ = true; }
    void SetChiral() { isChiral_ = true; }
    void SetAtomID(std::string const& s) { atomID_ = s; }
    void SetUnique(std::string const& s) { unique_ = s; Nduplicated_ = 0; }
    
    void SetNotMapped() { isMapped_ = false; }
    void SetNotComplete() { complete_ = false; }
    void SetNotChiral() { isChiral_ = false; }
  private:
    static const char AtomicElementChar[];
    bool isChiral_;
    bool isMapped_;
    bool complete_;
    std::string atomID_;
    std::string unique_;
    int Nduplicated_;
    char name_;
};
#endif
