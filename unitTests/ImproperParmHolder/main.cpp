// Unit test for ImproperParmHolder class
#include <cstdio>
#include "Parm/ImproperParmHolder.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

/*static void printImp(DihedralType const& imp) {
  printf("%i %i %i %i\n", imp.A1(), imp.A2(), imp.A3(), imp.A4());
}*/

int main() {
  using namespace Cpptraj::Parm;
  TypeNameHolder ip0(4);
  ip0.AddName("N*");
  ip0.AddName("CX");
  ip0.AddName("CT");
  ip0.AddName("HN");

  TypeNameHolder ip1(4);
  ip1.AddName("O");
  ip1.AddName("HO");
  ip1.AddName("CT");
  ip1.AddName("N");
  //ip1.SortImproperByAlpha("X");

  TypeNameHolder ip1a(4);
  ip1a.AddName("N");
  ip1a.AddName("O");
  ip1a.AddName("CT");
  ip1a.AddName("HO");

  ImproperParmHolder IP0;
  Cpptraj::Parm::RetType ret = IP0.AddParm( ip1, DihedralParmType( 3.0, 1.0, 0.0 ), false );
  if (ret == Cpptraj::Parm::ERR) return Err("Could not add improper parameter");
  bool found;
  DihedralParmArray impropers = IP0.FindParam( ip1a, found );
  if (!found) return Err("Could not find improper parameter (no wildcards).");

  TypeNameHolder ip2(4);
  ip2.AddName("X");
  ip2.AddName("X");
  ip2.AddName("CT");
  ip2.AddName("O");

  ImproperParmHolder IP;
  ret = IP.AddParm( ip2, DihedralParmType( 2.0, 1.0, 3.14159/2.0 ), false );
  if (ret == Cpptraj::Parm::ERR) return Err("Could not add improper parameter");
  //ip2.SortImproperByAlpha("X");
  
  impropers = IP.FindParam( ip1, found);
  if (found) return Err("Improper parameter search before wildcard added failed.");

  IP.SetWildcard('X');
  impropers = IP.FindParam( ip1, found );
  if (!found) return Err("Improper parameter search with wildcard match failed.");

  impropers = IP.FindParam( ip0, found );
  if (found) return Err("Improper parameter search found something when it should not have.");

  DihedralType imp( 0, 1, 2, 3, -1 );
  IP.ReorderImproper( imp, ImproperParmHolder::O_013 );
  if (imp.A1() != 0 || imp.A2() != 1 || imp.A3() != 2 || imp.A4() != 3)
    return Err("Improper reorder failed (O_013).");
  IP.ReorderImproper( imp, ImproperParmHolder::O_031 );
  if (imp.A1() != 0 || imp.A2() != 3 || imp.A3() != 2 || imp.A4() != 1)
    return Err("Improper reorder failed (O_031).");
  imp = DihedralType( 0, 1, 2, 3, -1 );
  IP.ReorderImproper( imp, ImproperParmHolder::O_103 );
  if (imp.A1() != 1 || imp.A2() != 0 || imp.A3() != 2 || imp.A4() != 3)
    return Err("Improper reorder failed (O_103).");
  imp = DihedralType( 0, 1, 2, 3, -1 );
  IP.ReorderImproper( imp, ImproperParmHolder::O_130 );
  //printImp(imp);
  if (imp.A1() != 1 || imp.A2() != 3 || imp.A3() != 2 || imp.A4() != 0)
    return Err("Improper reorder failed (O_130).");
  imp = DihedralType( 0, 1, 2, 3, -1 );
  IP.ReorderImproper( imp, ImproperParmHolder::O_301 );
  //printImp(imp);
  if (imp.A1() != 3 || imp.A2() != 0 || imp.A3() != 2 || imp.A4() != 1)
    return Err("Improper reorder failed (O_301).");
  imp = DihedralType( 0, 1, 2, 3, -1 );
  IP.ReorderImproper( imp, ImproperParmHolder::O_310 );
  //printImp(imp);
  if (imp.A1() != 3 || imp.A2() != 1 || imp.A3() != 2 || imp.A4() != 0)
    return Err("Improper reorder failed (O_310).");

  return 0;
}
