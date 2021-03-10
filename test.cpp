#include "O3Triangulation.h"
#include "O3Tetrahedron.h"

#include <iostream>



int main()
{
  O3Triangulation M;
  O3Tetrahedron *T = M.newTetrahedron();
  T->join(O3Tetrahedron::v, T);
  T->join(O3Tetrahedron::f0, T);
  T->join(O3Tetrahedron::e, T);
  std::cout << M.isoSig() << "\n";
  //M.printDihedralAngles();
  if (M.anglesAreValid())
    std::cout << "ANGLES OK\n";

  O3Triangulation N;
  O3Tetrahedron *T1 = N.newTetrahedron();
  O3Tetrahedron *T2 = N.newTetrahedron();
  T1->join(O3Tetrahedron::v, T2);
  T1->join(O3Tetrahedron::e, T2);
  T1->join(O3Tetrahedron::f1, T2);
  T1->join(O3Tetrahedron::f0, T2);
  std::cout << N.isoSig() << "\n";
  if (N.anglesAreValid())
    std::cout << "Angles OK\n";
  //N.printDihedralAngles();
  
}
