#include "O3Triangulation.h"
#include "O3Tetrahedron.h"

#include <iostream>
#include <set>
#include <string>

void recurse(O3Triangulation *M, int max_tets, std::set<std::string> *already_seen, std::set<std::string> *result);

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
  N.printDihedralAngles();

  O3Triangulation O;
  O.newTetrahedron();
  std::set<std::string> already_seen, result;
  
  recurse(&O, 2, &already_seen, &result);

  for(std::set<std::string>::const_iterator it = result.begin();
      it != result.end();
      ++it) {
    std::cout << *it << std::endl;
  }
}

void recurse(O3Triangulation *M, int max_tets, std::set<std::string> *already_seen, std::set<std::string> *result)
{
  const std::string isoSig = M->isoSig();

  if (already_seen->find(isoSig) != already_seen->end()) {
    // We've already encountered this triangulation.
    return;
  }


  already_seen->insert(isoSig);

  // Get the open faces in the triangulation.
  std::vector<std::pair<int, int>> openFaces = M->getOpenFaces();
  
  // If there are no open faces, we're done with this triangulation.
  if (openFaces.empty()) {
    // If the dihedral angles are in the permissible set, add this triangulation to our list.
    std::cout << "No open faces for " << isoSig << "\n";
    M->printDihedralAngles();
    if (M->anglesAreValid()) {
      result->insert(isoSig);
    }
    return;
  }

  //std::cout << "See a triangulation with " << M->size() << " O3Tets. The max is " << max_tets << "\n";
  std::cout << "The list of open faces is: \n";
  for (auto iter = openFaces.begin(); iter != openFaces.end(); iter++) {
    std::cout << (*iter).first << ":" << (*iter).second << ", ";
  }
  std::cout << "\n";
  if (M->size() < max_tets) {
    // Copy the triangulation.
    O3Triangulation N(*M);

    // Create a new tetrahedron and glue it to one of the existing open faces.
    // We make the arbitrary choice of choosing the first open face in the list.
    O3Tetrahedron *T1 = N.tetrahedron(openFaces.front().first);
    int face1 = openFaces.front().second;
    O3Tetrahedron *T2 = N.newTetrahedron();

    T1->join(face1, T2);

    recurse(&N, max_tets, already_seen, result);
  }

  // Iterate over possible gluings for the chosen face.
  // We've already fixed a face, so there is at most one gluing possible for each tetrahedron.
  for (int i=0; i < M->size(); i++) {
    int face1 = openFaces.front().second;
    // Only copy the triangulation if the correct face is open in the other tetrahedron.
    if (M->tetrahedron(i)->isOpen(face1)) {
      // Copy the triangulation
      O3Triangulation N(*M);
      // Get the fixed face.
      O3Tetrahedron *T1 = N.tetrahedron(openFaces.front().first);
      // For each tet, try to glue face1 of T1 to tet. If not possible,
      // the join method will not alter the triangulation. So when we
      // recurse, it will return immediately as the isoSig has already been seen.
      // this might make things slow.
      O3Tetrahedron *T2 = N.tetrahedron(i);
      T1->join(face1, T2);

      recurse(&N, max_tets, already_seen, result);
    }
  }
}
