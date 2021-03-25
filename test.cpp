#include "O3Triangulation.h"
#include "O3Tetrahedron.h"

#include <iostream>
#include <set>
#include <string>

void recurse(O3Triangulation *M, int max_tets, std::set<std::string> *already_seen, std::set<std::string> *result, int level);

int main()
{
  /*
  O3Triangulation M;
  O3Tetrahedron *T = M.newTetrahedron();
  T->join(O3Tetrahedron::v, T);
  T->join(O3Tetrahedron::f0, T);
  T->join(O3Tetrahedron::e, T);
  std::cout << "isoSig = " << M.reginaIsoSig() << "\n";
  std::cout << "O3isoSig = " << M.O3isoSig() << "\n";
  if (M.anglesAreValid())
    std::cout << "ANGLES OK\n";

  O3Triangulation N;
  O3Tetrahedron *T1 = N.newTetrahedron();
  O3Tetrahedron *T2 = N.newTetrahedron();
  T1->join(O3Tetrahedron::v, T2);
  T1->join(O3Tetrahedron::e, T2);
  T1->join(O3Tetrahedron::f1, T2);
  T1->join(O3Tetrahedron::f0, T2);
  std::cout << "isoSig = " << N.reginaIsoSig() << "\n";
  std::cout << "O3isoSig = " << N.O3isoSig() << "\n";
  if (N.anglesAreValid())
    std::cout << "Angles OK\n";
  */

  /*
  O3Triangulation M;
  O3Tetrahedron *T0 = M.newTetrahedron();
  O3Tetrahedron *T1 = M.newTetrahedron();
  O3Tetrahedron *T2 = M.newTetrahedron();

  T0->join(O3Tetrahedron::f0, T0);
  T0->join(O3Tetrahedron::v, T0);
  T0->join(O3Tetrahedron::e, T1);

  T1->join(O3Tetrahedron::v, T2);
  T1->join(O3Tetrahedron::f0, T2);

  T2->join(O3Tetrahedron::e, T2);
  T2->join(O3Tetrahedron::f0, T1);

  //std::cout << "destination sequence: " << M.dest
  std::cout << "O3IsoSig: " << M.O3isoSig() << "\n";
  */

  O3Triangulation O;
  O.newTetrahedron();
  std::set<std::string> already_seen, result;
  
  recurse(&O, 3, &already_seen, &result, 0);

  for(std::set<std::string>::const_iterator it = result.begin();
      it != result.end();
      ++it) {
    std::cout << *it << std::endl;
  }

}

void recurse(O3Triangulation *M, int max_tets, std::set<std::string> *already_seen, std::set<std::string> *result, int level)
{
  std::string padding = "";
  for (int i = 0; i < level; i++) {
    padding += "  ";
  }

  //const std::string isoSig = M->reginaIsoSig();
  const std::string isoSig = M->O3isoSig();
  //std::cout << "O3isoSig: " << isoSig << "\n";

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
    std::cout << "No open faces for " << M->reginaIsoSig() << "\n";
    std::cout << "O3IsoSig = " << M->O3isoSig() << "\n";
    if (M->anglesAreValid()) {
      std::cout << "ANGLES OK\n";
      result->insert(isoSig);
    }
    return;
  }

  /*std::cout << padding << "The list of open faces is: \n" << padding;
  for (auto iter = openFaces.begin(); iter != openFaces.end(); iter++) {
    std::cout << (*iter).first << ":" << (*iter).second << ", ";
  }
  std::cout << "\n";
  */

  int fixedTetIndex = openFaces.front().first;  
  if (M->size() < max_tets) {
    // Try all four ways of adding a new tetrahedron to the chosen tetrahedron.
    for (int f = 0; f <= 3; f++) {
      // If face f is open in the fixed tetrahedron
      if (M->tetrahedron(fixedTetIndex)->isOpen(f)) {
	// Copy the triangulation.
	O3Triangulation N(*M);

	// Create a new tetrahedron T2 and glue it to face f of an existing tetrahedron T1.
	// We choose T1 arbitrarily to be the first tetrahedron in the list of open face pairs.
	O3Tetrahedron *T1 = N.tetrahedron(openFaces.front().first);
	O3Tetrahedron *T2 = N.newTetrahedron();
	T1->join(f, T2);
	recurse(&N, max_tets, already_seen, result, level+1);
      }
    }
  }
  
  // Iterate over possible faces to glue in the fixed tetrahedron.
  for (int f=0; f <= 3; f++) {
    if (M->tetrahedron(fixedTetIndex)->isOpen(f)) {
      // Iterate over all tetrahedra T and glue to face f of the fixed tetrahedron if possible.
      for (int i=0; i < M->size(); i++) {
	// Only copy the triangulation if the correct face is open in the other tetrahedron.
	if (M->tetrahedron(i)->isOpen(f)) {
	  // Copy the triangulation
	  O3Triangulation N(*M);
	  // Get the fixed face.
	  O3Tetrahedron *T1 = N.tetrahedron(fixedTetIndex);
	  O3Tetrahedron *T2 = N.tetrahedron(i);
	  T1->join(f, T2);
	  
	  recurse(&N, max_tets, already_seen, result, level+1);
	}
      }
    }
  }
}
