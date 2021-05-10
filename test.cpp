#include "O3Triangulation.h"
#include "O3Tetrahedron.h"

#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include <chrono>
using namespace std::chrono;

void recurse(O3Triangulation *M, int max_tets, std::set<std::string> *already_seen, std::set<std::string> *result, int level);

int main()
{
  /*
  O3Triangulation O;
  O.newTetrahedron();
  std::set<std::string> already_seen, result;

  auto start = high_resolution_clock::now();
  recurse(&O, 10, &already_seen, &result, 0);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<minutes>(stop-start);
  std::cout << "Execution time: " << duration.count() << " minutes.\n";
  */

  // Read in a list of O3IsoSigs and construct the corresponding triangulations.
  std::vector<O3Triangulation> trigs;
  std::ifstream file("../data.txt");
  std::string sig;

  while (std::getline(file, sig)) {
    O3Triangulation M(sig);
    trigs.push_back(M);
    if (sig.compare(M.O3isoSig()) != 0) {
      std::cout << "C'est probleme!\n";
    }
  }

  for (auto trig = trigs.begin(); trig != trigs.end(); trig++) {
    std::cout << trig->O3isoSig() << "\n";
    trig->computeCuspCrossSections();
    std::cout << "\n";
    
  }
  
}

void recurse(O3Triangulation *M, int max_tets, std::set<std::string> *already_seen, std::set<std::string> *result, int level)
{
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
    if (M->anglesAreValid()) {
      result->insert(isoSig);
      std::cout << isoSig << "\n" << M->reginaIsoSig() << "\n";
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

	// Only copy the triangulation if the correct face f2 is open in the other tetrahedron.
	int f2 = f;
	if (f == 1) {
	  f2 = 2;
	}
	else if (f == 2) {
	  f2 = 1;
	}
	
	if (M->tetrahedron(i)->isOpen(f2)) {
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
