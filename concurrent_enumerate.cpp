#include "O3Triangulation.h"
#include "O3Tetrahedron.h"
#include "threadpool.h"
#include "threadSafeSet.h"

#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <chrono>
#include <memory>

using namespace std::chrono;
template <
    class result_t   = std::chrono::milliseconds,
    class clock_t    = std::chrono::steady_clock,
    class duration_t = std::chrono::milliseconds
>
auto since(std::chrono::time_point<clock_t, duration_t> const& start)
{
    return std::chrono::duration_cast<result_t>(clock_t::now() - start);
}


void recurse(std::shared_ptr<O3Triangulation> M, int max_tets, int level);
void printIsoSig(const std::vector<int> &vec);

// Define the global variables already_seen, result, and threadpool.
ThreadSafeSet<std::vector<int>> already_seen;
ThreadSafeSet<std::vector<int>> result;
Threadpool threadpool;

int main(int argc, char **argv)
{
  /* ENUMERATION */

  if (argc < 1 || argc > 2) {
    std::cout << "Enter the max number of tetrahedra.\n";
    return -1;
  }
  int N = std::atoi(argv[1]);

  auto start = std::chrono::steady_clock::now();
  
  std::shared_ptr<O3Triangulation> O = std::make_shared<O3Triangulation>();
  O->newTetrahedron();

  {
    Threadpool::RunningContext runningContext(threadpool);

    recurse(O, N, 0);
  }

  std::vector<std::set<std::vector<int>>> structured_result(N);
  for (auto it = result.get_set_unsafe().begin(); it != result.get_set_unsafe().end(); it++) {
    int sz = it->size() / 4;
    structured_result.at(sz-1).insert(*it);
  }

  for (int i = 0; i < N; i++) {
    //std::cout << result[i].size() << " triangulations with " << (i+1) << " tetrahedra:\n";
    for (auto it = structured_result[i].begin(); it != structured_result[i].end(); it++) {
      printIsoSig(*it);
    }
  }

  std::cout << "Elapsed(ms)=" << since(start).count() << std::endl;

}

void recurse(std::shared_ptr<O3Triangulation> M, int max_tets, int level)
{
  // If there are closed edges with invalid angles, throw away this triangulation.
  if (not M->checkClosedEdges()) {
    return;
  }

  const std::vector<int> isoSig = M->O3isoSig();

  // Attempt to insert the isoSig to already_seen.
  if (not already_seen.insert(isoSig).second) {
    // If we get here, the isoSig was already in already_seen, so we've seen this triangulation before.
    return;
  }
      
  // Get the open faces in the triangulation.
  std::vector<std::pair<int, int>> openFaces = M->getOpenFaces();

  
  // If there are no open faces, we're done with this triangulation.
  if (openFaces.empty()) {
    // If we arrive here, we are guaranteed that all dihedral angles are valid.
    result.insert(isoSig);
    return;
  }

  int fixedTetIndex = openFaces.front().first;
  
  if (M->size() < max_tets) {
    // Try all four ways of adding a new tetrahedron to the chosen tetrahedron.
    for (int f = 0; f <= 3; f++) {
      // If face f is open in the fixed tetrahedron
      if (M->tetrahedron(fixedTetIndex)->isOpen(f)) {
	// Copy the triangulation.
	std::shared_ptr<O3Triangulation> N = std::make_shared<O3Triangulation>(*M);

	// Create a new tetrahedron T2 and glue it to face f of an existing tetrahedron T1.
	// We choose T1 arbitrarily to be the first tetrahedron in the list of open face pairs.
	O3Tetrahedron *T1 = N->tetrahedron(fixedTetIndex);
	O3Tetrahedron *T2 = N->newTetrahedron();
	T1->join(f, T2);

	// Maybe the problem is we're passing a pointer to the new thread - and this thread
	// is continuing on and the pointer is destroyed when we go out of scope of this block.
	// To fix this I'll first try to just pass copies of regina triangulations
	// if that fails, I will rewrite a lightweight triangulation class.
	threadpool.scheduleOrExecute(std::bind(&recurse, N, max_tets, level+1));
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
	  std::shared_ptr<O3Triangulation> N = std::make_shared<O3Triangulation>(*M);
	  // Get the fixed face.
	  O3Tetrahedron *T1 = N->tetrahedron(fixedTetIndex);
	  O3Tetrahedron *T2 = N->tetrahedron(i);
	  T1->join(f, T2);

	  threadpool.scheduleOrExecute(std::bind(&recurse, N, max_tets, level+1));
	}
      }
    }
  }
}

void printIsoSig(const std::vector<int> &vec)
{
  int count = 0;
  for (int i = 0; i < vec.size(); i++) {
    std::cout << vec[i];
    count = count + 1;
    if (count == 4) {
      std::cout << " ";
      count = 0;
    }
    else {
      std::cout << ",";
    }
  }
  std::cout << "\n";
}
