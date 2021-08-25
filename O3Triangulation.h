#ifndef O3TRIANGULATION_H
#define O3TRIANGULATION_H

#include <vector>
#include <string>
#include <engine.h>
#include <triangulation/dim3.h>
#include <triangulation/dim2.h>

class O3Tetrahedron;

class O3Triangulation {
 public:
  O3Triangulation();
  O3Triangulation(const O3Triangulation& M);
  O3Triangulation(const std::vector<int>& isoSig);
  //O3Triangulation& operator=(const O3Triangulation& M);
  ~O3Triangulation();

  O3Tetrahedron* newTetrahedron();
  O3Tetrahedron* tetrahedron(int index);
  std::vector<std::pair<int,int>> getOpenFaces();
  bool anglesAreValid();
  std::vector<std::pair<int, int>> getDihedralAngles();
  bool compareDihedralAngles(O3Triangulation* M);
  void printDihedralAngles();

  std::vector<int> O3isoSig() { return minimalDestinationSequence(); }
  std::string O3isoSigString();
  std::string reginaIsoSig() { return trig.isoSig(); }
  int size() { return trig.size()/2; }
  int true_size() { return trig.size(); }

  bool checkClosedEdges();

  // The following two methods will produce undefined behavior for partially closed triangulations (ones that still have open faces)
  void computeCuspCrossSections();
  //void computeSingularLocus();
  void printSingularLocus(std::string directory);

  std::vector<int> destinationSequence();
 private:
  const int edge_labels[4][4] = { {0, 3, 2, 6},
				  {3, 0, 3, 2},
				  {2, 3, 0, 2},
				  {6, 2, 2, 0} };
  const std::vector<std::pair<int, int>> admissible_angles = { std::pair<int,int>(3,1),
							       std::pair<int,int>(2,1),
							       std::pair<int,int>(3,2),
							       std::pair<int,int>(2,2),
							       std::pair<int,int>(3,3),
							       std::pair<int,int>(2,4),
							       std::pair<int,int>(3,6),
							       std::pair<int,int>(6,2),
							       std::pair<int,int>(6,3),
							       std::pair<int,int>(6,4),
							       std::pair<int,int>(6,6),
							       std::pair<int,int>(6,12) };
  regina::Triangulation<3> trig;
  std::vector<O3Tetrahedron*> tets;

  std::vector<std::vector<std::vector<int>>> singular_graph;
  void computeSingularLocus();
  void print_graph();
  int get_nonsingular_vertex();
  std::vector<int> adjacent_vertices(int vertex);
  int valence(int vertex);
  bool is_loop(int vertex);
  

  std::vector<int> minimalDestinationSequence();
};

#endif
