#ifndef O3TETRAHEDRON_H
#define O3TETRAHEDRON_H

#include <engine.h>
#include <triangulation/dim3.h>

class O3Tetrahedron {
 public:
  O3Tetrahedron(regina::Triangulation<3>* trig);
  O3Tetrahedron(regina::Triangulation<3>* trig, int s0Index, int s1Index);

  void join(int myFace, O3Tetrahedron* you);
  bool isOpen(int f);
  bool eIdentified() { return adjacentSimplex(O3Tetrahedron::e) == index(); }
  //std::pair<int, int> dihedralAngle(int v1, int v2);

  std::vector<std::pair<int, int>> dihedralAngles();
  int adjacentSimplex(int f);

  int underlyingIndex(int tet);
  int index() { return s0->index()/2; }

  void setTriangleIndex(int index) { triangle_index = index; }
  int getTriangleIndex() { return triangle_index; }

  static const int v = 0;
  static const int f0 = 1;
  static const int f1 = 2;
  static const int e = 3;

 private:

  regina::Perm<4> glue; // the identity permutation.
  regina::Simplex<3> *s0, *s1;
  regina::Triangulation<3>* trig;
  //std::vector<std::pair<regina::Edge<3>*, int>> edges;
  //std::unordered_map<regina::Edge<3>*, int> edge_labels;

  std::pair<regina::Simplex<3>*, int> face(int f);
  //  regina::Edge<3> edge(int v1, int v2);

  int triangle_index;
  
};

#endif
