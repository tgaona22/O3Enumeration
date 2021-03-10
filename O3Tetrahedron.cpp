#include "O3Tetrahedron.h"
#include "O3Triangulation.h"

O3Tetrahedron::O3Tetrahedron(regina::Triangulation<3>* trig) : trig(trig)
{
  s0 = trig->newSimplex();
  s1 = trig->newSimplex();
  s0->join(1, s1, glue);
}

// used when copying a triangulation
O3Tetrahedron::O3Tetrahedron(regina::Triangulation<3>* trig, int s0Index, int s1Index) : trig(trig)
{
  s0 = trig->simplex(s0Index);
  s1 = trig->simplex(s1Index);
}

int O3Tetrahedron::underlyingIndex(int tet)
{
  if (tet == 0) {
    return s0->index();
  }
  if (tet == 1) {
    return s1->index();
  }
  return -1;
}

void O3Tetrahedron::join(int myFace, O3Tetrahedron* you)
{
  // For a given face (v, f0, f1, e) of an O3Tet and another O3Tet,
  // there is a unique way to pair faces.
  // This function does nothing if the corresponding face in you is not open.
  if (myFace == v) {
    // v always glues to v.
    // If corresponding face of you is not open, do nothing.
    if (you->isOpen(v)) {
      if (you == this) {
	// Glue v to itself.
	s0->join(0, s1, glue);
      }
      else {
	// Glue v to v in the other tetrahedron
	s0->join(0, you->s1, glue);
	(you->s0)->join(0, s1, glue);
      }
    }
  }
  else if (myFace == f0) {
    // f0 always glues to f1.
    if (you->isOpen(f1)) {
      s1->join(2, you->s0, glue);
    }
  }
  else if (myFace == f1) {
    // f1 always glues to f0.
    if (you->isOpen(f0)) {
      s0->join(2, you->s1, glue);
    }
  }
  else if (myFace == e) {
    // e always glues to e.
    if (you->isOpen(e)) {
      if (you == this) {
	// Glue e to itself.
	s0->join(3, s1, glue);
      }
      else {
	// Glue e to e in the other tetrahedron
	s0->join(3, you->s1, glue);
	(you->s0)->join(3, s1, glue);
      }
    }
  }
}

bool O3Tetrahedron::isOpen(int f)
{
  std::pair<regina::Simplex<3>*, int> pair = face(f);
  if (pair.first->adjacentSimplex(pair.second)) {
    return false;
  }
  return true;
}

std::pair<regina::Simplex<3>*, int> O3Tetrahedron::face(int f)
{
  switch(f) {
  case v:
    return std::make_pair(s0, 0);
  case f0:
    return std::make_pair(s1, 2);
  case f1:
    return std::make_pair(s0, 2);
  case e:
    return std::make_pair(s0, 3);
  }
  return std::make_pair(nullptr, 0);
}

//std::pair<int, int> O3Tetrahedron::dihedralAngle(int v1, int v2)

