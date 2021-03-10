#include "O3Triangulation.h"
#include "O3Tetrahedron.h"

O3Triangulation::O3Triangulation()
{
  
}

O3Triangulation::O3Triangulation(const O3Triangulation& M) : trig(M.trig)
{
  //trig = M.trig; // copy the triangulation
  for (auto iter = M.tets.begin(); iter != M.tets.end(); iter++) {
    tets.push_back(new O3Tetrahedron(&trig, (*iter)->underlyingIndex(0), (*iter)->underlyingIndex(1)));; // make new copies of the O3Tetrahedra. 
  }
}

/*O3Triangulation& O3Triangulation::operator=(const O3Triangulation& M)
{
  if (&M != this) {
    for (auto iter = tets.begin(); iter != tets.end(); iter++) {
      delete *iter;
    }
    tets.clear();

    trig(M.trig); // reassign the triangulation
    for (auto iter = M.tets.begin(); iter != M.tets.end(); iter++) {
      tets.push_back(new O3Tetrahedron(&trig, (*iter)->underlyingIndex(0), (*iter)->underlyingIndex(1)));
    }
  }
  return *this;
  }*/

O3Triangulation::~O3Triangulation()
{
  for (auto iter = tets.begin(); iter != tets.end(); iter++) {
    delete *iter;
  }
}

O3Tetrahedron* O3Triangulation::newTetrahedron()
{
  O3Tetrahedron* tet = new O3Tetrahedron(&trig);
  tets.push_back(tet);
  return tet;
}

O3Tetrahedron* O3Triangulation::tetrahedron(int index)
{
  if (index >= 0 && index < tets.size()) {
    return tets.at(index);
  }
  return nullptr;
}

// Returns a list of pairs (i,j) such that
// tets[i]->face(j) will be open.
std::vector<std::pair<int,int>> O3Triangulation::getOpenFaces()
{
  std::vector<std::pair<int, int>> indices;
  for (int i = 0; i < tets.size(); i++) {
    for (int j = 0; j <= 3; j++) {
      if (tets.at(i)->isOpen(j)) {
	indices.push_back(std::make_pair(i,j));
      }
    }
  }
  return indices;
}

void O3Triangulation::printDihedralAngles()
{
  for (auto edge = trig.edges().begin(); edge != trig.edges().end(); edge++) {
    regina::Perm<4> p = (*edge)->embedding(0).vertices();
    int degree = (*edge)->degree();
    std::cout << "(" << edge_labels[p[0]][p[1]] << "," << degree << ")\n";
  }
}

bool O3Triangulation::anglesAreValid()
{
  for (auto edge = trig.edges().begin(); edge != trig.edges().end(); edge++) {
    regina::Perm<4> p = (*edge)->embedding(0).vertices();
    int degree = (*edge)->degree();
    //std::cout << "(" << edge_labels[p[0]][p[1]] << "," << degree << ")\n";
    if (std::find(admissible_angles.begin(), admissible_angles.end(), std::make_pair(edge_labels[p[0]][p[1]], degree)) == admissible_angles.end()) {
      return false;
    }
  }
  return true;
}
