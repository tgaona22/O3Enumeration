#include "O3Triangulation.h"
#include "O3Tetrahedron.h"

#include <algorithm>
#include <set>

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

// Construct a triangulation from Mark's signature.
// We assume that O3IsoSig is a string of length 4*n containing only digits from 0 to n-1. 
O3Triangulation::O3Triangulation(const std::vector<int>& isoSig)
{
  if (isoSig.size() % 4 == 0) {
    int n = isoSig.size() / 4;

    for (int i = 0; i < n; i++) {
      tets.push_back(new O3Tetrahedron(&trig));
    }

    for (int j = 0; j < n; j++) {
      // For each tetrahedron, iterate over its faces and glue according to the O3IsoSig.
      // Namely, face f of tetrahedron j gets glued to tetrahedron whose index is O3IsoSig[4j + f].
      O3Tetrahedron *T = tets.at(j);
      for (int f = 0; f <= 3; f++) {
	int t = isoSig.at(4*j + f);
	T->join(f, tets.at(t));
      }
    }
  }
}

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
  std::vector<std::pair<int, int>> angles = getDihedralAngles();
  for (auto pair = angles.begin(); pair != angles.end(); pair++) {
    std::cout << "(" << (*pair).first << "," << (*pair).second << "), ";
  }
  std::cout << "\n";
}

bool O3Triangulation::anglesAreValid()
{
  std::vector<std::pair<int, int>> angles = getDihedralAngles();
  for (auto pair = angles.begin(); pair != angles.end(); pair++) {
    if (std::find(admissible_angles.begin(), admissible_angles.end(), *pair) == admissible_angles.end()) {
      return false;
    }
  }
  return true;
}

// Returns true iff M and this have the same set of edge labels and degrees.
bool O3Triangulation::compareDihedralAngles(O3Triangulation *M)
{
  std::vector<std::pair<int,int>> angles = getDihedralAngles();
  return std::is_permutation(angles.begin(), angles.end(), M->getDihedralAngles().begin());
}

std::vector<std::pair<int, int>> O3Triangulation::getDihedralAngles()
{
  std::vector<std::pair<int, int>> angles;
  for (auto edge = trig.edges().begin(); edge != trig.edges().end(); edge++) {
    regina::Perm<4> p = (*edge)->embedding(0).vertices();
    int degree = (*edge)->degree();
    angles.push_back(std::make_pair(edge_labels[p[0]][p[1]], degree));
  }
  return angles;
}

// Returns a list of 4n integers between 0 and n inclusive
// such the kth element of the list, which can be written as k = 4j + r
// is the index of the tetrahedron glued to face r of tetrahedron j.
// If this face is open, the kth element is labeled n (which is not the index of any tetrahedron,
// since they are indexed from 0 to n-1.)
std::vector<int> O3Triangulation::destinationSequence()
{
  std::vector<int> destSeq;
  // Iterate over all tetrahedra
  const int n = size();
  for (int i = 0; i < n; i++) {
    O3Tetrahedron *T = tetrahedron(i);
    // Iterate over all faces of T
    for (int f = 0; f <= 3; f++) {
      // Find the index of the O3tet glued to T at face f.
      int adj = T->adjacentSimplex(f);
      if (adj != -1) {
	destSeq.push_back(adj);
      }
      else {
	destSeq.push_back(n);
      }
    }
  }
  return destSeq;
}

/* Given a triangulation with N tetrahedra, construct the N destination sequences
   associated to the N canonical relabelings and return the lexicographically minimal one. */
std::vector<int> O3Triangulation::minimalDestinationSequence()
{
  std::vector<int> destSeq = destinationSequence();
  std::vector<int> minSeq = destSeq;
  
  const int n = size();
  
  if (n == 1) {
    return minSeq;
  }
    
  for (int i = 0; i < n; i++) {
    std::vector<int> L = {i};

    // Initialize D0 to contain T(i,0), T(i,1), T(i,2), T(i,3) from the destination sequence.
    std::vector<int> D0;
    for (int j = 0; j <= 3; j++) {
      D0.push_back(destSeq.at(4*i + j));
    }

    while (L.size() < n) {

      // Choose the first k in D0 such that k != n and k is not in L.
      int index = 0;
      int k = D0.at(index);
      while (k == n || std::find(L.begin(), L.end(), k) != L.end()) {
	index++;
	k = D0.at(index);
      }

      L.push_back(k);
      // Add T(k,0), T(k,1), T(k,2), T(k,3) to D0
      for (int j = 0; j <= 3; j++) {
	D0.push_back(destSeq.at(4*k + j));
      }
    }
    // When the above terminates, D0 should be a list of length 4n and L should be represent a permutation in S_n.

    // Get the permutation sigma in S_n, so sigma(j) = L(j).

    std::vector<int> sigmaInverse;
    for (int k = 0; k < n; k++) {
      // sigmaInverse(k) = index j such that k = l_j
      int j = 0;
      while (k != L.at(j)) {
	j++;
      }
      sigmaInverse.push_back(j);
    }
    // sigmaInverse(n) = n by definition
    sigmaInverse.push_back(n);

    // DI is the destination sequence for the relabeling determined by L.
    std::vector<int> DI;
    for (int k = 0; k < 4*n; k++) {
      DI.push_back(sigmaInverse.at(D0.at(k)));
    }

    // Compare DI to the current minSeq
    if (std::lexicographical_compare(DI.begin(), DI.end(), minSeq.begin(), minSeq.end())) {
      // DI < minSeq
      minSeq = DI;
    }
  }
  return minSeq;
}

void O3Triangulation::computeCuspCrossSections()
{
  // First, we need to sort ideal vertices into equivalence classes.
  // I'll store this as a list of sets of tetrahedra indices, where tetrahedra in the
  // same list have their ideal vertices identified.
  //std::vector<std::set<int>> cusps;
  std::vector<std::vector<int>> cusps;

  // Add the indices of all tetrahedra to a list of "unexamined" tetrahedra.
  std::vector<int> unexamined;
  for (int i = 0; i < size(); i++) {
    unexamined.push_back(i);
  }

  while (!unexamined.empty()) {
    std::vector<int> cusp;
    // Pop the first element off the list of unexamined and add it to a new cusp.
    cusp.push_back(unexamined.at(0));

    int cusp_iter = 0;
    while (cusp_iter < cusp.size()) {
      O3Tetrahedron *T = tetrahedron(cusp.at(cusp_iter));
      // remove the current tetrahedron from the list of unexamined.
      // something new I learned about c++ - remove doesn't actually delete anything...
      // it just moves them to the back of the array and gives you an iterator.
      unexamined.erase(std::remove(unexamined.begin(), unexamined.end(), cusp.at(cusp_iter)), unexamined.end());

      // look at what is glued to faces f0, f1, and e of this tetrahedron
      // add those indices to this cusp, if not already present.
      for (int f = 1; f <= 3; f++) {
	int adj_index = T->adjacentSimplex(f);
	auto found_iter = std::find(cusp.begin(), cusp.end(), adj_index);
	if (found_iter == cusp.end()) {
	  cusp.push_back(adj_index);
	}
      }
      // look at the next tetrahedron in the cusp.
      cusp_iter++;
    }

    cusps.push_back(cusp);
  }
  
  class Triangle {
  public:
    regina::Simplex<2> *t0, *t1;
    int tet_index;

    Triangle(regina::Simplex<2> *t0, regina::Simplex<2> *t1, int i) : t0(t0), t1(t1), tet_index(i) {}
    // returns the angle as a multiple of pi/6.
    int angle(int v)
    {
      if (v == 0) {
	return 3;
      }
      else if (v == 1) {
	return 2;
      }
      else {
	return 1;
      }
    }

    // 0 = f0, 1 = f1, 2 = e.
    bool isOpen(int f) {
      if (f == 0) {
	return !t0->adjacentSimplex(0);
      }
      else if (f == 1) {
	return !t1->adjacentSimplex(0);
      }
      else { // (f == 2)
	return !t0->adjacentSimplex(2);
      }
    }
  };
  
  std::vector<regina::Triangulation<2>> crossSections;
  // Construct a triangulation for the cross section of each cusp.
  int cusp_index = 0;
  for (auto cusp = cusps.begin(); cusp != cusps.end(); cusp++) {

    regina::Triangulation<2> crossSection;
    std::vector<Triangle> triangles;

    // Create a triangle for each tetrahedron.
    int tri_index = 0;
    regina::Perm<3> id;
    for (auto t = cusp->begin(); t != cusp->end(); t++) {
      O3Tetrahedron *T = tetrahedron(*t);

      Triangle tri(crossSection.newSimplex(), crossSection.newSimplex(), T->index());

      //join these two 2,3,6 triangles along an edge so they both make a (3,3,3) triangle.
      tri.t0->join(1, tri.t1, id);

      if (T->eIdentified()) {
	// join face e to itself
	tri.t0->join(2, tri.t1, id);
      }
      if (T->fIdentified()) {
	// join f0 to f1
	tri.t0->join(0, tri.t1, id);
      }
	
      T->setTriangleIndex(tri_index);
      tri_index = tri_index + 1;
      triangles.push_back(tri);
    }

    // For each triangle, make the appropriate gluings.
    for (int i = 0; i < triangles.size(); i++) {
      Triangle tri1 = triangles[i];
      O3Tetrahedron *T1 = tetrahedron(tri1.tet_index);

      for (int f = 1; f <= 3; f++) {
	// Only make a gluing if face f-1 is open in t1.
	if (tri1.isOpen(f-1)) {
	  
	  O3Tetrahedron *T2 = tetrahedron(T1->adjacentSimplex(f));
	  Triangle tri2 = triangles[T2->getTriangleIndex()];

	  // join f0 in tri1 to f1 in tri2
	  if (f == 1 && !tri2.t1->adjacentSimplex(0)) {
	    tri1.t0->join(0, tri2.t1, id);
	  }
	  // join f1 in tri1 to f0 in tri2
	  if (f == 2 && !tri2.t0->adjacentSimplex(0)) {
	    tri1.t1->join(0, tri2.t0, id);
	  }
	  // join e in tri1 to e in tri2
	  if (f == 3 && !tri2.t1->adjacentSimplex(2) && !tri2.t0->adjacentSimplex(2)) {
	    tri1.t0->join(2, tri2.t1, id);
	    tri1.t1->join(2, tri2.t0, id);
	  }
	}
	
      }
    }

    crossSections.push_back(crossSection);
    // After the above loop terminates, crossSection is a 2D triangulation which has been glued up according to
    // how the original 3D triangulation is glued up.

    // Next step is to determine all the cone points and their cone angles. 

    std::vector<int> conePoints;
    const regina::FaceList<2, 0>& vertices = crossSection.faces<0>();

    for (int i = 0; i < vertices.size(); i++) {
      regina::Face<2, 0> *vertex = vertices[i];
      // The resulting angle around the vertex will be degree * pi/6
      int degree = 0;
      // For each embedding of the vertex
      for (auto emb = vertex->begin(); emb != vertex->end(); emb++) {
	// find its triangle
	regina::Simplex<2> *t = emb->simplex();
	Triangle tri = *std::find_if(triangles.begin(), triangles.end(),
				     [t](Triangle Tri) -> bool {
				       return (Tri.t0 == t || Tri.t1 == t);
				     });
	// find the angle around the vertex
	degree = degree + tri.angle(emb->face());
      }

      // the degree should be either 2, 4, 6, or 12. The first three correspond to cone points of order 6, 3, and 2.
      if (degree <= 6) {
	int order = 12/degree;
	conePoints.push_back(order);
      }
    }
    
    std::cout << "\tCusp #" << cusp_index << ": ";
    for (auto entry = cusp->begin(); entry != cusp->end(); entry++) {
      std::cout << *entry << " ";
    }
    std::cout << "\n";

    std::cout << "\tCone points: ";

    std::sort(conePoints.begin(), conePoints.end());
    for (int i = 0; i < conePoints.size(); i++) {
      std::cout << conePoints[i] << " ";
    }
    std::cout << "\n";

    cusp_index++;
    
  }

  
}

bool O3Triangulation::checkClosedEdges()
{
  std::vector<regina::Face<3,1>*> closedEdges;
  const regina::FaceList<3,1>& edges = trig.faces<1>();
  // Iterate through the edges of the triangulation.
  for (int i = 0; i < edges.size(); i++) {

    bool isClosed = true;
    regina::Face<3,1> *edge = edges[i];

    // Iterate over all embeddings of the edge, checking that the faces adjacent to the edge
    // are glued. If all are, the edge is closed. If one is not, the edge is open.
    auto embedding = edge->begin();
    while (isClosed && embedding != edge->end()) {
      const regina::Simplex<3>* tet = embedding->simplex();
      const regina::Perm<4>& verts = embedding->vertices();
      // The faces adjacent to this edge are verts[2] and verts[3]. 
      if (not tet->adjacentSimplex(verts[2]) || not tet->adjacentSimplex(verts[3])) {
	isClosed = false;
      }
      embedding++;
    }

    const regina::Perm<4>& p = edge->embedding(0).vertices();
    int degree = edge->degree();
    // If the edge is closed, check that its angle is permissible.
    if (isClosed) {
      std::pair<int, int> angle = std::make_pair(edge_labels[p[0]][p[1]], degree);
      if (std::find(admissible_angles.begin(), admissible_angles.end(), angle) == admissible_angles.end()) {
	return false;
      }
    }
    else {
      // If the edge is open, check that its degree does not exceed the maximum possible for that edge label.
      if (degree > 2*edge_labels[p[0]][p[1]]) {
	return false;
      }	
    }

  }
  return true;
}
