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
O3Triangulation::O3Triangulation(std::string O3IsoSig)
{
  if (O3IsoSig.length() % 4 == 0) {
    int n = O3IsoSig.length() / 4;

    for (int i = 0; i < n; i++) {
      tets.push_back(new O3Tetrahedron(&trig));
    }

    for (int j = 0; j < n; j++) {
      // For each tetrahedron, iterate over its faces and glue according to the O3IsoSig.
      // Namely, face f of tetrahedron j gets glued to tetrahedron whose index is O3IsoSig[4j + f].
      O3Tetrahedron *T = tets.at(j);
      for (int f = 0; f <= 3; f++) {
	char t = O3IsoSig.at(4*j + f);
	T->join(f, tets.at(std::atoi(&t)));
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

// Helper function for sorting destination sequences
bool compareSequences(std::vector<int> a, std::vector<int> b)
{
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

// See Mark's algorithm
std::vector<int> O3Triangulation::minimalDestinationSequence()
{
  std::vector<int> destSeq = destinationSequence();
  std::vector<std::vector<int>> destSeqs;

  const int n = size();
  
  if (n == 1) {
    return destSeq;
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
    destSeqs.push_back(DI);
  }

  // Now we return the lexicographically minimal element of the two destination sequences.
  std::vector<int> isoSig = *std::min_element(destSeqs.begin(), destSeqs.end(), compareSequences);
  return isoSig;
}

std::string O3Triangulation::O3isoSig()
{
  std::vector<int> vec = minimalDestinationSequence();
  std::string sig = "";
  for (auto iter = vec.begin(); iter != vec.end(); iter++) {
    sig += std::to_string(*iter);
  }
  return sig;
}

void O3Triangulation::computeCuspCrossSections()
{
  // First, we need to sort ideal vertices into equivalence classes.
  // I'll store this as a list of sets of tetrahedra indices, where tetrahedra in the
  // same list have their ideal vertices identified.
  std::vector<std::set<int>> cusps;

  for (int t = 0; t < size(); t++) {
    // is there a set containing this index?
    // if not, create one.
    // if there is, we will be adding to that set.

    bool addToList = true;
    std::vector<std::set<int>>::iterator cusp;
    for (auto set = cusps.begin(); set != cusps.end(); set++) {
      if (set->find(t) != set->end()) {
	cusp = set;
	addToList = false;
      }
    }
    if (addToList) {
      std::set<int> set;
      cusps.push_back(set);
      cusp = cusps.end() - 1;
    }

    // for tetrahedron t, add to the class indices of the tetrahedra glued to faces e, f0, f1 of tet t.
    O3Tetrahedron *T = tetrahedron(t);
    for (int f = 1; f <= 3; f++) {
      cusp->insert(T->adjacentSimplex(f));
    }

  }

  // Testing purposes, print # of cusps and list of tetrahedra in each equivalence class.
  std::cout << "Number of cusps: " << cusps.size() << "\n";
  for (auto set = cusps.begin(); set != cusps.end(); set++) {
    std::cout << "{ ";
    for (auto elt = set->begin(); elt != set->end(); elt++) {
      std::cout << *elt << " ";
    }
    std::cout << "}\n";
  }

  struct Triangle {
    regina::Simplex<2> *triangle;
    int tet_index;
  };
  
  std::vector<regina::Triangulation<2>> crossSections;
  // Construct a triangulation for the cross section of each cusp.
  for (auto cusp = cusps.begin(); cusp != cusps.end(); cusp++) {

    regina::Triangulation<2> crossSection;
    std::vector<Triangle> triangles;

    // Create a triangle for each tetrahedron.
    int tri_index = 0;
    for (auto t = cusp->begin(); t != cusp->end(); t++) {
      O3Tetrahedron *T = tetrahedron(*t);
      Triangle Tri;
      Tri.triangle = crossSection.newSimplex();
      Tri.tet_index = T->index();
      T->setTriangleIndex(tri_index);
      tri_index = tri_index + 1;
      triangles.push_back(Tri);
    }

    // For each triangle, make the appropriate gluings.
    // the gluing map in each case is (0,1,2) -> (1,0,2).
    regina::Perm<3> gluing(0, 1);
    for (int i = 0; i < triangles.size(); i++) {
      Triangle t1 = triangles[i];
      O3Tetrahedron *T1 = tetrahedron(t1.tet_index);

      for (int f = 1; f <= 3; f++) {
	O3Tetrahedron *T2 = tetrahedron(T1->adjacentSimplex(f));
	Triangle t2 = triangles[T2->getTriangleIndex()];
	// Don't do anything in the case that T1 has e glued to itself.
	if (!(f == 3 && T2 == T1)) {
	  t1.triangle->join(f-1, t2.triangle, gluing);
	}
      }
    }

    crossSections.push_back(crossSection);
    // After the above loop terminates, crossSection is a 2D triangulation which has been glued up according to
    // how the original 3D triangulation is glued up, except for the case when a tetrahedron has face e glued to itself.

    // Next step is to determine all the cone points and their cone angles. 

    // Iterate over the vertices of the triangulation. The cone angle of a vertex of degree k is
    // k * 2Pi/6. The possibilities for k should be 1, 2, 3, or 6.
    std::vector<int> conePoints;
    const regina::FaceList<2, 0>& vertices = crossSection.faces<0>();
    std::vector<regina::Face<2,0>*> alreadySeenVertices;
    for (int i = 0; i < vertices.size(); i++) {
      regina::Face<2, 0> *vertex = vertices[i];
      // if we haven't encountered this vertex yet...
      if (std::find(alreadySeenVertices.begin(), alreadySeenVertices.end(), vertex) == alreadySeenVertices.end()) {
	// we have seen the vertex
	alreadySeenVertices.push_back(vertex);

	int degree = vertex->degree();
	
	// For each embedding of the vertex
	for (auto emb = vertex->begin(); emb != vertex->end(); emb++) {

	  // find its triangle, so we can access the associated tetrahedron.
	  regina::Simplex<2> *t = emb->simplex();
	  Triangle tri = *std::find_if(triangles.begin(), triangles.end(),
				    [t](Triangle Tri) -> bool {
				      return (Tri.triangle == t);
				    });
	  // if that tetrahedron has face e glued to itself
	  O3Tetrahedron *T = tetrahedron(tri.tet_index);
	  if (T->adjacentSimplex(O3Tetrahedron::e) == T->index()) {
	    // the vertex is either f0 or f1 (0 or 1)
	    regina::Face<2,0> *toIdentify;
	    if (emb->face() == 0) {
	      toIdentify = tri.triangle->face<0>(1);
	    }
	    else {
	      toIdentify = tri.triangle->face<0>(0);
	    }
	    // add the number of embeddings of this vertex which in the real triangulation
	    // would be identified with the current vertex. This unfortunate kludge is necessary
	    // since I can't glue an edge to itself in a regina triangulation.
	    if (toIdentify != vertex) {
	      degree = degree + toIdentify->degree();
	      alreadySeenVertices.push_back(toIdentify);
	    }
	  }
	}

	
	conePoints.push_back(degree);
	
      }
    }
    // Also, we produce a cone point of order 2 whenever face e is identified to itself in a tetrahedron.
    for (int i = 0; i < size(); i++) {
      O3Tetrahedron *T = tetrahedron(i);
      if (T->adjacentSimplex(O3Tetrahedron::e) == i) {
	//conePoints.push_back(2);
	conePoints.push_back(3);
      }
    }

    std::cout << "Cone points: ";
    for (int i = 0; i < conePoints.size(); i++) {
      std::cout << conePoints[i] << " ";
    }
    std::cout << "\n";
    
  }

  
}
