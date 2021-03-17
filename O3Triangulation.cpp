#include "O3Triangulation.h"
#include "O3Tetrahedron.h"

#include <algorithm>

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
  for (int i = 0; i < n; i++) {
    std::vector<int> L = {i};

    // Initialize D0 to contain T(i,0), T(i,1), T(i,2), T(i,3) from the destination sequence.
    std::vector<int> D0;
    for (int j = 0; j <= 3; j++) {
      D0.push_back(destSeq.at(4*i + j));
    }

    // while size(L) < n... let k = size(L)
    // at termination of each iteration size(L) increases by one.
    for (int k = 1; k < n; k++) {

      // Choose the first index l in D0 such that l != n and l is not in L.
      int l = 0;
      while (l == n || std::find(L.begin(), L.end(), l) != L.end()) {
	l++;
      }

      L.push_back(l);
      // Add T(l,0), T(l,1), T(l,2), T(l,3) to D0
      for (int j = 0; j <= 3; j++) {
	D0.push_back(destSeq.at(4*l + j));
      }
    }
    // When the above terminates, D0 should be a list of length 4n and L should be represent a permutation in S_n.

    // Get the permutation sigma in S_n, so sigma(j) = L(j).

    // I would like to use regina's Perm. But it is implemented as a template class, depending on an int n
    // which must be determined at compile time. Here our n is the size of the triangulation which is not until runtime.
    //regina::Perm sigma(L.data());
    //regina::Perm sigmaInverse = sigma.inverse();
    
    std::vector<int> sigmaInverse;
    for (int k = 0; k < n; k++) {
      // sigmaInverse(x) = index j such that x = l_j
      int j = 0;
      while (k != L.at(j)) {
	j++;
      }
      sigmaInverse.push_back(j);
    }

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
