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

    // I would like to use regina's Perm. But it is implemented as a template class, depending on an int n
    // which must be determined at compile time. Here our n is the size of the triangulation which is known not until runtime.
    //regina::Perm sigma(L.data());
    //regina::Perm sigmaInverse = sigma.inverse();
    
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

    /*
    std::cout << "L: ";
    for (auto it = L.begin(); it != L.end(); it++) {
      std::cout << *it << ", ";
    }
    std::cout << "\n";


    std::cout << "SigmaInverse: ";
    for (auto it = sigmaInverse.begin(); it != sigmaInverse.end(); it++) {
      std::cout << *it << ", ";
    }
    std::cout << "\n";
    */

    /*
    std::cout << "D0: ";
    for (auto it = D0.begin(); it != D0.end(); it++) {
      std::cout << *it << ", ";
    }
    std::cout << "\n";
    */
      
    
    // DI is the destination sequence for the relabeling determined by L.
    std::vector<int> DI;
    //    std::cout << "DI: ";
    for (int k = 0; k < 4*n; k++) {
      DI.push_back(sigmaInverse.at(D0.at(k)));
      //std::cout << sigmaInverse.at(D0.at(k)) << ", ";
    }
    //std::cout << "\n";
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
