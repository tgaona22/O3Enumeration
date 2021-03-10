#include <engine.h>
#include <census/census.h>
#include <triangulation/dim3.h>
#include <packet/container.h>

#include <cstring>

using namespace regina;

void printGluingInfo(Triangulation<3>* trig);

// The result of checking the number of tetrahedra adjacent to an edge
typedef enum {
	      EdgeValid,    // The number is less than 6 for an open edge or exactly 6
	      EdgeFixed,    // The number is 6 for an open edge that was closed up
	      EdgeInvalid   // The number is less than 6 for a closed edge or greater 6
} FixEdgeResult;

// Some handy typedefs
typedef std::deque<regina::EdgeEmbedding<3>>::const_iterator
EdgeEmbeddingIterator;

// An open face of triangulation encoded as index of the tetrahedron and index
// of the vertex opposite to that face.
typedef std::pair<int, int> TetIndexAndFace;

// Container of isomorphism signatures.
typedef std::set<std::string> IsoSigContainer;


// Given a closed edge, return a pair of empty edge embeddings.
//
// Given an open edge, return a pair of edge embeddings for this edge that
// correspond to the two open faces adjacent to this edge.
// The two open faces are lying on the tetrahedra of the two edge embeddings.
// For the first edge embedding, this faces is opposite to the vertex of
// the tetrahedron indexed by edgeEmbedding->getVertices()[2].
// Similar, for the second edge embedding and edgeEmbedding->getVertices()[3].

std::pair<EdgeEmbedding<3>, EdgeEmbedding<3>>
getOpenEdgeEmbeddingPair(const Edge<3> *edge)
{
  EdgeEmbedding<3> embed2;
  EdgeEmbedding<3> embed3;

  // Iterate through the edge embeddings
  for (EdgeEmbeddingIterator edgeEmbedding = edge->begin();
       edgeEmbedding != edge->end();
       ++edgeEmbedding) {
    
    const Simplex<3>* tet = edgeEmbedding->simplex();
    const Perm<4>& verts = edgeEmbedding->vertices();
    
    // No tet glued to verts[2], this is the first embedding
    if (not tet->adjacentSimplex(verts[2])) {
      embed2 = *edgeEmbedding;
    }
    // No tet glued to verts[3], this is the second embedding
    if (not tet->adjacentSimplex(verts[3])) {
      embed3 = *edgeEmbedding;
    }
  }

  return std::make_pair(embed2, embed3);
}

std::pair<EdgeEmbedding<3>, EdgeEmbedding<3>>
getGluingPair(const Edge<3> *edge) {

  // Iterate through the edge embeddings
  for (EdgeEmbeddingIterator firstEmbedding = edge->begin();
       firstEmbedding != edge->end();
       ++firstEmbedding) {
    
    const Simplex<3>* tet = firstEmbedding->simplex();
    const Perm<4>& verts = firstEmbedding->vertices();
    
    // No tet glued to verts[2], this is the first edge embedding
    if (not tet->adjacentSimplex(verts[2])) {
      // Loop through embeddings again, looking for a valid pairing with the first embedding.
      for (EdgeEmbeddingIterator secondEmbedding = edge->begin();
	   secondEmbedding != edge->end();
	   ++secondEmbedding) {

	const Simplex<3>* tet2 = secondEmbedding->simplex();
	const Perm<4>& verts2 = secondEmbedding->vertices();

	// If the second edge embedding has an open face
	if (not tet2->adjacentSimplex(verts2[3])) {
	  // Check that the tetrahedra are different or the tetrahedra are the same but the open faces are different.
	  if (tet != tet2 or verts[2] != verts2[3]) {
	    return std::make_pair(*firstEmbedding, *secondEmbedding);
	  }
	}
      }
    }
  }
  
  // If we get here, no valid pair was found and we return an empty pair of embeddings;
  EdgeEmbedding<3> empty;
  return std::make_pair(empty, empty);
}

// Take an edge embedding pair such as return from getOpenEdgeEmbeddingPair
// and glue the two open faces.
void
glueOpenEdgeEmbeddingPair(const std::pair<EdgeEmbedding<3>, EdgeEmbedding<3>>& pair)
{
  Simplex<3> *tet2   = pair.first.simplex();
  const Perm<4> &verts2 = pair.first.vertices();
  Simplex<3> *tet3   = pair.second.simplex();
  const Perm<4> &verts3 = pair.second.vertices();

  // The open faces for a pair returned from getOpenEdgeEmbeddingPair
  // are verts2[2] and verts3[3] respectively. This permutation sends
  // verts2[0] and verts2[1] to verts3[0] and verts3[1] to preserve the
  // edge and then exchange the two faces. See the documentation for
  // Simplex::join and Simplex::adjacentGluing for details on how
  // the permutation determines the face pairing.
  
  const Perm<4> perm(verts2[0], verts3[0],
		    verts2[1], verts3[1],
		    verts2[2], verts3[3],
		    verts2[3], verts3[2]);

  tet2->join(verts2[2], tet3, perm);
}

// Iterate through the edges of the triangulation.
// Check the number of edge embeddings (number of not necessarily distinct
// tetrahedra touching this edge) and whether the edge is valid.
// If there is an edge with the wrong number or an invalid edge,
// return EdgeInvalid.
// If an edge needed to be glued up, return EdgeFixed to indicate to the caller
// that the process needs to be restarted.
FixEdgeResult
fixOneEdge(Triangulation<3> *trig)
{
  // Iterate through the edges
  for (Triangulation<3>::EdgeIterator edge = trig->edges().begin();
       edge != trig->edges().end();
       ++edge) {

    // Return invalid if the orientations along the edge don't match.
    // This would result in a non-manifold because the link around the
    // center of the edge would be a RP2.
    if (not (*edge)->isValid()) {
      return EdgeInvalid;
    }
	
    // Determine whether the edge is open (i.e., unglued faces touch it,
    // the pair holds the corresponding edge embeddings) or closed.
    std::pair<EdgeEmbedding<3>, EdgeEmbedding<3>> embeddingPair =
      getOpenEdgeEmbeddingPair(*edge);

    // If the edge is closed but less than 6 edge embeddings,
    // return invalid.
    if ((*edge)->degree() < 6
	and not embeddingPair.first.simplex()) {
      return EdgeInvalid;
    }

    // If the edge has exactly 6 edge embeddings and is open,
    // close the edge. Indicate this to caller.
    if ((*edge)->degree() == 6
	and embeddingPair.first.simplex()) {
      // If we are given a bad pair, where tets and faces are the same, we need to find another pair.
      if (embeddingPair.first.simplex() == embeddingPair.second.simplex()
	  and embeddingPair.first.vertices()[2] == embeddingPair.second.vertices()[3]) {
	embeddingPair = getGluingPair(*edge);
	// We may not be able to find a good pair, in which case we return EdgeInvalid.
	if (not embeddingPair.first.simplex()) {
	  return EdgeInvalid;
	}
      }
      // If we get here, we have a valid edge embedding pair and can glue them. 
      glueOpenEdgeEmbeddingPair(embeddingPair);
      return EdgeFixed;
    }

    // If the edge has too many edge embedding, invalid.
    if ((*edge)->degree() > 6) {
      return EdgeInvalid;
    }
  }

  // All edges passed: valid!
  return EdgeValid;
}

// Glues up all edges that have 6 edge embeddings. It does so iteratively,
// i.e., if one gluing results in another edge having 6 edge embeddings, that
// edge will also be glued and so on.
// If the result is not valid, return invalid.
FixEdgeResult
fixEdges(Triangulation<3> *trig)
{
  while(true) {
    FixEdgeResult result = fixOneEdge(trig);
    if (result != EdgeFixed) {
      return result;
    }
  }

  // Should never reach here, needed to make compile.
  return EdgeInvalid;
}

// Get an open face of the triangulation.
// The algorithm is quicker if we pick one with an edge that is already 
// close to be glued up, i.e., that has a lot of edge embeddings.
TetIndexAndFace getOpenFace(Triangulation<3> *trig)
{
  // The highest number of embeddings for an edge encountered so far.
  int best_number_of_embeddings = 0;
  // The open face for which this occured
  TetIndexAndFace result;

  // Iterate through edges
  for (Triangulation<3>::EdgeIterator edge = trig->edges().begin();
       edge != trig->edges().end();
       ++edge) {

    // If this edge is better than what we saw so far.
    if ((*edge)->degree() > best_number_of_embeddings) {
	
      // Try to get one of the open faces
      EdgeEmbedding<3> edgeEmbedding =
	getOpenEdgeEmbeddingPair(*edge).first;

      // If this edge is open
      if (edgeEmbedding.simplex()) {
	// Set the new result
	best_number_of_embeddings = (*edge)->degree();
	  result = TetIndexAndFace(
				   edgeEmbedding.simplex()->index(),
				   edgeEmbedding.vertices()[2]);
      }
    }
  }
  
  return result;
}

// Just get all open faces of the triangulation
std::vector<TetIndexAndFace> getAllOpenFaces(Triangulation<3> *trig)
{
  std::vector<TetIndexAndFace> result;

  int i = 0; // Index of tetrahedron

  // Iterate through all tets
  for (Triangulation<3>::TetrahedronIterator tet = trig->simplices().begin();
       tet != trig->simplices().end();
       ++tet, ++i) {

    // Iterate through all faces
    for (int face = 0; face < 4; face++) {
      // Check that face is open, add to result
      if (not (*tet)->adjacentTetrahedron(face)) {
	result.push_back(TetIndexAndFace(i, face));
      }
    }
  }

  return result;
}

// The main function: it recurses to build all the triangulations
// M: the partial triangulation build so far, will be modified by gluing up
// the edges.
// max_tets: the maximal number of tetrahedra considered
// orientable: whether to return orientable or nonorientable triangulations
// already_seen: keeps track of isomorphism signature already handled
// result: result will be written in here
void recurse(Triangulation<3> *M, int max_tets, bool orientable,
	     IsoSigContainer *already_seen, IsoSigContainer *result)
{
  // Glue up edges that have 6 edge embeddings.
  // Check that result is valid, otherwise bail.
  if (fixEdges(M) == EdgeInvalid) {
    return;
  }

  // Get the isomorphism signature of the current triangulation
  const std::string isoSig = M->isoSig();

  // Bail if already handled earlier
  if (already_seen->find(isoSig) != already_seen->end()) {
    return;
  }
  // Mark as already handled
  already_seen->insert(isoSig);
    
  // Get all open faces
  const std::vector<TetIndexAndFace> tetIndexAndFaces2 = getAllOpenFaces(M);

  // If no open faces, we are done
  if (tetIndexAndFaces2.empty()) {
    // Add to results if orientability matches
    // (we don't even produce intermediate nonorientable results
    // when orientable is set)
    if (orientable or not M->isOrientable()) {
      result->insert(isoSig);
	printGluingInfo(M);
	}
    return;
  }
    
  // Get a good choice of a face
  const TetIndexAndFace tetIndexAndFace1 = getOpenFace(M);

  // If maximal number of tetrahedra has not be reached yet, glue
  // a new tetrahedron
  if (M->size() < max_tets) {
    // Make a copy of triangulation
    Triangulation<3> N(*M);
    // Odd permutation so that orientation of the two glued
    // tetrahedra matches
    static const Perm<4> perm(1,0,2,3);

    // Get the chosen face (in the copy)
    Simplex<3>* tet1 = N.simplex(tetIndexAndFace1.first);
    // Create a new tetrahedron
    Simplex<3>* tet2 = N.newSimplex();
    // Glue
    tet1->join(tetIndexAndFace1.second, tet2, perm);

    // And recurse
    recurse(&N, max_tets, orientable,
	    already_seen, result);
  }

  // Iterate through all open faces
  for (std::vector<TetIndexAndFace>::const_iterator
	 tetIndexAndFace2 = tetIndexAndFaces2.begin();
       tetIndexAndFace2 != tetIndexAndFaces2.end();
       ++tetIndexAndFace2) {
    // If this face is not the chosen face
    if (tetIndexAndFace1 != *tetIndexAndFace2) {
      // Try all permutations...
      for (int i = 0; i < 24; i++) {
	const Perm<4> perm(Perm<4>::Sn[i]);
	// Pick only permutations where the one open face is glued
	// to the other open face
	if (perm[tetIndexAndFace1.second] == tetIndexAndFace2->second) {
	  // If orientable, consider only odd permutations
	  if ((not orientable) or perm.sign() == -1) {
	    // Make a copy of the triangulation
	    Triangulation<3> N(*M);
	    // Get the chosen face (in the copy)
	    Simplex<3>* tet1 =
	      N.simplex(tetIndexAndFace1.first);
	    // Get this face (in the copy)
	    Simplex<3>* tet2 =
	      N.simplex(tetIndexAndFace2->first);
	    
	    // Glue
	    tet1->join(tetIndexAndFace1.second, tet2, perm);

	    // And recurse
	    recurse(&N, max_tets, orientable,
		    already_seen, result);
	  }
	}
      }
    }
  }
}

// Print usage
void printUsage()
{
  std::cerr << "Usage: genIsomoSigsOfTetrahedralTriangulations NUM_TETS [n|o]"
	    << std::endl;
  exit(1);
}

void printGluingInfo(Triangulation<3>* trig)
{
  std::vector<std::size_t> fVector = trig->fVector();
  std::cout << "Found an open edge with 6 embeddings. The triangulation has faces of the following dimensions:\n";
  for (int i=0; i<fVector.size(); i++) {
    std::cout << "\tDim " << i << ": " << fVector.at(i) << "\n";
  }
  std::cout << "The gluing data for each simplex is: \n";
  for (Triangulation<3>::SimplexIterator simps = trig->simplices().begin();
       simps != trig->simplices().end();
       simps++) {
    std::cout << "Tet " << (*simps)->index() << ":\n";
    for (int i = 0; i < 4; i++) {
      std::cout << "\tFace " << i;
      if ((*simps)->adjacentSimplex(i)) {
	std::cout << " adjacent to Tet " << (*simps)->adjacentSimplex(i)->index() << " via " << (*simps)->adjacentGluing(i).str() << "\n";
      }
      else {
	std::cout << " has no adjacent Tet.\n";
      }
    }
  }
}

int main(int argc, char** argv)
{
  if (argc < 2 or argc > 3) {
    printUsage();
  }

  // Get number of tets
  const int number_tets = std::atoi(argv[1]);

  // Fallback to orientable
  bool orientable = true;

  // Parse orientability
  if (argc == 3) {
    if (std::strcmp(argv[2], "n") == 0) {
      orientable = false;
    } else if (std::strcmp(argv[2], "o") != 0) {
      printUsage();
    }
  }

  // Containers for isomorphism signatures of partial triangulations
  // already seen and for the result
  std::set<std::string> already_seen, result;

  // Initialize a triangulation with one unglued tetrahedron
  Triangulation<3> N;
  N.newTetrahedron();

  // Start the recursion
  recurse(&N, number_tets, orientable, &already_seen, &result);

  // And output the result, std::set's iterator will go
  // through this in lexicographic order.
  for(std::set<std::string>::const_iterator it = result.begin();
      it != result.end();
      ++it) {
    std::cout << *it << std::endl;
  }
}
