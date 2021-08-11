#include "O3Triangulation.h"
#include "O3Tetrahedron.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

void printIsoSig(const std::vector<int> &vec);
std::vector<int> readIsoSig(std::string sig);

int main(int argc, char **argv)
{  
  // Read in a file containing a list of minimal destination sequences.
  if (argc < 1 || argc > 2) {
    std::cout << "Enter the path of a file containing the list of destination sequences.\n";
    return -1;
  }
  
  std::ifstream data(argv[1]);
  std::string sig;
  while (data) {
    std::getline(data, sig);
    std::vector<int> destSeq = readIsoSig(sig);
    
    O3Triangulation T(destSeq);
    std::cout << T.size() << " tetrahedra: ";
    printIsoSig(T.O3isoSig());
    std::cout << "\n";
    T.computeCuspCrossSections();
  }
}

void printIsoSig(const std::vector<int> &vec)
{
  int count = 0;
  for (int i = 0; i < vec.size(); i++) {
    std::cout << vec[i] << ",";
    count = count + 1;
    if (count == 4) {
      std::cout << " ";
      count = 0;
    }
  }
}

std::vector<int> readIsoSig(std::string sig)
{
  std::vector<int> seq;
  std::stringstream ss(sig);
  for (int i; ss >> i;) {
    seq.push_back(i);
    if (ss.peek() == ',' || ss.peek() == ' ' || ss.peek() == '\n')
      ss.ignore();
  }
  return seq;
}


