#include "O3Triangulation.h"
#include "O3Tetrahedron.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

void printIsoSig(const std::vector<int> &vec);
std::vector<int> readIsoSig(std::string sig);
std::string paddedNum(int n);

int main(int argc, char **argv)
{  
  // Read in a file containing a list of minimal destination sequences.
  if (argc < 1 || argc > 2) {
    std::cout << "Enter the path of a file containing the list of destination sequences.\n";
    return -1;
  }
  
  std::ifstream data(argv[1]);
  std::string sig;

  std::vector<std::vector<O3Triangulation>> trigs;
  const int max_tets = 36;
  for (int i = 0; i < max_tets; i++) {
    std::vector<O3Triangulation> dummy;
    trigs.push_back(dummy);
  }

  while (data) {
    std::getline(data, sig);
    std::vector<int> destSeq = readIsoSig(sig);

    int num_tets = destSeq.size() / 4;
    O3Triangulation T(destSeq);
    //printIsoSig(T.O3isoSig()); std::cout << "\n";
    if (T.size() >= 1) {
      trigs[num_tets-1].push_back(T);
    }
  }

  std::cout << "Done reading\n";

  for (int n = 0; n < max_tets; n++) {
    for (int j = 0; j < trigs[n].size(); j++) {
      O3Triangulation T = trigs[n][j];
      std::cout << T.size() << " tetrahedra: ";
      printIsoSig(T.O3isoSig());
      std::cout << "\n";
      T.computeCuspCrossSections();
      
      std::string filename = "../singloci/" + paddedNum(n+1) + "t" + paddedNum(j);
      T.printSingularLocus(filename);
    }
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

std::string paddedNum(int n)
{
  std::string result;
  if (0 <= n && n < 10) {
    result = "0" + std::to_string(n);
  }
  else {
    result = std::to_string(n);
  }
  return result;
}
