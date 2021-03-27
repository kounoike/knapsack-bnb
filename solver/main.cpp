#include "problem.hpp"

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << argv[0]  << " <problem file> <output csvfile>" << std::endl;
    return -1;
  }
  try {
    KnapsackProblem p(argv[1]);
    p.solve();
    p.printBase();
    p.writeCsv(argv[2]);
    // p.print();
    // p.printSorted();
  } catch(std::exception &ex) {
    std::cerr << ex.what() << std::endl;
    return -1;
  }

  return 0;
}