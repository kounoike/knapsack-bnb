#include "problem.hpp"
#include <chrono>

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << argv[0]  << " <problem file> <output csvfile>" << std::endl;
    return -1;
  }
  try {
    auto start = std::chrono::system_clock::now();
    KnapsackProblem p(argv[1]);
    auto endRead = std::chrono::system_clock::now();
    p.solve();
    auto endSolve = std::chrono::system_clock::now();
    p.printBase();
    p.writeCsv(argv[2]);
    auto endOutput = std::chrono::system_clock::now();

    std::cout
      << "Read time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endRead - start).count() << std::endl
      << "Solve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endSolve - endRead).count() << std::endl
      << "Output time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endOutput - endSolve).count() << std::endl;

    // p.print();
    // p.printSorted();
  } catch(std::exception &ex) {
    std::cerr << ex.what() << std::endl;
    return -1;
  }

  return 0;
}