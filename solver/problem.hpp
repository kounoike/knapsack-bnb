#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <queue>
#include <chrono>

typedef unsigned int T;

class KnapsackProblem;

class Subproblem {
public:
  int num;
  std::vector<bool> isFixedToZero;
  std::vector<bool> isFixedToOne;
  T partialProfit;
  T partialWeight;
  int criticalIndex;
  int depth;
  T upperBound;
  const KnapsackProblem* problem;
  bool isSolved = false;
  bool isInfeasible = false;
  std::vector<int> solution;
  T lowerBound;

  Subproblem(const KnapsackProblem& p);
  Subproblem(const Subproblem& parent, int fixValue);
  void toSolution();
  void improveUpperBound();
  void updatePinning();
};

class KnapsackProblem {
public:
  T capacity;
  std::vector<int> ids;
  std::vector<T> unsortedWeights;
  std::vector<T> unsortedProfits;
  std::vector<int> sortMap;
  std::vector<T> weights;
  std::vector<T> profits;
  std::vector<double> ratios;
  std::vector<int> solution = {0};
  T solutionValue;
  bool solved = false;

  int criticalItemIndex;
  T partialProfit;
  T partialWeight;
  std::vector<int> heuristicSolution = {0};
  T heuristicValue;
  std::vector<int> bestSolution;
  T bestValue;
  T U1;
  T U2;
  T U3;
  T U;
  std::vector<int> J0;
  std::vector<int> J1;
  std::vector<int> F;
  int num;

  // solving status
  int numBranchAndBound = 0;
  std::vector<int> numBranchesVariable;

  KnapsackProblem(int n) :
    num(n),
    ids(n),
    unsortedProfits(n),
    unsortedWeights(n),
    sortMap(n),
    weights(n),
    profits(n),
    ratios(n),
    solution(n),
    heuristicSolution(n),
    solved(false),
    numBranchesVariable(n)
  {}
  KnapsackProblem(std::string filename);

  void calcU1();
  void reduction();
  
  void preprocess() {
    calcU1();
    reduction();
  }

  void solve() {
    // solve this problem
    // preprocess();
    // branchAndBound();

    // solve Fundamental Problem
    solveFundamental();
  }

  void solveFundamental();
  void branchAndBound();
  void printBase();
  void print();
  void printSorted();
  void writeCsv(std::string filename);
};