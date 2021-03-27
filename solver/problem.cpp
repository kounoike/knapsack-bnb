#include "problem.hpp"

Subproblem::Subproblem(const KnapsackProblem& p): problem(&p), num(p.num), isFixedToZero(p.num), isFixedToOne(p.num), depth(0) {
  partialProfit = 0;
  partialWeight = 0;
  criticalIndex = p.criticalItemIndex;
  for (int i = 0; i < p.criticalItemIndex; ++i) {
      partialProfit += p.profits[i];
      partialWeight += p.weights[i];
  }
  upperBound = p.U;
  lowerBound = p.bestValue;
}

Subproblem::Subproblem(const Subproblem& parent, int fixValue):
  num(parent.num),
  isFixedToOne(parent.isFixedToOne),
  isFixedToZero(parent.isFixedToZero),
  partialProfit(parent.partialProfit),
  partialWeight(parent.partialWeight),
  criticalIndex(parent.criticalIndex),
  depth(parent.depth + 1),
  problem(parent.problem)
{
  if (fixValue == 0) {
    isFixedToZero[criticalIndex] = true;
    for (int j = criticalIndex + 1; j < num; ++j) {
      if (isFixedToOne[j] || isFixedToZero[j]) {
        continue;
      } else if (partialWeight + problem->weights[j] < problem->capacity) {
        partialWeight += problem->weights[j];
        partialProfit += problem->profits[j];
      } else if (partialWeight + problem->weights[j] == problem->capacity) {
        // std::cout << "pittari(0) at " << j << " criticalIndex: " << criticalIndex << std::endl;;
        partialWeight += problem->weights[j];
        partialProfit += problem->profits[j];
        criticalIndex = j + 1;
        isSolved = true;
        toSolution();
        // std::cout << "  solution value: " << partialProfit << std::endl;
        return;
      } else {
        criticalIndex = j;
        upperBound = partialProfit + std::floor((problem->capacity - partialWeight) * problem->ratios[criticalIndex]);
        return;
      }
    }
    // 全部入った場合
    isSolved = true;
    criticalIndex = num;
    upperBound = partialProfit;
    toSolution();
    // std::cout << "can contain all items." << std::endl;
    return;
  } else if (fixValue == 1) {
    isFixedToOne[criticalIndex] = true;
    partialProfit += problem->profits[criticalIndex];
    partialWeight += problem->weights[criticalIndex];
    for (int j = criticalIndex - 1; j >= 0; --j) {
      if (isFixedToOne[j] || isFixedToZero[j]) {
        continue;
      } else if (partialWeight - problem->weights[j] > problem->capacity) {
        partialProfit -= problem->profits[j];
        partialWeight -= problem->weights[j];
      } else if (partialWeight - problem->weights[j] == problem->capacity) {
        partialProfit -= problem->profits[j];
        partialWeight -= problem->weights[j];
        isSolved = true;
        // std::cout << "pittari(1) at " << j << " criticalIndex: " << criticalIndex << std::endl;;
        criticalIndex = j;
        upperBound = partialProfit;
        toSolution();
        // std::cout << "  solution value: " << partialProfit << std::endl;
        return;
      } else {
        partialProfit -= problem->profits[j];
        partialWeight -= problem->weights[j];
        criticalIndex = j;
        upperBound = partialProfit + std::floor((problem->capacity - partialWeight) * problem->ratios[j]);
        return;
      }
    }
    // 全部外しても容量オーバーの場合
    upperBound = 0;
    // std::cout << "infeasible subproblem" << std::endl;
    isInfeasible = true;
  }
}

void Subproblem::toSolution() {
  partialProfit = 0;
  solution.resize(num);
  for (int i = 0; i < num; ++i) {
    if (isFixedToZero[i]) {
      solution[i] = 0;
    } else if (isFixedToOne[i]) {
      solution[i] = 1;
      partialProfit += problem->profits[i];
    } else if (i < criticalIndex) {
      solution[i] = 1;
      partialProfit += problem->profits[i];
    } else {
      solution[i] = 0;
    }
  }
}

void Subproblem::improveUpperBound() {
  // T U2_0 = upperBound;
  // T U2_1 = upperBound;
  // for (int i = criticalIndex + 1; i < num; ++i) {
  //   if (!isFixedToZero[i] && !isFixedToOne[i]) {
  //     U2_0 = partialProfit + std::floor((problem->capacity - partialWeight) * problem->ratios[i]);
  //     break;
  //   }
  // }
  // for (int i = criticalIndex - 1; i >= 0; --i) {
  //   if (!isFixedToZero[i] && !isFixedToOne[i]) {
  //     U2_1 = partialProfit + std::floor(problem->profits[criticalIndex] - (problem->weights[criticalIndex] - (problem->capacity - partialWeight)) * problem->ratios[i]);
  //     break;
  //   }
  // }
  // T U2 = std::max(U2_0, U2_1);
  // if (U2 < upperBound) {
  //   upperBound = U2;
  // }

  T tmpp1 = partialProfit + problem->profits[criticalIndex];
  T tmpw1 = partialWeight + problem->weights[criticalIndex];
  T U3_1 = upperBound;
  for (int i = criticalIndex - 1; i >= 0; --i) {
    if (isFixedToOne[i] || isFixedToZero[i]) {
      continue;
    }
    if (tmpw1 - problem->weights[i] <= problem->capacity) {
      tmpp1 -= problem->profits[i];
      tmpw1 -= problem->weights[i];
      U3_1 = tmpp1 + std::floor((problem->capacity - tmpw1) * problem->ratios[i]);
      break;
    } else {
      tmpp1 -= problem->profits[i];
      tmpw1 -= problem->weights[i];
    }
  }

  T U3_0 = upperBound;
  T tmpp0 = partialProfit;
  T tmpw0 = partialWeight;
  for (int i = criticalIndex + 1; i < num; ++i) {
    if (isFixedToZero[i] || isFixedToOne[i]) {
      continue;
    }
    if (tmpw0 + problem->weights[i] > problem->capacity) {
      U3_0 = tmpp0 + std::floor((problem->capacity - tmpw0) * problem->ratios[i]);
      break;
    } else {
      tmpp0 += problem->profits[i];
      tmpw0 += problem->weights[i];
    }
  }
  T U3 = std::max(U3_0, U3_1);
  if (U3 < upperBound) {
    // std::cout << "improved: " << upperBound << " -> " << U3 << " : " << upperBound - U3 << std::endl;
    upperBound = U3;
  }
}

void Subproblem::updatePinning() {
  if (lowerBound < problem->bestValue) {
    lowerBound = problem->bestValue;
  } else {
    return;
  }
  T lb = lowerBound;
  int numUpdate = 0;
  for (int j = 0; j < criticalIndex; ++j) {
    if (isFixedToZero[j] || isFixedToOne[j]) continue;
    T sigma0j = num;
    T sump = partialProfit - problem->profits[j];
    T sumw = partialWeight - problem->weights[j];
    for (int k = criticalIndex; k < num; ++k) {
      if (sumw + problem->weights[k] > problem->capacity) {
        sigma0j = k;
        break;
      } else {
        sump += problem->profits[k];
        sumw += problem->weights[k];
      }
    }
    lb = std::max(lb, sump);
    T ub;
    if (sigma0j < num) {
      ub = sump + std::floor((problem->capacity - sumw) * problem->ratios[sigma0j]);
    } else {
      ub = sump;
    }
    if (ub < lb) {
      isFixedToOne[j] = true;
      ++numUpdate;
    }
  }
  for (int j = criticalIndex; j < num; ++j) {
    if (isFixedToZero[j] || isFixedToOne[j]) continue;
    int sigma1s = -1;
    T sump = partialProfit + problem->profits[j];
    T sumw = partialWeight + problem->weights[j];
    for (int k = criticalIndex; k >= 0; --k) {
      if (sumw - problem->weights[k] <= problem->capacity) {
        sump -= problem->profits[k];
        sumw -= problem->weights[k];
        sigma1s = k;
        break;
      } else {
        sump -= problem->profits[k];
        sumw -= problem->weights[k];
      }
    }
    lb = std::max(lb, sump);
    T ub;
    if (sigma1s > -1) {
      ub = sump + std::floor((problem->capacity - sumw) * problem->ratios[sigma1s]);
      if (ub < lb) {
        isFixedToZero[j] = true;
        ++numUpdate;
      }
    }
  }
  if (numUpdate > 0) {
    partialProfit = 0;
    partialWeight = 0;
    for (int i = 0; i < num; ++i) {
      if (isFixedToOne[i]) {
        partialProfit += problem->profits[i];
        partialWeight += problem->weights[i];
      }
    }
    for (int i = 0; i < num; ++i) {
      if (isFixedToZero[i] || isFixedToOne[i]) {
        continue;
      } else if (partialWeight + problem->weights[i] > problem->capacity) {
        criticalIndex = i;
        return;
      } else if (partialWeight + problem->weights[i] == problem->capacity){
        partialProfit += problem->profits[i];
        partialWeight += problem->weights[i];
        isSolved = true;
        upperBound = partialProfit;
        toSolution();
        return;
      } else {
        partialProfit += problem->profits[i];
        partialWeight += problem->weights[i];
      }
    }
    isSolved = true;
    upperBound = partialProfit;
    toSolution();
  }
  // std::cout << "update " << numUpdate << std::endl;
}

KnapsackProblem::KnapsackProblem(std::string filename) {
  std::ifstream ifs(filename);
  ifs >> num;
  KnapsackProblem p(num);
  ids.swap(p.ids);
  unsortedProfits.swap(p.unsortedProfits);
  unsortedWeights.swap(p.unsortedWeights);
  sortMap.swap(p.sortMap);
  weights.swap(p.weights);
  profits.swap(p.profits);
  ratios.swap(p.ratios);
  solution.swap(p.solution);
  heuristicSolution.swap(p.heuristicSolution);
  numBranchesVariable.swap(p.numBranchesVariable);

  for (int i = 0; i < num; ++i) {
    ifs >> ids[i] >> unsortedProfits[i] >> unsortedWeights[i];
  }
  ifs >> capacity;

  // sort
  std::vector<int> map(num);
  for (int i = 0; i < num; ++i) {
    map[i] = i;
  }
  std::sort(map.begin(), map.end(), [this](int a, int b){
    double ratioa = (double)this->unsortedProfits[a]/(double)this->unsortedWeights[a];
    double ratiob = (double)this->unsortedProfits[b]/(double)this->unsortedWeights[b];
    if (std::abs(ratioa - ratiob) > 1e-10) {
      return ratioa > ratiob;
    } else {
      return this->unsortedWeights[a] > this->unsortedWeights[b];
    }
  });
  for (int i = 0; i < num; ++i) {
    int idx = map[i];
    sortMap[i] = idx;
    profits[i] = unsortedProfits[idx];
    weights[i] = unsortedWeights[idx];
    ratios[i] = (double)profits[i] / (double)weights[i];
  }
}


void KnapsackProblem::calcU1() {
  T currentWeight = 0;
  T currentProfit = 0;
  for (int j = 0; j < num; ++j) {
    if (currentWeight + weights[j] < capacity) {
      currentWeight += weights[j];
      currentProfit += profits[j];
      heuristicSolution[j] = 1;
    } else if (currentWeight + weights[j] == capacity) {
      // jを入れるとちょうど
      solutionValue = currentWeight + weights[j];
      criticalItemIndex = -1;
      for (int i = 0; i <= j; ++i) {
        solution[i] = 1;
      }
      solved = true;
      return;
    } else {
      // jで入らなくなる
      T remainCapacity = capacity - currentWeight;
      U1 = currentProfit + std::floor(remainCapacity * ratios[j]);
      partialProfit = currentProfit;
      partialWeight = currentWeight;
      
      T U2_0;
      T U2_1;
      if (j < num - 1) {
        U2_0 = currentProfit + std::floor(remainCapacity * ratios[j + 1]);
      } else {
        U2_0 = currentProfit;
      }
      if (j > 0) {
        U2_1 = currentProfit + std::floor(profits[j] - (weights[j] - remainCapacity) * ratios[j - 1]);
      } else {
        U2_1 = currentProfit;
      }
      U2 = std::max(U2_0, U2_1);

      // U3
      int sigma0s = -1;
      int sigma1s = -1;
      T tmpSumW1 = 0;
      for (int k = 0; k < num; ++k) {
        if (tmpSumW1 + weights[k] > capacity - weights[j]) {
          sigma1s = k;
          break;
        } else {
          tmpSumW1 += weights[k];
        }
      }
      T tmpSumW0 = 0;
      for (int k = 0; k < num; ++k) {
        if (k != j) {
          if (tmpSumW0 + weights[k] > capacity) {
            sigma0s = k;
            break;
          } else {
            tmpSumW0 += weights[k];
          }
        }
      }
      T U3_0_p = 0;
      T U3_0_w = 0;
      for (int i = 0; i < sigma0s; ++i) {
        if (i != j) {
          U3_0_p += profits[i];
          U3_0_w += weights[i];
        }
      }
      T U3_0 = U3_0_p + std::floor((capacity - U3_0_w) * ratios[sigma0s]);
      
      T U3_1_p = 0;
      T U3_1_w = 0;
      for (int i = 0; i < sigma1s; ++i) {
        U3_1_p += profits[i];
        U3_1_w += weights[i];
      }
      T U3_1 = profits[j] + U3_1_p + std::floor((capacity - weights[j] - U3_1_w) * ratios[sigma1s]);
      U3 = std::max(U3_0, U3_1);

      U = std::min(std::min(U1, U2), U3);
      criticalItemIndex = j;
      for (int i = j + 1; i < num; ++i) {
        if (currentWeight + weights[i] <= capacity) {
          heuristicSolution[i] = 1;
          currentWeight += weights[i];
          currentProfit += profits[i];
        }
      }
      heuristicValue = currentProfit;
      return;
    }
  }
  // 全部入る
  solved = true;
  solutionValue = 0;
  for (int i = 0; i < num; ++i) {
    solution[i] = 1;
    solutionValue += profits[i];
  }
  U = solutionValue;
}

void KnapsackProblem::reduction() {
  if (solved) return;
  // T lb = 0;
  // for (int j = 0; j < criticalItemIndex; ++j) {
  //   lb += profits[j];
  // }
  T lb = heuristicValue;
  for (int j = 0; j < criticalItemIndex; ++j) {
    T sigma0j = num;
    T sump = partialProfit - profits[j];
    T sumw = partialWeight - weights[j];
    for (int k = criticalItemIndex; k < num; ++k) {
      if (sumw + weights[k] > capacity) {
        sigma0j = k;
        break;
      } else {
        sump += profits[k];
        sumw += weights[k];
      }
    }
    lb = std::max(lb, sump);
    T ub;
    if (sigma0j < num) {
      ub = sump + std::floor((capacity - sumw) * ratios[sigma0j]);
    } else {
      ub = sump;
    }
    if (ub < lb) {
      J1.push_back(j);
    } else {
      F.push_back(j);
    }
  }
  for (int j = criticalItemIndex; j < num; ++j) {
    int sigma1s = -1;
    T sump = partialProfit + profits[j];
    T sumw = partialWeight + weights[j];
    for (int k = criticalItemIndex; k >= 0; --k) {
      if (sumw - weights[k] <= capacity) {
        sump -= profits[k];
        sumw -= weights[k];
        sigma1s = k;
        break;
      } else {
        sump -= profits[k];
        sumw -= weights[k];
      }
    }
    lb = std::max(lb, sump);
    T ub;
    if (sigma1s > -1) {
      ub = sump + std::floor((capacity - sumw) * ratios[sigma1s]);
      if (ub < lb) {
        J0.push_back(j);
      } else {
        F.push_back(j);
      }
    } else {
      F.push_back(j);
    }
  }
}

void KnapsackProblem::solveFundamental() {
  auto startPreprocess = std::chrono::system_clock::now();
  preprocess();
  auto endPreprocess = std::chrono::system_clock::now();
  KnapsackProblem fundamentalProblem(F.size());
  T fixedProfit = 0;
  T fixedWeight = 0;
  for (auto idx: J1) {
    fixedProfit += profits[idx];
    fixedWeight += weights[idx];
  }
  for (int i = 0; i < F.size(); ++i) {
    int idx = F[i];
    fundamentalProblem.profits[i] = profits[idx];
    fundamentalProblem.weights[i] = weights[idx];
    fundamentalProblem.ratios[i] = ratios[idx];
    fundamentalProblem.capacity = capacity - fixedWeight;
  }
  auto endConvert = std::chrono::system_clock::now();
  std::cout << "fixedProfit: " << fixedProfit << std::endl;
  std::cout << "fixedWeight: " << fixedWeight << std::endl;
  fundamentalProblem.calcU1();
  fundamentalProblem.branchAndBound();
  auto endBranchAndBound = std::chrono::system_clock::now();
  numBranchAndBound = fundamentalProblem.numBranchAndBound;
  for (auto idx: J0) {
    solution[idx] = 0;
  }
  for (auto idx: J1) {
    solution[idx] = 1;
  }
  for (int i = 0; i < F.size(); ++i) {
    solution[F[i]] = fundamentalProblem.solution[i];
    numBranchesVariable[F[i]] = fundamentalProblem.numBranchesVariable[i];
  }
  solutionValue = fundamentalProblem.solutionValue + fixedProfit;
  solved = true;
}

void KnapsackProblem::branchAndBound() {
  bestValue = heuristicValue;
  bestSolution = heuristicSolution;
  Subproblem rootSpb(*this);
  // auto compare = [](const Subproblem& a, const Subproblem& b){return a.upperBound < b.upperBound;};
  auto compare = [](const Subproblem& a, const Subproblem& b){if (a.depth != b.depth) return a.depth < b.depth; else return a.upperBound < b.upperBound;};
  std::priority_queue<Subproblem, std::vector<Subproblem>, decltype(compare)> pool{compare};
  pool.push(rootSpb);
  std::cout << "num: " << num << std::endl;
  std::cout << "initial Best Value: " << bestValue << std::endl;
  while(!pool.empty()) {
    Subproblem spb = pool.top();
    pool.pop();
    if (spb.upperBound <= bestValue) {
      // std::cout << "pruned" << std::endl;
      continue;
    }
    // std::cout << spb.upperBound << " " << bestValue << " " << spb.upperBound - bestValue << std::endl;
    ++numBranchAndBound;
    ++numBranchesVariable[spb.criticalIndex];
    spb.updatePinning();
    if (numBranchAndBound % 10000 == 0) {
      std::cout << "iteration: " << numBranchAndBound
        << " pool: " << pool.size()
        << " ub: " << spb.upperBound
        << " best: " << bestValue
        << " gap: " << spb.upperBound - bestValue
        << " depth: " << spb.depth
        << std::endl;
    }
    Subproblem spb0(spb, 0);
    Subproblem spb1(spb, 1);
    if (spb0.isSolved) {
      if (spb0.partialProfit > bestValue) {
        bestSolution = spb0.solution;
        bestValue = spb0.partialProfit;
        std::cout << "spb0 updated " << bestValue << std::endl;
      } else {
        // std::cout << "spb0 not updated " << bestValue << std::endl;
      }
    } else if (!spb0.isInfeasible) {
      spb0.improveUpperBound();
      pool.push(spb0);
      // std::cout << "spb0 pushed" << std::endl;
    } else {
      // std::cout << "spb0 infeasible" << std::endl;
    }
    if (spb1.isSolved) {
      if (spb1.partialProfit > bestValue) {
        bestSolution = spb1.solution;
        bestValue = spb1.partialProfit;
        std::cout << "spb1 updated " << bestValue << std::endl;
      } else {
        // std::cout << "spb1 not updated " << bestValue << std::endl;
      }
    } else if (!spb1.isInfeasible) {
      spb1.improveUpperBound();
      pool.push(spb1);
      // std::cout << "spb1 pushed" << std::endl;
    } else {
      // std::cout << "spb1 infeasible" << std::endl;
    }
  }
  solution = bestSolution;
  solutionValue = bestValue;
  solved = true;
}

void KnapsackProblem::printBase() {
  std::cout << "solved: " << (solved ? "true" : "false") << std::endl;
  std::cout << "num: " << num << std::endl;
  std::cout << "capacity: " << capacity << std::endl;
  std::cout << "U1: " << U1 << std::endl;
  std::cout << "U2: " << U2 << std::endl;
  std::cout << "U3: " << U3 << std::endl;
  std::cout << "U: " << U << std::endl;
  std::cout << "CriticalItemIndex: " << criticalItemIndex << std::endl;
  std::cout << "|J0|: " << J0.size() << std::endl;
  std::cout << "|J1|: " << J1.size() << std::endl;
  std::cout << "|F|: " << F.size() << std::endl;
  std::cout << "HeuristicValue: " << heuristicValue << std::endl;
  std::cout << "Gap: " << U - heuristicValue << std::endl;
  if (solved) {
    std::cout << "SolutionValue: " << solutionValue << std::endl;
    T val = 0;
    for (int i = 0; i < num; ++i) {
      if (solution[i] == 1) {
        val += profits[i];
      }
    }
    std::cout << "sum of solution: " << val << std::endl;
    std::cout << "num B&B: " << numBranchAndBound << std::endl;
  }
}

void KnapsackProblem::print() {
  printBase();
  std::cout << num << " " << capacity << std::endl;
  for (int i = 0; i < num; ++i) {
    std::cout << ids[i] << " " << unsortedProfits[i] << " " << unsortedWeights[i] << std::endl;
  }
}
void KnapsackProblem::printSorted() {
  printBase();
  std::cout << num << " " << capacity << std::endl;
  for (int i = 0; i < num; ++i) {
    std::cout << sortMap[i] << " " << profits[i] << " " << weights[i] << " " << ratios[i] << " " << solution[i] << std::endl;
  }
}

void KnapsackProblem::writeCsv(std::string filename) {
  std::cout << filename << std::endl;
  for (int i = 0; i < 10; ++i) {
    std::cout << sortMap[i] << " " << profits[i] << " " << weights[i] << std::endl;
  }
  std::ofstream ofs(filename);
  ofs << "id,profit,weight,ratio,x,gx,J0,J1,F,numBranches" << std::endl;
  std::sort(J0.begin(), J0.end());
  for (auto idx: J0) {
    ofs << ids[sortMap[idx]] << ","
      << profits[idx] << ","
      << weights[idx] << ","
      << ratios[idx] << ","
      << solution[idx] << ","
      << heuristicSolution[idx] << ","
      << 1 << ","
      << 0 << ","
      << 0 << ","
      << 0 << std::endl;
  }
  std::sort(J1.begin(), J1.end());
  for (auto idx: J1) {
    ofs << ids[sortMap[idx]] << ","
      << profits[idx] << ","
      << weights[idx] << ","
      << ratios[idx] << ","
      << solution[idx] << ","
      << heuristicSolution[idx] << ","
      << 0 << ","
      << 1 << ","
      << 0 << ","
      << 0 << std::endl;
  }
  std::sort(F.begin(), F.end());
  for (auto idx: F) {
    ofs << ids[sortMap[idx]] << ","
      << profits[idx] << ","
      << weights[idx] << ","
      << ratios[idx] << ","
      << solution[idx] << ","
      << heuristicSolution[idx] << ","
      << 0 << ","
      << 0 << ","
      << 1 << ","
      << numBranchesVariable[idx] << std::endl;
  }
}