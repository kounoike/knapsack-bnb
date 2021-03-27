from ortools.algorithms import pywrapknapsack_solver
import sys
import time

with open(sys.argv[1]) as f:
  num = int(f.readline())
  p = []
  w = []
  for idx in range(num):
    id_, p_, w_ = map(int, f.readline().split())
    p.append(p_)
    w.append(w_)
  cap = int(f.readline())

solver = pywrapknapsack_solver.KnapsackSolver(
  pywrapknapsack_solver.KnapsackSolver.KNAPSACK_MULTIDIMENSION_BRANCH_AND_BOUND_SOLVER,
  'KnapsackExample'
)

solver.Init(p, [w], [cap])
start = time.time()
computed_value = solver.Solve()
end = time.time()
packed_items = []
packed_weights = []
total_weight = 0
print('Total value =', computed_value)
# for i in range(len(p)):
#     if solver.BestSolutionContains(i):
#         packed_items.append(i)
#         packed_weights.append(w[i])
#         total_weight += w[i]
print('Total weight:', total_weight)
# print('Packed items:', packed_items)
# print('Packed_weights:', packed_weights)
print('Problem solved in %f seconds' % (end - start))
