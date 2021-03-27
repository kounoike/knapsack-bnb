import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def draw_problem(p, w):
  fig = plt.figure(figsize=(8, 8))
  ax = fig.add_subplot(1, 1, 1, xlabel="weight", ylabel="profit")
  ax.scatter(w, p, s=1)
  ax.set_xlim(0, ax.get_xlim()[1])
  ax.set_ylim(0, ax.get_ylim()[1])
  fig.savefig("problem.png")

with open(sys.argv[1]) as f:
  num = int(f.readline())
  p = np.empty(num)
  w = np.empty(num)
  for idx in range(num):
    id_, p_, w_ = map(int, f.readline().split())
    p[idx] = p_
    w[idx] = w_

draw_problem(p, w)