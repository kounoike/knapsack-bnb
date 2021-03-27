import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, s=1)
    print("scatter done.")

    # now determine nice limits by hand:
    binwidth = 100
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(0, lim + binwidth, binwidth)
    print("drawing histx")
    ax_histx.hist(x, bins=bins)
    print("drawing histy")
    ax_histy.hist(y, bins=bins, orientation='horizontal')

def draw_profit_weight(df):
  x = df.weight.to_numpy()
  y = df.profit.to_numpy()

  # definitions for the axes
  left, width = 0.1, 0.65
  bottom, height = 0.1, 0.65
  spacing = 0.005


  rect_scatter = [left, bottom, width, height]
  rect_histx = [left, bottom + height + spacing, width, 0.2]
  rect_histy = [left + width + spacing, bottom, 0.2, height]

  # start with a square Figure
  fig = plt.figure(figsize=(8, 8))

  ax = fig.add_axes(rect_scatter)
  ax_histx = fig.add_axes(rect_histx, sharex=ax)
  ax_histy = fig.add_axes(rect_histy, sharey=ax)
  ax.set_xlabel("weight")
  ax.set_ylabel("profit")

  # use the previously defined function
  scatter_hist(x, y, ax, ax_histx, ax_histy)

  fig.savefig("profit_weight.png")

def draw_pinned(df):
  x = df.weight.to_numpy()
  y = df.profit.to_numpy()
  color = df.id
  color[df.J0 == 1] = "#ff0000"
  color[df.J1 == 1] = "#0000ff"
  color[df.F == 1] = "#00ff00"
  fig = plt.figure(figsize=(8, 8))
  plt.xlabel("weight")
  plt.ylabel("profit")
  plt.scatter(x, y, c=color, s=1)
  fig.savefig("pinned.png")

def draw_branches(df):
  x = df.weight.to_numpy()
  y = df.profit.to_numpy()
  color = df.id
  color[df.J0 == 1] = "#000000"
  color[df.J1 == 1] = "#000000"
  color[df.F == 1] = "#00ff00"
  size = df.numBranches
  size[df.F == 0] = 1

  fig = plt.figure(figsize=(8, 8))
  plt.xlabel("weight")
  plt.ylabel("profit")
  plt.scatter(x, y, c=color, s=size)
  fig.savefig("branches.png")

def draw_solution(df):
  x = df.weight.to_numpy()
  y = df.profit.to_numpy()
  color = df.id
  color[df.x == 0] = "#000000"
  color[df.x == 1] = "#ff0000"
  color[df.gx == 1] = "#0000ff"
  color[df.gx + df.x == 2] = "#880088"

  fig = plt.figure(figsize=(8, 8))
  plt.xlabel("weight")
  plt.ylabel("profit")
  plt.scatter(x, y, c=color, s=2)
  fig.savefig("solution.png")

df = pd.read_csv(sys.argv[1])
print(df)
draw_profit_weight(df)
draw_pinned(df)
draw_branches(df)
draw_solution(df)