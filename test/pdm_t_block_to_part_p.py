#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import argparse
import matplotlib.pyplot as plt
import Pypdm.Pypdm as PDM
from Pypdm.Pypdm import npy_pdm_gnum_dtype

# MPI
comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size

# Command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--filename", default="/scratchm/khoogvel/workspace/build/paradigm/impi_gcc12/test/pdm_t_block_to_part_p_d/output.txt") # circulente

parser.add_argument("-n", "--n_iter", default=1)

parser.add_argument("-s", "--size", default=10)
parser.add_argument("-ssf", "--size_shift", default=3) # quasi-diagonal
parser.add_argument("-csf", "--center_shift", default=0) # circulente

parser.add_argument("-r", "--randomize", action="store_true") # random partition

args = parser.parse_args()

# Size (TO DO: size +/- 20% rand ?)
size = int(args.size)

# Block
block_distribution = np.array([i*size for i in range(n_rank+1)]).astype(PDM.npy_pdm_gnum_dtype)
block_data = np.array([1 for i in range(size)]).astype(np.intc)

# Part
n_part = 1

size_shift   = int(args.size_shift)
center_shift = int(args.center_shift)

randomize = args.randomize

if randomize:
  part_ln_to_gn = np.random.randint(1, size * n_rank + 1, size, dtype=PDM.npy_pdm_gnum_dtype)
else:

  beg = ((i_rank + center_shift)%n_rank)*size - size_shift + 1
  if i_rank == 0:
    beg = ((i_rank + center_shift)%n_rank)*size + 1

  end = ((i_rank + center_shift)%n_rank)*size + size_shift + size
  if i_rank == (n_rank-1):
    end = ((i_rank + center_shift)%n_rank)*size + size
  part_ln_to_gn = np.random.randint(beg, end+1, size, dtype=PDM.npy_pdm_gnum_dtype)

# Iterations
n_iter = args.n_iter
for i in range(n_iter):

  # MPI Barrier
  comm.Barrier()

  # Create BTP
  btp = PDM.BlockToPart(block_distribution,
                        comm,
                        [part_ln_to_gn],
                        n_part)

  # Output communication graph
  if i == 0:
    btp.comm_graph_dump("output.txt")

  # Exchange
  part_stride, part_data = btp.exchange_field(block_data)

  # MPI Barrier
  comm.Barrier()

# Output timings
PDM.time_per_step_dump(comm,
                       "output.txt")

# Post-process "create"
if i_rank == 0:
  filename = args.filename

  times = {}
  comm_graph = {}

  times["binary_search"]   = {}
  times["binary_search"]["min_elaps"]  = []
  times["binary_search"]["mean_elaps"] = []
  times["binary_search"]["max_elaps"]  = []
  times["binary_search"]["min_cpu"]  = []
  times["binary_search"]["mean_cpu"] = []
  times["binary_search"]["max_cpu"]  = []

  times["create_exchange"] = {}
  times["create_exchange"]["min_elaps"]  = []
  times["create_exchange"]["mean_elaps"] = []
  times["create_exchange"]["max_elaps"]  = []
  times["create_exchange"]["min_cpu"]  = []
  times["create_exchange"]["mean_cpu"] = []
  times["create_exchange"]["max_cpu"]  = []

  times["data_exchange"]   = {}
  times["data_exchange"]["min_elaps"]  = []
  times["data_exchange"]["mean_elaps"] = []
  times["data_exchange"]["max_elaps"]  = []
  times["data_exchange"]["min_cpu"]  = []
  times["data_exchange"]["mean_cpu"] = []
  times["data_exchange"]["max_cpu"]  = []

  with open(filename, "r") as f:
    for line in f:
      values = line.split()

      if values[0] == "i_rank":
        i_rank = int(values[1])
        comm_graph[i_rank] = {}

      if values[0] == "node":
        node = int(values[1])
        comm_graph[i_rank]["node"] = node

      if values[0] == "n_send":
        n_send = np.zeros(n_rank)
        for j_rank in range(n_rank):
          n_send[j_rank] = int(values[j_rank + 1])
        comm_graph[i_rank]["n_send"] = n_send

      if values[0] == "binary_search":
        times["binary_search"]["min_elaps"].append(values[2])
        times["binary_search"]["mean_elaps"].append(values[3])
        times["binary_search"]["max_elaps"].append(values[4])
        times["binary_search"]["min_cpu"].append(values[6])
        times["binary_search"]["mean_cpu"].append(values[7])
        times["binary_search"]["max_cpu"].append(values[8])

      if values[0] == "create_exchange":
        times["create_exchange"]["min_elaps"].append(values[2])
        times["create_exchange"]["mean_elaps"].append(values[3])
        times["create_exchange"]["max_elaps"].append(values[4])
        times["create_exchange"]["min_cpu"].append(values[6])
        times["create_exchange"]["mean_cpu"].append(values[7])
        times["create_exchange"]["max_cpu"].append(values[8])


      if values[0] == "data_exchange":
        times["data_exchange"]["min_elaps"].append(values[2])
        times["data_exchange"]["mean_elaps"].append(values[3])
        times["data_exchange"]["max_elaps"].append(values[4])
        times["data_exchange"]["min_cpu"].append(values[6])
        times["data_exchange"]["mean_cpu"].append(values[7])
        times["data_exchange"]["max_cpu"].append(values[8])

  # Print
  print(times)

  # Plot communication graph
  mat = np.zeros((n_rank, n_rank), dtype=int)

  for j_rank in range(n_rank):
    mat[:,j_rank] = comm_graph[j_rank]["n_send"][:]

  fig, ax = plt.subplots(figsize=(4,4), dpi=300)

  ax.imshow(mat, cmap=plt.colormaps["Blues"])

  ax.set_aspect("equal")
  ax.set_xlabel("rank (block)")
  ax.set_ylabel("rank (part)")

  # fichier image
  output_name = "/stck/khoogvel/workspace/comm_graph.svg"

  plt.savefig(output_name, bbox_inches="tight")

  plt.close()

  # Plot distance graph (TO DO: dans le bon sens?)
  mat = np.zeros((n_rank, n_rank), dtype=int)

  for j_rank in range(n_rank):
    for k_rank in range(n_rank):
      if comm_graph[j_rank]["n_send"][k_rank] > 0:
        distance = 0 # no exchange
        if (k_rank != j_rank):
          if comm_graph[j_rank]["node"]-comm_graph[k_rank]["node"] == 0:
            distance = 2 # same node
          else:
            distance = 3 # different node
        else:
          distance = 1 # same rank
        mat[k_rank,j_rank] = distance

  fig, ax = plt.subplots(figsize=(4,4), dpi=300)

  ax.imshow(mat, cmap=plt.colormaps["Blues"])

  ax.set_aspect("equal")
  ax.set_xlabel("rank (block)")
  ax.set_ylabel("rank (part)")

  # fichier image
  output_name = "/stck/khoogvel/workspace/distance_graph.svg"

  plt.savefig(output_name, bbox_inches="tight")

  plt.close()

