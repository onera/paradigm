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

parser.add_argument("-f", "--filename", default="/scratchm/khoogvel/workspace/build/paradigm/impi_gcc12/test/pdm_t_block_to_part_p_d/btp_create.txt") # circulente

parser.add_argument("-s", "--size", default=10)
parser.add_argument("-ssf", "--size_shift", default=3) # quasi-diagonal
parser.add_argument("-csf", "--center_shift", default=0) # circulente

parser.add_argument("-r", "--randomize", action="store_true") # random partition

args = parser.parse_args()

# Size
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

# Create BTP
btp = PDM.BlockToPart(block_distribution,
                      comm,
                      [part_ln_to_gn],
                      n_part)

# Exchange
part_stride, part_data = btp.exchange_field(block_data) # utiliser le inplace ?

# Post-process "create"
if i_rank == 0:
  filename = args.filename

  dico = {}
  with open(filename, "r") as f:
    for line in f:
      values = line.split()

      if values[0] == "n_rank":
        n_rank = int(values[1])

      if values[0] == "i_rank":
        i_rank = int(values[1])
        dico[i_rank] = []

      if values[0] == "node":
        node = int(values[1])
        dico[i_rank].append(node)

      if values[0] == "elaps":
        elaps = float(values[1])
        dico[i_rank].append(elaps)

      if values[0] == "cpu":
        cpu = float(values[1])
        dico[i_rank].append(cpu)

      if values[0] == "n_send":
        n_send = np.zeros(n_rank)
        for j_rank in range(n_rank):
          n_send[j_rank] = int(values[j_rank + 1])
        dico[i_rank].append(n_send)

  # Plot communication graph
  mat = np.zeros((n_rank, n_rank), dtype=int)

  for j_rank in range(n_rank):
    mat[:,j_rank] = dico[j_rank][3][:]

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
      if dico[j_rank][3][k_rank] > 0:
        distance = 0 # no exchange
        if (k_rank != j_rank):
          if abs(dico[j_rank][0]-dico[k_rank][0]) == 0:
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

