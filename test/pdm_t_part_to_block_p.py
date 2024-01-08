#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import argparse
import Pypdm.Pypdm as PDM
from Pypdm.Pypdm import npy_pdm_gnum_dtype

# MPI
comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size

# Command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--filename", default="output.txt")

parser.add_argument("-n", "--n_iter", default=10)

parser.add_argument("-s", "--size", default=10)
parser.add_argument("-sp", "--shift_percentage", default=0) # quasi-diagonal (ex : 10% == 0.1)

parser.add_argument("-r", "--randomize", action="store_true") # random partition
parser.add_argument("-w", "--weight", action="store_true") # weights to one
parser.add_argument("-d", "--distrib", action="store_true") # user distrib for true diagonal mode

args = parser.parse_args()

# Size
size = int(args.size)

# Part
n_part = 1

shift_percentage = float(args.shift_percentage)
shift_size = n_rank * shift_percentage
shift_size = int(shift_size)

randomize = args.randomize

if randomize:
  part_ln_to_gn = np.random.randint(1, size * n_rank + 1, size, dtype=PDM.npy_pdm_gnum_dtype)
else:
  beg = ((i_rank)%n_rank)*size - shift_size*size + 1
  beg = max(1, min(beg, n_rank*size))
  end = ((i_rank)%n_rank)*size + shift_size*size + size
  end = max(1, min(end, n_rank*size))
  part_ln_to_gn = np.random.randint(beg, end+1, size, dtype=PDM.npy_pdm_gnum_dtype)

part_data   = np.array([1 for i in range(size)]).astype(np.intc)

weight = args.weight

part_weight = None
if (weight):
  part_weight = [np.array([1. for i in range(size)]).astype(np.double)]

# Block
distrib = args.distrib

block_distribution = None
if (distrib):
  block_distribution = np.array([i*size for i in range(n_rank+1)]).astype(PDM.npy_pdm_gnum_dtype)

filename = args.filename

# Iterations
n_iter = int(args.n_iter)
for i in range(n_iter):

  # MPI Barrier
  comm.Barrier()

  # Create PTB
  ptb = PDM.PartToBlock(comm,
                        [part_ln_to_gn],
                        pWeight=part_weight,
                        partN=n_part,
                        t_distrib=0,
                        t_post=0,
                        userDistribution=block_distribution)

  # Output communication graph
  if i == 0:
    ptb.comm_graph_dump(filename)

  # Exchange
  block_stride, block_data = ptb.exchange_field([part_data])

  # MPI Barrier
  comm.Barrier()

# Output timings
PDM.ptb_time_per_step_dump(comm,
                           filename)

if i_rank == 0:
  print("End :)")
