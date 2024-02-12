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

parser.add_argument("-w", "--write", action="store_true") # Activate performance file write

parser.add_argument("-f", "--filename", default="output.txt")

parser.add_argument("-n", "--n_iter", default=10)

parser.add_argument("-s", "--size", default=10)
parser.add_argument("-sp", "--shift_percentage", default=0) # quasi-diagonal (ex : 10% == 0.1)

parser.add_argument("-r", "--randomize", action="store_true") # random partition

args = parser.parse_args()

# Size
size = int(args.size)

# Block
block_distribution = np.array([i*size for i in range(n_rank+1)]).astype(PDM.npy_pdm_gnum_dtype)
block_data = np.array([1 for i in range(size)]).astype(np.intc)

# Part
n_part = 1

shift_percentage = float(args.shift_percentage)
shift_size = n_rank * shift_percentage
shift_size = int(shift_size)

randomize = args.randomize

if randomize:
  part_ln_to_gn = np.random.randint(1, size * n_rank + 1, size, dtype=PDM.npy_pdm_gnum_dtype)
else:
  beg = i_rank*size - shift_size*size
  end = i_rank*size + shift_size*size + size

  part_ln_to_gn = np.random.randint(beg, end, size, dtype=PDM.npy_pdm_gnum_dtype)

  for i in range(len(part_ln_to_gn)):
    g = part_ln_to_gn[i]
    if (g >= 0) and (g < n_rank*size):
      part_ln_to_gn[i] = g + 1
    elif g < 0:
      part_ln_to_gn[i] = n_rank*size + g + 1
    elif g >= n_rank*size:
      part_ln_to_gn[i] = g - n_rank*size + 1

filename = args.filename

# Iterations
n_iter = int(args.n_iter)
for i in range(n_iter):

  # MPI Barrier
  comm.Barrier()

  # Create BTP
  btp = PDM.BlockToPart(block_distribution,
                        comm,
                        [part_ln_to_gn],
                        n_part)

  # Output communication graph
  if args.write:
    if i == 0:
      btp.comm_graph_dump(filename)

  # Exchange
  part_stride, part_data = btp.exchange_field(block_data)

  # MPI Barrier
  comm.Barrier()

# Output timings
if args.write:
  PDM.btp_time_per_step_dump(comm,
                             filename)

if i_rank == 0:
  print("End :)")
