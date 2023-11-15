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

parser.add_argument("-n", "--n_iter", default=1)

parser.add_argument("-s", "--size", default=10)
parser.add_argument("-ssf", "--size_shift", default=10) # quasi-diagonal (percentage)
parser.add_argument("-csf", "--center_shift", default=0) # circulente

parser.add_argument("-r", "--randomize", action="store_true") # random partition

args = parser.parse_args()

# Size (TO DO: size +/- 20% rand ?)
size = int(args.size)

# Part
n_part = 1

p_size_shift = int(args.size_shift)
size_shift   = size * (p_size_shift/100)
size_shift   = int(size_shift)
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

part_data = np.array([1 for i in range(size)]).astype(np.intc)

# Block
block_distribution = np.array([i*size for i in range(n_rank+1)]).astype(PDM.npy_pdm_gnum_dtype)

filename = args.filename

# Iterations
n_iter = int(args.n_iter)
for i in range(n_iter):

  # MPI Barrier
  comm.Barrier()

  # Create PTB
  btp = PDM.PartToBlock(comm,
                        part_ln_to_gn,
                        None,
                        n_part)

  # Exchange
  block_stride, block_data = btp.exchange_field(part_data)

  # MPI Barrier
  comm.Barrier()
