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

  print(i_rank)
  print(part_ln_to_gn)

# Create BTP
btp = PDM.BlockToPart(block_distribution,
                      comm,
                      [part_ln_to_gn],
                      n_part)

# Exchange
part_stride, part_data = btp.exchange_field(block_data)

# Rank per node -> TO DO : mettre dans btp.c
shared_comm = comm.Split_type(MPI.COMM_TYPE_SHARED)

i_shared_rank = shared_comm.rank

data = shared_comm.gather(i_rank,root=0)

if i_shared_rank == 0:
   print("group {} has ranks {}\n".format(i_rank, data))

