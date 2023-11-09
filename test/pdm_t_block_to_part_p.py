#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import random
import sys
import Pypdm.Pypdm as PDM
from Pypdm.Pypdm import npy_pdm_gnum_dtype

# MPI
comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size

# Size
size = 100
size_per_rank = size // n_rank

# Partition distribution
diagonal = 0
circulante = 1
aleatoire = 0
compact_sort = 0

# Block
block_distribution = np.array([i*size_per_rank for i in range(n_rank+1)]).astype(PDM.npy_pdm_gnum_dtype)
block_data = np.array([1 for i in range(size_per_rank)]).astype(np.intc)

# Part
n_part = 1

if diagonal:
  part_ln_to_gn = np.array([i_rank*size_per_rank + i + 1 for i in range(size_per_rank)]).astype(PDM.npy_pdm_gnum_dtype)

if circulante:
  part_ln_to_gn = np.array([((i_rank+1)%n_rank)*size_per_rank + i + 1 for i in range(size_per_rank)]).astype(PDM.npy_pdm_gnum_dtype)

if aleatoire: # int not long int
  part_ln_to_gn = np.array([random.randint(1, n_rank*size_per_rank) for i in range(size_per_rank)]).astype(PDM.npy_pdm_gnum_dtype)

if compact_sort:
  # to do
  # rank 0 recoit les valeurs de tous les ranks
  # fait le trie
  # envoie aux ranks leur gnum
  alpha = 0.4
  part_value = i_rank + alpha * (np.random.rand(size_per_rank)* 2 - 1)
  part_ln_to_gn = np.argsort(part_value)
  part_ln_to_gn = i_rank*size_per_rank + part_ln_to_gn.astype(PDM.npy_pdm_gnum_dtype) + 1
  print(part_ln_to_gn)

# Create BTP
btp = PDM.BlockToPart(block_distribution,
                      comm,
                      [part_ln_to_gn],
                      n_part)

# Exchange
part_stride, part_data = btp.exchange_field(block_data)

# Rank per node
shared_comm = comm.Split_type(MPI.COMM_TYPE_SHARED)

i_shared_rank = shared_comm.rank

data = shared_comm.gather(i_rank,root=0)

if i_shared_rank == 0:
   print("group {} has ranks {}\n".format(i_rank, data))
