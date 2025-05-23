#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM

comm = MPI.COMM_WORLD

i_rank = comm.rank

if i_rank == 0:
  n_part = 2
  n_entity = [4, 5]
  entity_ln_to_gn = [np.array([1, 5, 6, 2]).astype(PDM.npy_pdm_gnum_dtype),
                     np.array([5, 3, 13, 10, 12]).astype(PDM.npy_pdm_gnum_dtype)]
  ref_part_bound_proc_idx = [np.array([0, 1, 4, 5]).astype(np.int32),
                             np.array([0, 1, 2, 4]).astype(np.int32)]
  ref_part_bound_part_idx = [np.array([0, 0, 1, 1, 4, 4, 5]).astype(np.int32),
                             np.array([0, 1, 1, 1, 2, 2, 4]).astype(np.int32)]
  ref_part_bound          = [np.array([2, 0, 2, 1, 1, 1, 2, 1, 4, 1, 2, 2, 2, 1, 2, 5, 3, 2, 1, 3]).astype(np.int32),
                             np.array([1, 0, 1, 2, 1, 1, 2, 5, 2, 2, 1, 4, 4, 2, 1, 1]).astype(np.int32)]
  ref_part_priority       = [np.array([0, 1, 0, 3]).astype(np.int32),
                             np.array([1, 5, -1, 1, -1]).astype(np.int32)]
elif i_rank == 1:
  n_part = 3
  n_entity = [2, 7, 1]
  entity_ln_to_gn = [np.array([4, 8]).astype(PDM.npy_pdm_gnum_dtype),
                     np.array([1, 2, 9, 4, 5, 8, 11]).astype(PDM.npy_pdm_gnum_dtype),
                     np.array([7]).astype(PDM.npy_pdm_gnum_dtype)]
  ref_part_bound_proc_idx = [np.array([0, 0, 2, 3]).astype(np.int32),
                             np.array([0, 4, 6, 7]).astype(np.int32),
                             np.array([0, 0, 0, 1]).astype(np.int32)]
  ref_part_bound_part_idx = [np.array([0, 0, 0, 0, 2, 2, 3]).astype(np.int32),
                             np.array([0, 3, 4, 6, 6, 6, 7]).astype(np.int32),
                             np.array([0, 0, 0, 0, 0, 0, 1]).astype(np.int32)]
  ref_part_bound          = [np.array([1, 1, 2, 4, 2, 1, 2, 6, 2, 2, 1, 7]).astype(np.int32),
                             np.array([1, 0, 1, 1, 2, 0, 1, 4, 5, 0, 1, 2, 5, 0, 2, 1, 4, 1, 1, 1, 6, 1, 1, 2, 6, 2, 1, 7]).astype(np.int32),
                             np.array([1, 2, 1, 5]).astype(np.int32)]
  ref_part_priority       = [np.array([2, 3]).astype(np.int32),
                             np.array([0, 3, -1, 2, 1, 3, -1]).astype(np.int32),
                             np.array([4]).astype(np.int32)]
elif i_rank == 2:
  n_part = 1
  n_entity = [7]
  entity_ln_to_gn = [np.array([10, 14, 6, 3, 7, 15, 8]).astype(PDM.npy_pdm_gnum_dtype)]
  ref_part_bound_proc_idx = [np.array([0, 3, 6, 6]).astype(np.int32)]
  ref_part_bound_part_idx = [np.array([0, 1, 3, 4, 5, 6, 6]).astype(np.int32)]
  ref_part_bound          = [np.array([3, 0, 1, 3, 4, 0, 2, 2, 1, 0, 2, 4, 7, 1, 1, 2, 7, 1, 2, 6, 5, 1, 3, 1]).astype(np.int32)]
  ref_part_priority       = [np.array([1, -1, 0, 5, 4, -1, 3]).astype(np.int32)]

part_distribution = np.array([0, 2, 5, 6]).astype(PDM.npy_pdm_gnum_dtype)

data = PDM.generate_entity_graph_comm(comm,
                                      part_distribution,
                                      None,
                                      n_part,
                                      n_entity,
                                      entity_ln_to_gn,
                                      None)

for i_part in range(n_part):
  part_bound_proc_idx = data[i_part]['np_part_bound_proc_idx']
  part_bound_part_idx = data[i_part]['np_part_bound_part_idx']
  part_bound          = data[i_part]['np_part_bound']
  part_priority       = data[i_part]['np_part_priority']
  np.array_equal(part_bound_proc_idx, ref_part_bound_proc_idx[i_part])
  np.array_equal(part_bound_part_idx, ref_part_bound_part_idx[i_part])
  np.array_equal(part_bound,          ref_part_bound         [i_part])
  np.array_equal(part_priority,       ref_part_priority      [i_part])
