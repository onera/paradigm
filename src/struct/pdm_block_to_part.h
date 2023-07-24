/*
 * \file
 * \author equemera
 *
 * \date April 14, 2016, 7:56 AM
 */

#ifndef PDM_BLOCK_TO_PART_H
#define	PDM_BLOCK_TO_PART_H

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct PDM_block_to_part_t
 * \brief  Block to partition redistribution
 *
 */

typedef struct _pdm_block_to_part_t PDM_block_to_part_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Reset global statistic
 *
 */

void
PDM_block_to_part_global_statistic_reset
(
 void
);


/**
 *
 * \brief Get global timer in part to block
 *
 * \param [in]   comm                 MPI communicator
 * \param [out]  min_exch_rank_send   Global min part of ranks used to send 
 * \param [out]  min_exch_rank_recv   Global min part of ranks used to receive 
 * \param [out]  max_exch_rank_send   Global max part of ranks used to send 
 * \param [out]  max_exch_rank_recv   Global max part of ranks used to receive
 * \param [out]  min_exch_data_send   Global min sent data for a rank 
 * \param [out]  min_exch_data_recv   Global min received data for a rank
 * \param [out]  max_exch_data_send   Global max sent data for a rank
 * \param [out]  max_exch_data_recv   Global max received data for a rank
 * 
 */

void
PDM_block_to_part_global_statistic_get
(
 PDM_MPI_Comm comm,
 int *min_exch_rank_send,
 int *min_exch_rank_recv,
 int *max_exch_rank_send,
 int *max_exch_rank_recv,
 unsigned long long *min_exch_data_send,
 unsigned long long *min_exch_data_recv,
 unsigned long long *max_exch_data_send,
 unsigned long long *max_exch_data_recv
);

/**
 *
 * \brief Get global timer in block to part
 *
 * \param [in]   comm              MPI communicator
 * \param [out]  min_elaps         Min elapsed time
 * \param [out]  max_elaps         Max elapsed time
 * \param [out]  min_cpu           Min cpu time
 * \param [out]  max_cpu           Max cpu time
 * \param [out]  min_elaps_create  Global min elapsed for create function
 * \param [out]  max_elaps_create  Global max elapsed for create function
 * \param [out]  min_cpu_create    Global min cpu for create function
 * \param [out]  max_cpu_create    Global max cpu for create function
 * \param [out]  min_elaps_exch    Global min elapsed for exch function
 * \param [out]  max_elaps_exch    Global max elapsed for exch function
 * \param [out]  min_cpu_exch      Global min cpu for exch function
 * \param [out]  max_cpu_exch      Global max cpu for exch function
 * 
 */

void
PDM_block_to_part_global_timer_get
(
 PDM_MPI_Comm comm,
 double       *min_elaps_create,
 double       *max_elaps_create,
 double       *min_cpu_create,
 double       *max_cpu_create,
 double       *min_elaps_exch,
 double       *max_elaps_exch,
 double       *min_cpu_exch,
 double       *max_cpu_exch
);

/**
 *
 * \brief Create a block to partitions redistribution
 *
 * \param [in]   block_distrib_idx Block distribution (size : \ref size of \ref comm + 1)
 * \param [in]   gnum_elt          Element global number (size : \ref n_part)
 * \param [in]   n_elt             Local number of elements (size : \ref n_part)
 * \param [in]   n_part            Number of partition
 * \param [in]   comm              MPI communicator
 *
 * \return   Initialized \ref PDM_block_to_part instance
 *
 */

PDM_block_to_part_t *
PDM_block_to_part_create
(
 const PDM_g_num_t   *block_distrib_idx,
 const PDM_g_num_t  **gnum_elt,
 const int           *n_elt,
 const int            n_part,
 const PDM_MPI_Comm   comm
);

/**
 *
 * \brief Create a block to partitions with sparse representation
 *
 * \param [in]   delt_gnum             Sorted list of gnum that represent the partial block, should be betwenn [1, N]
 * \param [in]   dn_elt                Number of element in the partial block
 * \param [in]   gnum_elt              Element global number (size : \ref n_part)
 * \param [in]   n_elt                 Local number of elements (size : \ref n_part)
 * \param [in]   n_part                Number of partition
 * \param [in]   comm                  MPI communicator
 *
 * \return   Initialized \ref PDM_block_to_part instance
 *
 */
PDM_block_to_part_t *
PDM_block_to_part_create_from_sparse_block
(
 const PDM_g_num_t     *delt_gnum,
 const int              dn_elt,
 const PDM_g_num_t    **gnum_elt,
 const int             *n_elt,
 const int              n_part,
 const PDM_MPI_Comm     comm
);


PDM_block_to_part_t *
PDM_block_to_part_create_cf
(
 const PDM_g_num_t    *block_distrib_idx,
 const PDM_g_num_t    **gnum_elt,
 const int            *n_elt,
 const int             n_part,
 const PDM_MPI_Fint    fcomm
);


/**
 *
 * \brief Initialize an exchange
 *
 * \param [in]   btp          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
 *                            Constant stride for \ref PDM_STRIDE_VAR
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Partition stride
 * \param [out]  part_data    Partition data
 *
 */

void
PDM_block_to_part_exch_in_place
(
 PDM_block_to_part_t *btp,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                 *block_stride,
 void                *block_data,
 int                **part_stride,
 void               **part_data
);


/**
 *
 * \brief Initialize an exchange
 * (part_stride and part_data are allocated in function)
 *
 * \param [in]   btp          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
 *                            Constant stride for \ref PDM_STRIDE_VAR
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Partition stride
 * \param [out]  part_data    Partition data
 *
 */

void
PDM_block_to_part_exch
(
 PDM_block_to_part_t *btp,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                 *block_stride,
 void                *block_data,
 int               ***part_stride,
 void              ***part_data
);


/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] btp  Block to part structure
 *
 * \return       NULL
 */

PDM_block_to_part_t *
PDM_block_to_part_free
(
 PDM_block_to_part_t *btp
);


/**
 *
 * \brief Return index in the block for a gnum
 *
 * \param [in] ptb         Part to block structure
 * \param [in] gNum        Global number
 *
 * \return  Index
 */

PDM_l_num_t
PDM_block_to_part_gnum_idx_get
(
 PDM_block_to_part_t *btp,
 PDM_g_num_t gNum
);


/**
 *
 * \brief Get the number of partitions
 *
 * \param [in] btp         Block to part structure
 *
 * \return  Number of partitions
 */

int
PDM_block_to_part_n_part_get
(
 PDM_block_to_part_t *btp
 );


/**
 *
 * \brief Get the number of elements in a given partition
 *
 * \param [in] btp         Block to part structure
 * \param [in] i_part      Id of current partition
 *
 * \return  Number of element in the current partition
 */

int
PDM_block_to_part_n_elt_get
(
 PDM_block_to_part_t *btp,
 const int            i_part
 );

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_BLOCK_TO_PART_H */
