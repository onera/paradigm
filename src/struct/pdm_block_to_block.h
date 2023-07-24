/*
 * \file
 * \author equemera
 *
 * \date January 18, 2018, 9:48 AM
 */

#ifndef PDM_BLOCK_TO_BLOCK_H
#define	PDM_BLOCK_TO_BLOCK_H

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
 * \struct PDM_block_to_block_t
 * \brief  Block to partition redistribution
 *
 */

typedef struct _pdm_block_to_block_t PDM_block_to_block_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a block to another block redistribution
 *
 * \param [in]   blockDistribIdx Block distribution (size : \ref size of \ref comm + 1)
 * \param [in]   gnum_elt        Element global number (size : \ref n_part)
 * \param [in]   n_elt           Local number of elements (size : \ref n_part)
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized \ref PDM_block_to_block instance
 *
 */

PDM_block_to_block_t *
PDM_block_to_block_create
(
 PDM_g_num_t   *block_distrib_ini_idx,
 PDM_g_num_t   *block_distrib_end_idx,
 PDM_MPI_Comm   comm
);


/**
 *
 * \brief Initialize an exchange (Variable stride is not yet available)
 *
 * \param [in]   btb          Block to part structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Constant stride
 * \param [in]   block_stride Stride for each block element for \ref PDM_STRIDE_VAR
 *                            Constant stride for \ref PDM_STRIDE_VAR
 * \param [in]   block_data   Block data
 * \param [out]  part_stride  Partition stride
 * \param [out]  part_data    Partition data
 *
 */

int
PDM_block_to_block_exch
(
 PDM_block_to_block_t *btb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                 *block_stride_ini,
 void                *block_data_ini,
 int                 *block_stride_end,
 void               **block_data_end
);

/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] btb  Block to part structure
 *
 * \return       NULL
 */

PDM_block_to_block_t *
PDM_block_to_block_free
(
 PDM_block_to_block_t *btb
);


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_BLOCK_TO_BLOCK_H */

