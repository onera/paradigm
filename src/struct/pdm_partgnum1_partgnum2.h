/*
 * File:   pdm_part_to_part.h
 * Author: equemera
 *
 * Created on April 14, 2016, 7:56 AM
 */

#ifndef PDM_PART_TO_PART_H
#define	PDM_PART_tO_pART_H

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
 * \struct PDM_part_to_part_t
 * \brief  Block to partition redistribution
 *
 */

typedef struct _pdm_part_to_part_t PDM_part_to_part_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a partitions to partitions redistribution
 *
 * \param [in]   gnum_elt1          Element global number (size : \ref n_part1)
 * \param [in]   n_elt1             Local number of elements (size : \ref n_part1)
 * \param [in]   n_part1            Number of partition
 * \param [in]   gnum_elt2          Element global number (size : \ref n_part2)
 * \param [in]   n_elt2             Local number of elements (size : \ref n_part2)
 * \param [in]   n_part2            Number of partition
 * \param [in]   comm              MPI communicator
 *
 * \return   Initialized \ref PDM_part_to_part instance
 *
 */

PDM_part_to_part_t *
PDM_part_to_part_create
(
 const PDM_g_num_t    **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t    **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int             n_part2,
 const PDM_MPI_Comm    comm
);


PDM_part_to_part_t *
PDM_part_to_part_create_cf
(
 const PDM_g_num_t    **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t    **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const PDM_MPI_Fint    fcomm
);


/**
 *
 * \brief Initialize an exchange
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   s_data        Data size
 * \param [in]   t_stride      Stride type
 * \param [in]   part1_stride  Partition 1 stride
 * \param [in]   part1_data    Partition 1 data
 * \param [out]  part2_stride  Partition 2 stride
 * \param [out]  part2_data    Partition 2 data
 *
 */

void
PDM_part_to_part_exch
(
 PDM_part_to_part_t *ptp,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                **part1_stride,
 void               **part1_data
 int                **part2_stride,
 void               **part2_data
);


/**
 *
 * \brief Initialize an exchange
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   s_data        Data size
 * \param [in]   t_stride      Stride type
 * \param [in]   part1_stride  Partition 1 stride
 * \param [in]   part1_data    Partition 1 data
 * \param [out]  part2_stride  Partition 2 stride
 * \param [out]  part2_data    Partition 2 data
 *
 */

void
PDM_part_to_part_exch_with_alloc
(
 PDM_part_to_part_t *ptp,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                **part1_stride,
 void               **part1_data
 int                ***part2_stride,
 void               ***part2_data
);


/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] ptp  Block to part structure
 *
 * \return       NULL
 */

PDM_part_to_part_t *
PDM_part_to_part_free
(
 PDM_part_to_part_t *ptp
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_tO_pART_H */
