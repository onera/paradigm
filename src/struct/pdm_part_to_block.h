/*
 * \file
 */

#ifndef __PDM_PART_TO_BLOCK_H__
#define __PDM_PART_TO_BLOCK_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_geom.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \enum PDM_part_to_block_distrib_t
 * \brief Type of distribution
 *
 */

typedef enum {

  PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC          = 0,  /*!< Distribute block on all processes */
  PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE = 1,  /*!< Distribute block on one processe pere node */
  PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE      = 2   /*!< Distribute block on part of nodes */

} PDM_part_to_block_distrib_t;


/**
 * \enum PDM_part_to_block_post_t
 * \brief Type of Post processing about blocks
 *
 */

typedef enum {

  PDM_PART_TO_BLOCK_POST_NOTHING       = 0,  /*!< No post processing     */
  PDM_PART_TO_BLOCK_POST_CLEANUP       = 1,  /*!< Cleanup multi-elements */
  PDM_PART_TO_BLOCK_POST_MERGE         = 2,  /*!< Merge multi-elements   */
  PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM = 3   /*!< TMP   */

} PDM_part_to_block_post_t;

/**
 * \enum PDM_part_to_block_distrib_t
 * \brief Type of distribution
 *
 */

/* typedef enum {

 PDM_PART_TO_BLOCK_STRIDE_CST = 0, */  /*!< Constant stride element */
 /* PDM_PART_TO_BLOCK_STRIDE_VAR = 1, */  /*!< Variable stride element */

/* } PDM_part_to_block_stride_t; */


/**
 * \struct PDM_part_to_block_t
 * \brief  Partitioning to block redistribution
 *
 */

typedef struct _pdm_part_to_block_t PDM_part_to_block_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief PDM_extents_conformize
 *
 * Correction extents to manage singular cases and di-symetrizes pb
 * eps = 1.e-3 is a standard value
 *
 */
void
PDM_extents_conformize(int    dim,
                       double extents[],
                       double eps);


/**
 *
 * \brief Reset global statistic
 *
 */

void
PDM_part_to_block_global_statistic_reset
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
PDM_part_to_block_global_statistic_get
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
 * \brief Get global timer in part to block
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
 * \param [out]  min_elaps_create2 Global min elapsed for create2 function
 * \param [out]  max_elaps_create2 Global max elapsed for create2 function
 * \param [out]  min_cpu_create2   Global min cpu for create2 function
 * \param [out]  max_cpu_create2   Global max cpu for create2 function
 * \param [out]  min_elaps_exch    Global min elapsed for exch function
 * \param [out]  max_elaps_exch    Global max elapsed for exch function
 * \param [out]  min_cpu_exch      Global min cpu for exch function
 * \param [out]  max_cpu_exch      Global max cpu for exch function
 *
 */

void
PDM_part_to_block_global_timer_get
(
 PDM_MPI_Comm comm,
 double       *min_elaps_create,
 double       *max_elaps_create,
 double       *min_cpu_create,
 double       *max_cpu_create,
 double       *min_elaps_create2,
 double       *max_elaps_create2,
 double       *min_cpu_create2,
 double       *max_cpu_create2,
 double       *min_elaps_exch,
 double       *max_elaps_exch,
 double       *min_cpu_exch,
 double       *max_cpu_exch
);

/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   partActiveNode  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   weight          Weight of elements (or NULL)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized PDM_part_to_block_t
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                         partActiveNode,
 PDM_g_num_t                 **gnum_elt,
 double                       **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
);


/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   partActiveNode  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   weight          Weight of elements (or NULL)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized PDM_part_to_block_t
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create_from_distrib
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        partActiveNode,
 PDM_g_num_t                 **gnum_elt,
 const PDM_g_num_t            *dataDistribIndex,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
);

PDM_part_to_block_t *
PDM_part_to_block_geom_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_part_geom_t               geom_kind,
 double                      **pvtx_coords,
 PDM_g_num_t                 **gnum_elt,
 int                         **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
);

/**
 *
 * \brief Return number of active ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Number of active ranks
 *
 */

int
PDM_part_to_block_n_active_ranks_get
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return active ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  active ranks
 *
 */

int *
PDM_part_to_block_active_ranks_get
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return if current rank is active
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  if current rank is active
 *
 */

int
PDM_part_to_block_is_active_rank
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return number of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Number of element in the current process
 *
 */

int
PDM_part_to_block_n_elt_block_get
(
 PDM_part_to_block_t *ptb
 );


/**
 *
 * \brief Return global numbers of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global numbers
 *
 */

PDM_g_num_t *
PDM_part_to_block_block_gnum_get
(
 PDM_part_to_block_t *ptb
);

/**
 *
 * \brief Return numbers of occurence of each gnum element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global numbers counter
 *
 */

int *
PDM_part_to_block_block_gnum_count_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 * \param [in]   part_stride  Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   part_data    partitioned data
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 * \return       Size of highest block
 *
 */

int
PDM_part_to_block_exch
(
 PDM_part_to_block_t       *ptb,
 size_t                     s_data,
 PDM_stride_t               t_stride,
 int                        cst_stride,
 int                      **part_stride,
 void                     **part_data,
 int                      **block_stride,
 void                     **block_data
);



void
PDM_part_to_block_reverse_exch
(
 PDM_part_to_block_t *ptb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                 *block_stride,
 void                *block_data,
 int               ***part_stride,
 void              ***part_data
);


/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 * \param [in]   part_stride  Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   part_data    partitioned data
 *
 * \return       Size of highest block
 *
 */

int
PDM_part_to_block_async_exch
(
 PDM_part_to_block_t       *ptb,
 size_t                     s_data,
 PDM_stride_t               t_stride,
 int                        cst_stride,
 int                      **part_stride,
 void                     **part_data
);

/**
 *
 * \brief Initialize a asynchronous data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 * \param [in]   part_stride  Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   part_data    partitioned data
 *
 *
 */
void
PDM_part_to_block_iexch
(
       PDM_part_to_block_t  *ptb,
 const PDM_mpi_comm_kind_t   k_comm,
       size_t                s_data,
       PDM_stride_t          t_stride,
       int                   cst_stride,
       int                 **part_stride,
       void                **part_data,
       int                 **block_stride,
       void                **block_data,
       int                  *request
);

/**
 *
 * \brief Finalize and post-treated a asynchronous data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   id of request to wait / post
 *
 * \return       Size of highest block
 *
 */
int
PDM_part_to_block_iexch_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
);


/**
 *
 * \brief Initialize a asynchronous data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 * \param [in]   block_stride Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   block_data   block data
 * \param [in]   part_stride  Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   part_data    partitioned data
 *
 *
 */
void
PDM_part_to_block_reverse_iexch
(
       PDM_part_to_block_t  *ptb,
 const PDM_mpi_comm_kind_t   k_comm,
       size_t                s_data,
       PDM_stride_t          t_stride,
       int                   cst_stride,
       int                  *block_stride,
       void                 *block_data,
       int                ***part_stride,
       void               ***part_data,
       int                  *request
);

/**
 *
 * \brief Finalize and post-treated a asynchronous data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   id of request to wait / post
 *
 * \return       Size of highest block
 *
 */
void
PDM_part_to_block_reverse_iexch_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
);

/**
 *
 * \brief Wait for an exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 *
 */
void
PDM_part_to_block_async_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
);

/**
 *
 * \brief Get the raw exchange buffer and stride and deallocate memory
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 */
int
PDM_part_to_block_asyn_get_raw
(
 PDM_part_to_block_t *ptb,
 int                  request_id,
 int                **block_stride,
 void               **block_data
);

/**
 *
 * \brief Get the raw exchange buffer and stride and deallocate memory
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 */
int
PDM_part_to_block_asyn_post_treatment
(
 PDM_part_to_block_t *ptb,
 int                  request_id,
 int                **block_stride,
 void               **block_data
);

/**
 *
 * \brief Free a part to block structure
 *
 * \param [inout] ptb         Part to block structure
 *
 * \return       NULL
 */

PDM_part_to_block_t *
PDM_part_to_block_free
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Return block distribution
 *
 * \param [in] ptb         Part to block structure
 *
 * \return  Distribution (size = communicator size + 1)
 */

PDM_g_num_t *
PDM_part_to_block_distrib_index_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Return processus destination
 *
 * \param [in] ptb         Part to block structure
 *
 * \return  Destination (size = sum of partition elements)
 */

PDM_l_num_t *
PDM_part_to_block_destination_get
(
 PDM_part_to_block_t *ptb
);


PDM_g_num_t*
PDM_part_to_block_adapt_partial_block_to_block
(
 PDM_part_to_block_t  *ptb,
 int                 **block_n,
 PDM_g_num_t           n_g_block
);



/**
 *
 * \brief Return global weights of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global weights
 *
 */

double **
PDM_part_to_block_global_weight_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Get number of MPI ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Number of MPI ranks
 *
 */

int
PDM_part_to_block_n_ranks_get
(
 PDM_part_to_block_t *ptb
);


/**
 *
 * \brief Return total number of element in the current process (summed over all partitions)
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Total number of element in the current process
 *
 */

int
PDM_part_to_block_n_elt_proc_get
(
 PDM_part_to_block_t *ptb
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_writer_PART_TO_BLOCK_H__ */
