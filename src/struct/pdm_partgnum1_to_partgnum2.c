/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_partgnum1_to_to_partgnum2.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm.h"
#include "pdm_timer.h"
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/


/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Static function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
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
 * \param [in]   send_to_gnum2_idx  Index of data to send to gnum2 from gnum1 
 *                                  (for each part size : \ref n_elt1+1) 
 * \param [in]   send_to_gnum2      Data to send to gnum2 from gnum1 
 * \param [in]   comm               MPI communicator
 *
 * \return   Initialized \ref PDM_partgnum1_to_partgnum2 instance
 *
 */

PDM_partgnum1_to_partgnum2_t *
PDM_partgnum1_to_partgnum2_create
(
 const PDM_g_num_t   **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t   **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **gnum1_to_gnum2_idx,
 const PDM_g_num_t   **gnum1_to_gnum2,
 const PDM_MPI_Comm    comm
)
{
  PDM_part_to_block_t *ptp =
    (PDM_part_to_block_t *) malloc (sizeof(PDM_part_to_block_t));

  ptp->n_part1                 = n_part1;    
  ptp->gnum_elt1               = gnum_elt1;
  ptp->n_elt1                  = n_elt1;
  ptp->n_part2                 = n_part2;                
  ptp->gnum_elt2               = gnum_elt2;              
  ptp->n_elt2                  = n_elt2;                
  ptp->gnum1_to_gnum2_idx      = gnum1_to_gnum2_idx;    
  ptp->gnum1_to_gnum2          = gnum1_to_gnum2;        
  ptp->comm                    = comm;

  ptp->gnum1_to_send_buffer    = NULL;  
  ptp->recv_buffer_to_gnum2    = NULL;  
  
  ptp->async_n_exch            = 0;           
  ptp->async_l_array           = 0;          
  ptp->async_s_data            = NULL;          
  ptp->async_cst_stride        = NULL;      
  ptp->async_send_request      = NULL;    
  ptp->async_recv_request      = NULL;    
  ptp->async_send_buffer       = NULL;     
  ptp->async_recv_buffer       = NULL;     
  ptp->async_n_send_buffer     = NULL;   
  ptp->async_i_send_buffer     = NULL;   
  ptp->async_n_recv_buffer     = NULL;   
  ptp->async_i_recv_buffer     = NULL;   
}


PDM_partgnum1_to_partgnum2_t *
PDM_partgnum1_to_partgnum2_create_cf
(
 const PDM_g_num_t    **gnum_elt1,
 const int            *n_elt1,
 const int             n_part1,
 const PDM_g_num_t    **gnum_elt2,
 const int            *n_elt2,
 const int             n_part2,
 const int           **gnum1_to_gnum2_idx,
 const PDM_g_num_t   **gnum1_to_gnum2,
 const PDM_MPI_Fint    fcomm
)
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(fcomm);
  return PDM_partgnum1_to_partgnum2_create (gnum_elt1, n_elt1, n_part1,
                                            gnum_elt2, n_elt2, n_part2,
                                            gnum1_to_gnum2_idx, gnum1_to_gnum2);
}


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
PDM_partgnum1_to_partgnum2_exch
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 const size_t                  s_data,
 const PDM_stride_t            t_stride,
 int                         **part1_stride,
 void                        **part1_data
 int                         **part2_stride,
 void                        **part2_data
)
{
  PDM_UNUSED (ptb);
  PDM_UNUSED (s_data);
  PDM_UNUSED (t_stride);
  PDM_UNUSED (part1_stride);
  PDM_UNUSED (part1_data);
  PDM_UNUSED (part2_stride);
  PDM_UNUSED (part2_data);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_exch not yet implemente\n");
}



/**
 *
 * \brief Initialize an exchange with allocation of result arrays
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
PDM_partgnum1_to_partgnum2_exch_with_alloc
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 const size_t               s_data,
 const PDM_stride_t         t_stride,
 int                      **part1_stride,
 void                     **part1_data
 int                     ***part2_stride,
 void                    ***part2_data
)
{
  PDM_UNUSED (ptb);
  PDM_UNUSED (s_data);
  PDM_UNUSED (t_stride);
  PDM_UNUSED (part1_stride);
  PDM_UNUSED (part1_data);
  PDM_UNUSED (part2_stride);
  PDM_UNUSED (part2_data);

  PDM_error(__FILE__, __LINE__, 0,
            "Error PDM_partgnum1_to_partgnum2_exch not yet implemente\n");
}


/**
 *
 * \brief Initialize a asynchronus issend
 *
 * \param [in]   ptp           Block to part structure
 * \param [in]   s_data        Data size
 * \param [in]   cst_stride    Constant stride
 * \param [in]   part1_data    Partition 1 data
 * \param [out]  request       Request
 *
 */

void
PDM_partgnum1_to_partgnum2_issend
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 const size_t                  s_data,
 const int                     cst_stride,
 void                        **part1_data,
 int                          *request
)
{
  
}


/**
 *
 * \brief Wait a asynchronus issend
 *
 * \param [in]  ptp           part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_partgnum1_to_partgnum2_issend_wait
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 int                        request
)
{
  
}


/**
 *
 * \brief Initialize a asynchronus irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  s_data        Data size
 * \param [in]  cst_stride    Constant stride
 * \param [in]  part1_data    Partition 2 data
 * \param [out] request       Request
 *
 */

void
PDM_partgnum1_to_partgnum2_irecv
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 const size_t               s_data,
 const int                  cst_stride,
 void                     **part2_data,
 int                       *request
)
{
  
}


/**
 *
 * \brief Initialize a asynchronus irecv
 *
 * \param [in]  ptp           Part to part structure
 * \param [in]  request       Request
 *
 */

void
PDM_partgnum1_to_partgnum2_irecv_wait
(
 PDM_partgnum1_to_partgnum2_t *ptp,
 int                        request
)
{
  
}


/**
 *
 * \brief Free a block to part structure
 *
 * \param [inout] ptp  Block to part structure
 *
 * \return       NULL
 */

PDM_partgnum1_to_partgnum2_t *
PDM_partgnum1_to_partgnum2_free
(
 PDM_partgnum1_to_partgnum2_t *ptp
)
{
  if (ptp == NULL) {
    return NULL;
  }

  if (ptp->gnum1_to_send_buffer != NULL) {
    free (ptp->gnum1_to_send_buffer);
  }
  if (ptp->recv_buffer_to_gnum2 != NULL) {  
    free (ptp->recv_buffer_to_gnum2);
  }  
  
  if (ptp->async_l_array != 0) {
    for (int i = 0; i < ptp->async_l_array; i++) {    
      free (ptp->async_s_data[i]); 
      free (ptp->async_cst_stride[i]);        
      free (ptp->async_send_request[i]);    
      free (ptp->async_recv_request[i]);    
      free (ptp->async_send_buffer[i]);      
      free (ptp->async_recv_buffer[i]);      
      free (ptp->async_n_send_buffer[i]);  
      free (ptp->async_i_send_buffer[i]);  
      free (ptp->async_n_recv_buffer[i]);  
      free (ptp->async_i_recv_buffer[i]);  
    }
    free (ptp->async_s_data); 
    free (ptp->async_cst_stride);        
    free (ptp->async_send_request);    
    free (ptp->async_recv_request);    
    free (ptp->async_send_buffer);      
    free (ptp->async_recv_buffer);      
    free (ptp->async_n_send_buffer);  
    free (ptp->async_i_send_buffer);  
    free (ptp->async_n_recv_buffer);  
    free (ptp->async_i_recv_buffer);  
  }

  ptp->async_n_exch            = 0;           
  ptp->async_l_array           = 0;          

  free(ptp);
  return NULL;
  
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
