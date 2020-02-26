/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_distant_neighbor.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_handles.h"
#include "pdm_mpi.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/**
 * \struct _distant_neighbor_t
 * \brief  Define a point merge structures
 *
 */
typedef struct  {
  PDM_MPI_Comm comm;             /*!< MPI communicator */
  int          n_part;           /*!< Number of partitions */
  const int   *n_entity;         /*!< Number of entities for each partition */
  const int  **neighbor_idx;     /*!< Indexes of candidate for each current part point
                                  *   (size = number of entities in the current part + 1) */
  const int  **neighbor_desc;    /*!< Candidates description (process,
                                  *                           part in the process,
                                  *                           entitiy number in the part) */
} _distant_neighbor_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

static PDM_Handles_t *_pdns   = NULL;

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */
static _distant_neighbor_t *
_get_from_id
(
 int  id
)
{
  _distant_neighbor_t *pdn = (_distant_neighbor_t *) PDM_Handles_get (_pdns, id);

  if (pdn == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_distant_neighbor error : Bad identifier\n");
  }

  return pdn;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Return an initialized \ref _distant_neighbor_t structure
 *
 * This function returns an initialized \ref _distant_neighbor_t structure
 *
 * \param [in]   comm          MPI communicator
 * \param [in]   n_part        Number of partitions
 * \param [in]   n_entity      Number of entities for each partition
 * \param [out]  neighbor_idx  Indexes of candidate for each current part point
 *                              (size = number of entities in the current part + 1)
 * \param [out]  neighbor_desc Candidates description (process,
 *                                                     part in the process,
 *                                                     entitiy in the part)
 *
 * \return      A new initialized \ref PDM_distant_neighbor structure
 *
 */
int
PDM_distant_neighbor_create
(
const PDM_MPI_Comm   comm,
const int            n_part,
const int           *n_entity,
const int          **neighbor_idx,
const int          **neighbor_desc
)
{
  printf(" PDM_distant_neighbor_create \n");

  if (_pdns == NULL) {
    _pdns = PDM_Handles_create (4);
  }

  _distant_neighbor_t *pdn = (_distant_neighbor_t *) malloc(sizeof(_distant_neighbor_t));

  int id = PDM_Handles_store (_pdns, pdn);

  pdn->comm          = comm;
  pdn->n_part        = n_part;
  pdn->n_entity      = n_entity;
  pdn->neighbor_idx  = neighbor_idx;
  pdn->neighbor_desc = neighbor_desc;

  /*
   *  Setup exchange protocol
   */


  return id;
}


/**
 *
 * \brief Free an distant negihtbor structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_distant_neighbor_free
(
 const int          id
)
{
  _distant_neighbor_t *pdn = _get_from_id (id);

  free (pdn);

  PDM_Handles_handle_free (_pdns, id, PDM_FALSE);

  const int n_ppm = PDM_Handles_n_get (_pdns);

  if (n_ppm == 0) {
    _pdns = PDM_Handles_free (_pdns);
  }

}


/**
 * \brief Exchange data between \ref _distant_neighbor_t structure
 * \param [in]   id          identifier of internal structre
 *  NB : On va commencer par des entiers en stride constantes
 *
 */
void
PDM_distant_neighbor_exch
(
 const int      id,
 size_t         s_data,
 PDM_stride_t   t_stride,
 int            cst_stride,
 const int    **send_entity_stride,
 const int    **send_entity_data,
       int    **recv_entity_stride,
       int    **recv_entity_data
)
{
  printf(" PDM_distant_neighbor_exchange \n");
  _distant_neighbor_t *pdn = _get_from_id (id);

  if(t_stride !=  PDM_STRIDE_CST) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_distant_neighbor_exch : STRIDE_CST is only availble \n");
    abort ();
  }


}

#ifdef __cplusplus
}
#endif /* __cplusplus */
