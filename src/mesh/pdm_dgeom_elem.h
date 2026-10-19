#ifndef __PDM_DGEOM_ELEM_H__
#define __PDM_DGEOM_ELEM_H__

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
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief  Compute the 3D center coordinates of distributed entities (entity1)
 *         from their descending connectivity to other distributed entities (entity2)
 *
 * \param [in]  dentity1_entity2_idx  Connectivity index (size = dn_entity1+1)
 * \param [in]  dentity1_entity2      Connectivity array
 * \param [in]  dn_entity1            Local number of entity1 on the current MPI rank
 * \param [in]  dentity2_distrib      Distribution of entity2
 * \param [out] dentity1_coord        Entity1 coordinate (size = 3 * dn_entity1)
 * \param [in]  dentity2_coord        Entity2 coordinate (size = 3 * dn_entity2)
 * \param [in]  comm                  MPI communicator
 *
 */
void
PDM_compute_center_from_descending_connectivity
(
  const int          *dentity1_entity2_idx,
  const PDM_g_num_t  *dentity1_entity2,
  const int           dn_entity1,
  const PDM_g_num_t  *dentity2_distrib,
        double       *dentity1_coord,
        double       *dentity2_coord,
        PDM_MPI_Comm  comm
);

/**
 * \brief Compute the characteristic length for each distributed vertex
 *        (defined as the minimum distance to its connected neighboring vertices)
 *
 * \param [in]  comm               MPI communicator
 * \param [in]  dn_face            Number of distributed faces
 * \param [in]  dn_edge            Number of distributed edges
 * \param [in]  dn_vtx             Number of distributed vertices
 * \param [in]  dface_vtx_idx      Connectivity index (size = dn_face+1 )
 * \param [in]  dface_vtx          Distributed face-vertex connectivity array containing global vertex IDs (used if dedge_vtx is NULL)
 * \param [in]  dedge_vtx          Distributed edge-vertex connectivity array containing global vertex IDs (2 vertices per edge, can be NULL)
 * \param [in]  dvtx_coord         Input 3D coordinates of vertices (size: 3 * dn_vtx)
 * \param [out] dchar_length_out   Allocated array containing the characteristic length (size: dn_vtx)
 *
 */
void
PDM_compute_vtx_characteristic_length
(
  PDM_MPI_Comm    comm,
  int             dn_face,
  int             dn_edge,
  int             dn_vtx,
  int            *dface_vtx_idx,
  PDM_g_num_t    *dface_vtx,
  PDM_g_num_t    *dedge_vtx,
  double         *dvtx_coord,
  double        **dchar_length_out
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DGEOM_ELEM_H__ */
