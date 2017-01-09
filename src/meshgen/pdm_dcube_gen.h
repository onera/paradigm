#ifndef __PDM_PART_DCUBE_H__
#define __PDM_PART_DCUBE_H__


#include "pdm.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if !defined (__hpux) && !defined (_AIX) 
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/**
 *
 * \brief Create a distributed cube
 *
 * \param [out]  id             dcube identifier
 * \param [in]   comm           Communicator
 * \param [in]   nVtxSeg        Number of vertices in segments
 * \param [in]   length         Segment length
 *
 */

void
PDM_part_dcube_init 
(
 int                *id,
 PDM_MPI_Comm        comm,
 const PDM_g_num_t  nVtxSeg, 
 const double        length
 );

void
PROCF (pdm_part_dcube_init, PDM_PART_DCUBE_INIT)  
(
 int                *id,
 const PDM_MPI_Fint *comm,
 const PDM_g_num_t  *nVtxSeg, 
 const double       *length
);


/**
 *
 * \brief Return distributed cube size
 *
 * \param [in]   id          dcube identifier
 * \param [out]  NFaceGroup  Number of faces groups
 * \param [out]  dNCell      Number of cells stored in this process 
 * \param [out]  dNFace      Number of faces stored in this process
 * \param [out]  dNVtx       Number of vertices stored in this process
 * \param [out]  dFaceVtxL   Length of dFaceVtx array
 * \param [out]  dFacegroupL Length of dFacegroup array
 *
 */

void
PDM_part_dcube_dim_get 
(
 int                id,
 int                *nFaceGroup,
 int                *dNCell,
 int                *dNFace,
 int                *dNVtx,
 int                *dFaceVtxL,
 int                *dFacegrouL
);


void
PROCF(pdm_part_dcube_dim_get, PDM_PART_DCUBE_DIM_GET) 
(
 int                *id,
 int                *nFaceGroup,
 int                *dNCell,
 int                *dNFace,
 int                *dNVtx,
 int                *dFaceVtxL,
 int                *dFacegroupL
);


/**
 *
 * \brief Return distributed cube data
 *
 * \param [in]  id            dcube identifier
 * \param [out] dFaceCell     Faces from cells connectivity (size = 2 * dNFace)
 * \param [out] dFaceVtxIdx   Faces from vertices connectivity index (size = dNface + 1)
 * \param [out] dFaceVtx      Faces from vertices connectivity (size = dFaceVtxL)
 * \param [out] dVtxCoord     Vertices coordinates (size = 3 * dNVtx)
 * \param [out] dFaceGroupIdx Faces groups index (size = NFaceGroup + 1)
 * \param [out] dFaceGroup    Faces groups (size = dFacegroupL)
 *
 */

void
PDM_part_dcube_data_get 
(
 int                 id,
 PDM_g_num_t      **dFaceCell,
 int               **dFaceVtxIdx, 
 PDM_g_num_t      **dFaceVtx,
 double            **dVtxCoord,
 int               **dFaceGroupIdx,
 PDM_g_num_t      **dFaceGroup 
); 


void 
PROCF (pdm_part_dcube_data_get, PDM_PART_DCUBE_DATA_GET)  
(
 int               *id,
 PDM_g_num_t      *dFaceCell,
 int               *dFaceVtxIdx, 
 PDM_g_num_t      *dFaceVtx,
 double            *dVtxCoord,
 int               *dFaceGroupIdx,
 PDM_g_num_t      *dFaceGroup 
); 


/**
 *
 * \brief Free a distributed cube
 *
 * \param [in]  id            dcube identifier
 *
 */

void
PDM_part_dcube_free
(
 int             id
 ); 

void 
PROCF (pdm_part_dcube_free, PDM_PART_DCUBE_FREE)
(
 int            *id
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PART_DCUBE_H__ */
