#ifndef PDM_PARA_OCTREE_CUH
#define	PDM_PARA_OCTREE_CUH

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __CUDACC__
#define MANAGED __managed__
#else
#define MANAGED
#endif

/*============================================================================
 * Type definitions
 *============================================================================*/


/*  */

/*============================================================================
 * Global variable
 *============================================================================*/

 //extern MANAGED PDM_Handles_t *_octrees;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create an octree structure
 *
 * \param [in]   n_point_cloud      Number of point cloud
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 * \param [in]   comm               MPI communicator
 *
 * \return     Identifier
 */

int
PDM_para_octree_create_GPU
(
 const int n_point_cloud,
 const int depth_max,
 const int points_in_leaf_max,
 const int build_leaf_neighbours,
 const PDM_MPI_Comm comm
);

/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_free_GPU
(
 const int          id
);


//test
#ifdef __CUDACC__
__global__
#endif
void
print_from_gpu
(
int id
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PARA_OCTREE_CUH */
