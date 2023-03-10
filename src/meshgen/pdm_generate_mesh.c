/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_generate_mesh.h"

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a partitionned sphere mesh (2D).
 *
 * \param [in]  elt_type    Mesh element type
 * \param [in]  order       Mesh element order
 * \param [in]  ho_ordering ?
 * \param [in]  radius      Radius of the sphere
 * \param [in]  center_x    x-coordinate of the sphere center
 * \param [in]  center_y    y-coordinate of the sphere center
 * \param [in]  center_z    z-coordinate of the sphere center
 * \param [in]  n_u         Number of points in longitude
 * \param [in]  n_v         Number of points in latitude
 * \param [in]  n_part      Number of mesh partitions
 * \param [in]  part_method Mesh partitionning method
 *
 * \return PDM_part_mesh_t or PDM_part_mesh_nodal_t
 *
 */

// TO DO : add communicator as argument

PDM_generate_mesh_sphere
(
 PDM_Mesh_nodal_elt_t   elt_type,
 int                    order,
 const char           **ho_ordering,
 const double           radius,
 const double           center_x,
 const double           center_y,
 const double           center_z,
 const PDM_g_num_t      n_u,
 const PDM_g_num_t      n_v,
 const int              n_part,
 const PDM_split_dual_t part_method
)
{
  if (elt_type == PDM_MESH_NODAL_TRIA3) {

    assert(order == 1);

    // why would I prefer using PDM_sphere_surf_gen_nodal ?
    // maybe if n_u != n_v put an if on that

    // generate distributed sphere mesh
    assert(n_u == n_v);
    PDM_MPI_Comm comm;
    PDM_dmesh_nodal_t *dmn = malloc(sizeof(PDM_dmesh_nodal_t));
    PDM_sphere_surf_icosphere_gen_nodal(comm, // TO ADD
                                        n_u,
                                        center_x,
                                        center_y,
                                        center_z,
                                        radius,
                                        &dmn);

    // generate partionned sphere mesh
    // TO DO :  do as in PDM_sphere_surf_icosphere_gen_part

    // TO DO : free dmn ?

  } else {
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for element type %d\n", (int) elt_type);
  }

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
