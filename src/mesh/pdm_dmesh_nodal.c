
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_quick_sort.h"
#include "pdm_geom_elem.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_para_graph_dual.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Maximum number of sections depending of section type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 * \def _compute_keys
 */
static
void
_compute_keys
(
const int          n_face_elt_tot,
const int         *dcell_face_vtx_idx,
const PDM_g_num_t *dcell_face_vtx,
      PDM_g_num_t *ln_to_gn,
      PDM_g_num_t  key_mod
)
{
  for(int i_face = 0; i_face < n_face_elt_tot; ++i_face ) {
    PDM_g_num_t key = 0;
    for(int idx = dcell_face_vtx_idx[i_face]; idx < dcell_face_vtx_idx[i_face+1]; ++idx) {
      key += dcell_face_vtx[idx];
    }
    // min_vtx =
    ln_to_gn[i_face] = key % key_mod;
  }
}


static void
_make_absolute_face_numbering(PDM_dmesh_nodal_t* dmesh_nodal)
{

  PDM_g_num_t n_face_proc = dmesh_nodal->dn_face;
  PDM_g_num_t beg_num_abs;

  PDM_MPI_Scan(&n_face_proc, &beg_num_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dmesh_nodal->pdm_mpi_comm);
  beg_num_abs -= n_face_proc;

  /** Compute the distribution of elements amont proc **/
  dmesh_nodal->face_distrib = (PDM_g_num_t *) malloc((dmesh_nodal->n_rank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dn_face = (PDM_g_num_t) dmesh_nodal->dn_face;
  PDM_MPI_Allgather((void *) &_dn_face,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&dmesh_nodal->face_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    dmesh_nodal->pdm_mpi_comm);

  // dmesh_nodal->face_distrib[0] = 1;
  dmesh_nodal->face_distrib[0] = 0;

  for (int i = 1; i < dmesh_nodal->n_rank+1; i++) {
    dmesh_nodal->face_distrib[i] +=  dmesh_nodal->face_distrib[i-1];
  }

  if (0 == 1) {
    printf("beg_num_abs::Face : "PDM_FMT_G_NUM" \n", beg_num_abs);
    printf("dmesh_nodal->face_distrib : "PDM_FMT_G_NUM,  dmesh_nodal->face_distrib[0]);
    for (int i = 1; i < dmesh_nodal->n_rank+1; i++) {
      printf(" "PDM_FMT_G_NUM, dmesh_nodal->face_distrib[i]);
    }
    printf("\n");
  }
}

/**
 *
 * \brief Free vtx structure
 *
 * \param[inout]  vtx    Vertices
 *
 * \return        NULL
 *
 */

static
PDM_DMesh_nodal_vtx_t *
_vtx_free
(
 PDM_DMesh_nodal_vtx_t *vtx
)
{
  if (vtx != NULL) {
    if (vtx->distrib != NULL) {
      free (vtx->distrib);
      vtx->distrib = NULL;
    }
    free (vtx);
  }
  return NULL;
}

/**
 *
 * \brief Free a standard section
 *
 * \param [inout]  _bloc_std    Standard section
 *
 * \return         Null
 *
 */

static
void
_section_std_free
(
PDM_DMesh_nodal_section_std_t *_section_std
)
{
  if (_section_std == NULL) {
    return;
  }

  if (_section_std->distrib != NULL) {
    free (_section_std->distrib);
    _section_std->distrib = NULL;
  }
}


/**
 *
 * \brief Free a polygon section
 *
 * \param [inout]  _bloc_poly2d    Polygon section
 *
 * \return         Null
 *
 */

static
void
_section_poly2d_free
(
PDM_DMesh_nodal_section_poly2d_t *_section_poly2d
)
{

  if (_section_poly2d == NULL) {
    return;
  }

  if (_section_poly2d->distrib != NULL) {
    free (_section_poly2d->distrib);
    _section_poly2d->distrib = NULL;
  }

}

/**
 *
 * \brief Free a polyhedron section
 *
 * \param [inout]  _section_poly3d    Polyhedron section
 *
 * \return         Null
 *
 */
static
void
_section_poly3d_free
(
PDM_DMesh_nodal_section_poly3d_t *_section_poly3d
)
{
  if (_section_poly3d == NULL) {
    return;
  }

  if (_section_poly3d->distrib != NULL) {
    free (_section_poly3d->distrib);
    _section_poly3d->distrib = NULL;
  }
}


/**
 *
 * \brief Free a list of standard section
 *
 * \param [inout]  sections    standard sections
 * \param [inout]  n_sections  Number of standard sections
 */
static
void
_sections_std_free
(
 PDM_DMesh_nodal_section_std_t **sections,
 int                             n_sections
)
{
  for(int i_section = 0; i_section < n_sections; ++i_section) {
    _section_std_free(sections[i_section]);
  }
  if(sections != NULL){
    free(sections);
    sections = NULL;
  }
}

/**
 *
 * \brief Free a list of polygon section
 *
 * \param [inout]  sections     Polygon sections
 * \param [inout]  n_sections   Number of polygon sections
 */
static
void
_sections_poly2d_free
(
 PDM_DMesh_nodal_section_poly2d_t **sections,
 int                                n_sections
)
{
  for(int i_section = 0; i_section < n_sections; ++i_section) {
    _section_poly2d_free(sections[i_section]);
  }
  if(sections != NULL){
    free(sections);
    sections = NULL;
  }
}


/**
 *
 * \brief Free a list of polyhedron section
 *
 * \param [inout]  sections    Polyhedron sections
 * \param [inout]  n_sections  Number of polyhedron sections
 */
static
void
_sections_poly3d_free
(
 PDM_DMesh_nodal_section_poly3d_t **sections,
 int                                n_sections
)
{
  for(int i_section = 0; i_section < n_sections; ++i_section) {
    _section_poly3d_free(sections[i_section]);
  }
  if(sections != NULL){
    free(sections);
    sections = NULL;
  }
}

/**
 *
 * \brief Initialize a mesh
 *
 * \param [inout]  mesh        Mesh
 * \param [in]     n_part      Number of partitions
 */

static void
_mesh_init
(
      PDM_dmesh_nodal_t  *dmesh_nodal,
const PDM_MPI_Comm        comm,
      int                 mesh_dimension,
      PDM_g_num_t         n_vtx,
      PDM_g_num_t         n_cell,
      PDM_g_num_t         n_face,
      PDM_g_num_t         n_edge
)
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size (comm, &n_rank);
  PDM_MPI_Comm_rank (comm, &i_rank);

  dmesh_nodal->pdm_mpi_comm             = comm;
  dmesh_nodal->n_rank                   = n_rank;
  dmesh_nodal->i_rank                   = i_rank;

  dmesh_nodal->mesh_dimension           = mesh_dimension;
  dmesh_nodal->n_cell_abs               = n_cell;
  dmesh_nodal->n_face_abs               = n_face;
  dmesh_nodal->n_edge_abs               = n_edge;
  dmesh_nodal->n_vtx_abs                = n_vtx;

  dmesh_nodal->vtx                      = malloc(sizeof(PDM_DMesh_nodal_vtx_t ));
  dmesh_nodal->vtx->_coords             = NULL;
  dmesh_nodal->vtx->distrib             = NULL;
  dmesh_nodal->vtx->n_vtx               = 0;

  dmesh_nodal->n_section_tot            = 0;
  dmesh_nodal->section_type             = NULL;
  dmesh_nodal->section_idx              = NULL;

  dmesh_nodal->n_section                = 0;
  dmesh_nodal->n_section_std            = 0;
  dmesh_nodal->n_section_poly3d         = 0;

  // dmesh_nodal->sections_id              = NULL;
  dmesh_nodal->sections_std             = NULL;
  dmesh_nodal->sections_poly3d          = NULL;

  // dmesh_nodal->sections_id_l1           = NULL;
  dmesh_nodal->n_section_l1             = 0;
  dmesh_nodal->n_section_std_l1         = 0;
  dmesh_nodal->n_section_poly2d_l1      = 0;
  dmesh_nodal->sections_std_l1          = NULL;
  dmesh_nodal->sections_poly2d_l1       = NULL;

  // dmesh_nodal->sections_id_l2           = NULL;
  dmesh_nodal->n_section_l2             = 0;
  dmesh_nodal->n_section_std_l2         = 0;
  dmesh_nodal->sections_std_l2          = NULL;

  dmesh_nodal->section_distribution     = NULL;

  dmesh_nodal->dn_cell                  = -1;
  dmesh_nodal->dcell_face               = NULL;
  dmesh_nodal->dcell_face_idx           = NULL;
  dmesh_nodal->cell_distrib             = NULL;

  dmesh_nodal->dn_face                  = -1;
  dmesh_nodal->_dface_vtx               = NULL;
  dmesh_nodal->_dface_vtx_idx           = NULL;
  dmesh_nodal->_dface_cell              = NULL;
  dmesh_nodal->face_distrib             = NULL;

}

/**
 *
 * \brief Update sections identifier list
 *
 * \param [inout]  mesh        Mesh
 */

// static void
// _update_sections_id
// (
// PDM_dmesh_nodal_t *mesh
// )
// {
//   // TODO : I think this function is useless, we never use section_id
//   PDM_UNUSED(mesh);
//   abort();
  // int n_sections = 0;

  // if (dmesh_nodal->sections_std != NULL) {
  //   n_sections += PDM_Handles_n_get (dmesh_nodal->sections_std);
  // }

  // if (dmesh_nodal->sections_poly2d != NULL) {
  //   n_sections += PDM_Handles_n_get (dmesh_nodal->sections_poly2d);
  // }

  // if (dmesh_nodal->sections_poly3d != NULL) {
  //   n_sections += PDM_Handles_n_get (dmesh_nodal->sections_poly3d);
  // }

  // if (dmesh_nodal->n_section < n_sections) {
  //   dmesh_nodal->sections_id = (int *) realloc(dmesh_nodal->sections_id, sizeof(int) * n_sections);
  // }

  // int k = 0;
  // if (dmesh_nodal->sections_std != NULL) {
  //   const int *id1 = PDM_Handles_idx_get (dmesh_nodal->sections_std);
  //   int n = PDM_Handles_n_get (dmesh_nodal->sections_std);
  //   for (int i = 0; i < n; i++) {
  //     dmesh_nodal->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_STD;
  //   }
  // }

  // if (dmesh_nodal->sections_poly2d != NULL) {
  //   const int *id1 = PDM_Handles_idx_get (dmesh_nodal->sections_poly2d);
  //   int n = PDM_Handles_n_get (dmesh_nodal->sections_poly2d);
  //   for (int i = 0; i < n; i++) {
  //     dmesh_nodal->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_POLY2D;
  //   }
  // }

  // if (dmesh_nodal->sections_poly3d != NULL) {
  //   const int *id1 = PDM_Handles_idx_get (dmesh_nodal->sections_poly3d);
  //   int n = PDM_Handles_n_get (dmesh_nodal->sections_poly3d);
  //   for (int i = 0; i < n; i++) {
  //     dmesh_nodal->sections_id[k++] = id1[i] + PDM_BLOCK_ID_BLOCK_POLY3D;
  //   }
  // }

  // dmesh_nodal->n_section = n_sections;
// }

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

PDM_dmesh_nodal_t*
PDM_DMesh_nodal_create
(
const PDM_MPI_Comm comm,
      int          mesh_dimension,
      PDM_g_num_t  n_vtx,
      PDM_g_num_t  n_cell,
      PDM_g_num_t  n_face,
      PDM_g_num_t  n_edge
)
{
  PDM_dmesh_nodal_t *mesh = (PDM_dmesh_nodal_t *) malloc (sizeof(PDM_dmesh_nodal_t));

  _mesh_init (mesh, comm, mesh_dimension, n_vtx, n_cell, n_face, n_edge);

  return (PDM_dmesh_nodal_t *) mesh;
}


/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  idx      Nodal mesh handle
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 * \return      NULL
 *
 */

void
PDM_DMesh_nodal_free
(
      PDM_dmesh_nodal_t *dmesh_nodal,
const int                partial
)
{
  printf("void PDM_DMesh_nodal_free \n");

  if (dmesh_nodal != NULL) {

    // if (dmesh_nodal->sections_id != NULL) {
    //   free (dmesh_nodal->sections_id);
    // }

    // dmesh_nodal->sections_id = NULL;

     _vtx_free(dmesh_nodal->vtx);

     _sections_std_free   (dmesh_nodal->sections_std   , dmesh_nodal->n_section_std   );
     _sections_poly3d_free(dmesh_nodal->sections_poly3d, dmesh_nodal->n_section_poly3d);

     _sections_std_free   (dmesh_nodal->sections_std_l1   , dmesh_nodal->n_section_std_l1   );
     _sections_poly2d_free(dmesh_nodal->sections_poly2d_l1, dmesh_nodal->n_section_poly2d_l1);

     _sections_std_free   (dmesh_nodal->sections_std_l2   , dmesh_nodal->n_section_std_l2   );


    /* free standard sections */

    // if (dmesh_nodal->sections_std != NULL) {
    //   int n_section_std = PDM_Handles_n_get (dmesh_nodal->sections_std);
    //   const int *list_ind = PDM_Handles_idx_get (dmesh_nodal->sections_std);

    //   while (n_section_std > 0) {
    //     PDM_DMesh_nodal_section_std_t *_bloc_std =
    //       (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (dmesh_nodal->sections_std, list_ind[0]);
    //     _section_std_free(_bloc_std);
    //     PDM_Handles_handle_free (dmesh_nodal->sections_std, list_ind[0], PDM_FALSE);
    //     n_section_std = PDM_Handles_n_get (dmesh_nodal->sections_std);
    //   }

    //   dmesh_nodal->sections_std = PDM_Handles_free (dmesh_nodal->sections_std);
    // }

    // /* Free polygon sections */

    // if (dmesh_nodal->sections_poly2d != NULL) {
    //   int n_section_poly2d = PDM_Handles_n_get (dmesh_nodal->sections_poly2d);
    //   const int *list_ind = PDM_Handles_idx_get (dmesh_nodal->sections_poly2d);

    //   while (n_section_poly2d > 0) {
    //     PDM_DMesh_nodal_section_poly2d_t *_bloc_poly2d =
    //       (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (dmesh_nodal->sections_poly2d, list_ind[0]);
    //     _section_poly2d_free(_bloc_poly2d);
    //     PDM_Handles_handle_free (dmesh_nodal->sections_poly2d, list_ind[0], PDM_FALSE);
    //     n_section_poly2d = PDM_Handles_n_get (dmesh_nodal->sections_poly2d);
    //   }

    //   dmesh_nodal->sections_poly2d = PDM_Handles_free (dmesh_nodal->sections_poly2d);
    // }

    // /* Free polyhedron sections */

    // if (dmesh_nodal->sections_poly3d != NULL) {
    //   int n_section_poly3d = PDM_Handles_n_get (dmesh_nodal->sections_poly3d);
    //   const int *list_ind = PDM_Handles_idx_get (dmesh_nodal->sections_poly3d);

    //   while (n_section_poly3d > 0) {
    //     PDM_DMesh_nodal_section_poly3d_t *_bloc_poly3d =
    //       (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (dmesh_nodal->sections_poly3d, list_ind[0]);
    //     _section_poly3d_free(_bloc_poly3d);
    //     PDM_Handles_handle_free (dmesh_nodal->sections_poly3d, list_ind[0], PDM_FALSE);
    //     n_section_poly3d = PDM_Handles_n_get (dmesh_nodal->sections_poly3d);
    //   }

    //   dmesh_nodal->sections_poly3d = PDM_Handles_free (dmesh_nodal->sections_poly3d);
    // }

    // if (dmesh_nodal->sections_id != NULL) {
    //   free (dmesh_nodal->sections_id);
    // }

    if(partial == 0){
      if (dmesh_nodal->dcell_face_idx != NULL) {
        free (dmesh_nodal->dcell_face_idx);
      }

      if (dmesh_nodal->dcell_face != NULL) {
        free (dmesh_nodal->dcell_face);
      }

      if (dmesh_nodal->cell_distrib != NULL) {
        free (dmesh_nodal->cell_distrib);
      }

      if (dmesh_nodal->_dface_vtx_idx != NULL) {
        free (dmesh_nodal->_dface_vtx_idx);
      }

      if (dmesh_nodal->_dface_vtx != NULL) {
        free (dmesh_nodal->_dface_vtx);
      }

      if (dmesh_nodal->face_distrib != NULL) {
        free (dmesh_nodal->face_distrib);
      }
    }

    free(dmesh_nodal);

  }
}


/**
 * \brief Define partition vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 *
 */

void
PDM_DMesh_nodal_coord_set
(
       PDM_dmesh_nodal_t *dmesh_nodal,
 const int                n_vtx,
 const PDM_real_t        *coords
)
{

  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  if (vtx->_coords != NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Vertices are already defined\n");
  }

  /* Mapping memoire */

  vtx->n_vtx   = n_vtx;
  vtx->_coords = coords;

  vtx->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmesh_nodal->n_rank + 1));

  PDM_g_num_t *_distrib = vtx->distrib + 1;
  _distrib[0] = 0;
  PDM_g_num_t _n_vtx = n_vtx;

  PDM_MPI_Scan (&_n_vtx, _distrib, 1, PDM__PDM_MPI_G_NUM,
                PDM_MPI_SUM, dmesh_nodal->pdm_mpi_comm);


}


/**
 * \brief  Return number of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Number of vertices
 *
 */

int
PDM_DMesh_nodal_n_vtx_get
(
  PDM_dmesh_nodal_t *dmesh_nodal
)
{

  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  return vtx->n_vtx;
}


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Coordinates of vertices
 *
 */

const double *
PDM_DMesh_nodal_vtx_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{

  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  return vtx->_coords;
}


/**
 * \brief  Return number of sections
 *
 * \param [in]  hdl            Nodal mesh handle
 *
 * \return  Number of sections
 *
 */

int
PDM_DMesh_nodal_n_section_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return dmesh_nodal->n_section;
}


/**
 * \brief  Return sections identifier
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 *
 * \return  Blocks identifier
 *
 */
// int *
// PDM_DMesh_nodal_sections_id_get
// (
//PDM_dmesh_nodal_t  *dmesh_nodal
// )
// {
//   PDM_dmesh_nodal_t * mesh =
//           (PDM_dmesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

//   if (dmesh_nodal == NULL) {
//     PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
//   }

//   return dmesh_nodal->sections_id;
// }


/**
 * \brief  Return type of section
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  Type of section
 *
 */

PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_section_type_get
(
PDM_dmesh_nodal_t  *dmesh_nodal,
const int   id_section
)
{
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(id_section);
  abort();

  // PDM_dmesh_nodal_t *mesh =
  //         (PDM_dmesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  PDM_Mesh_nodal_elt_t t_elt = PDM_MESH_NODAL_POLY_3D;

  // if (dmesh_nodal == NULL) {
  //   PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  // }

  // if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

  //   t_elt = PDM_MESH_NODAL_POLY_3D;
  //   PDM_DMesh_nodal_section_std_t *section =
  //           (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (dmesh_nodal->sections_std, id_section);

  //   if (section == NULL) {
  //     PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
  //   }

  //   t_elt = section->t_elt;
  // }

  // else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

  //   t_elt = PDM_MESH_NODAL_POLY_2D;

  // }

  // else {

  //   t_elt = PDM_MESH_NODAL_POLY_3D;

  // }

  return t_elt;
}

/**
 * \brief  Add a new section to the current mesh
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  st_free_data   Status of Release of the memory
 *                             when the section is destroyed
 * \param [in]  id_section       Block identifier
 *
 * \return Block identifier
 *
 */

int
PDM_DMesh_nodal_section_add
(
      PDM_dmesh_nodal_t    *dmesh_nodal,
const PDM_Mesh_nodal_elt_t  t_elt
)
{



  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  // > Le rangement depend de mesh_dimension -> L'ordre importe

  int id_section = dmesh_nodal->n_section_tot++;

  dmesh_nodal->section_type = realloc(dmesh_nodal->section_type, sizeof(PDM_section_type_t) * dmesh_nodal->n_section_tot);
  dmesh_nodal->section_idx  = realloc(dmesh_nodal->section_idx , sizeof(int               ) * dmesh_nodal->n_section_tot);

  switch (t_elt) {

  case PDM_MESH_NODAL_POINT    :
  case PDM_MESH_NODAL_BAR2     :
    {
      dmesh_nodal->section_type[dmesh_nodal->n_section_tot-1] = PDM_SECTION_TYPE_STD1D;
      dmesh_nodal->section_idx [dmesh_nodal->n_section_tot-1] = dmesh_nodal->n_section_std_l2;

      dmesh_nodal->n_section_l2++;
      dmesh_nodal->n_section_std_l2++;

      dmesh_nodal->sections_std_l2 = realloc(dmesh_nodal->sections_std_l2, dmesh_nodal->n_section_std_l2 * sizeof(PDM_DMesh_nodal_section_std_t));
      dmesh_nodal->sections_std_l2[dmesh_nodal->n_section_std_l2-1] = malloc( sizeof(PDM_DMesh_nodal_section_std_t) );
      dmesh_nodal->sections_std_l2[dmesh_nodal->n_section_std_l2-1]->t_elt   = t_elt;
      dmesh_nodal->sections_std_l2[dmesh_nodal->n_section_std_l2-1]->n_elt   = -1;
      dmesh_nodal->sections_std_l2[dmesh_nodal->n_section_std_l2-1]->_connec = NULL;
      dmesh_nodal->sections_std_l2[dmesh_nodal->n_section_std_l2-1]->distrib = NULL;

    }
    break;

  case PDM_MESH_NODAL_TRIA3    :
  case PDM_MESH_NODAL_QUAD4    :
    {
      dmesh_nodal->section_type[dmesh_nodal->n_section_tot-1] = PDM_SECTION_TYPE_STD2D;
      dmesh_nodal->section_idx [dmesh_nodal->n_section_tot-1] = dmesh_nodal->n_section_std_l1;

      dmesh_nodal->n_section_l1++;
      dmesh_nodal->n_section_std_l1++;

      dmesh_nodal->sections_std_l1 = realloc(dmesh_nodal->sections_std_l1, dmesh_nodal->n_section_std_l1 * sizeof(PDM_DMesh_nodal_section_std_t));
      dmesh_nodal->sections_std_l1[dmesh_nodal->n_section_std_l1-1] = malloc( sizeof(PDM_DMesh_nodal_section_std_t) );
      dmesh_nodal->sections_std_l1[dmesh_nodal->n_section_std_l1-1]->t_elt   = t_elt;
      dmesh_nodal->sections_std_l1[dmesh_nodal->n_section_std_l1-1]->n_elt   = -1;
      dmesh_nodal->sections_std_l1[dmesh_nodal->n_section_std_l1-1]->_connec = NULL;
      dmesh_nodal->sections_std_l1[dmesh_nodal->n_section_std_l1-1]->distrib = NULL;

    }
    break;

  case PDM_MESH_NODAL_TETRA4   :
  case PDM_MESH_NODAL_PYRAMID5 :
  case PDM_MESH_NODAL_PRISM6   :
  case PDM_MESH_NODAL_HEXA8    :
    {
      if(dmesh_nodal->mesh_dimension != 3){
        PDM_error(__FILE__, __LINE__, 0, "You cannot specify a 3D standard elements if your meshes not 3D %d\n",
                  dmesh_nodal->mesh_dimension);
        abort();
      }
      dmesh_nodal->section_type[dmesh_nodal->n_section_tot-1] = PDM_SECTION_TYPE_STD3D;
      dmesh_nodal->section_idx [dmesh_nodal->n_section_tot-1] = dmesh_nodal->n_section_std;

      dmesh_nodal->n_section++;
      dmesh_nodal->n_section_std++;

      dmesh_nodal->sections_std = realloc(dmesh_nodal->sections_std, dmesh_nodal->n_section_std * sizeof(PDM_DMesh_nodal_section_std_t));
      dmesh_nodal->sections_std[dmesh_nodal->n_section_std-1] = malloc( sizeof(PDM_DMesh_nodal_section_std_t) );
      dmesh_nodal->sections_std[dmesh_nodal->n_section_std-1]->t_elt   = t_elt;
      dmesh_nodal->sections_std[dmesh_nodal->n_section_std-1]->n_elt   = -1;
      dmesh_nodal->sections_std[dmesh_nodal->n_section_std-1]->_connec = NULL;
      dmesh_nodal->sections_std[dmesh_nodal->n_section_std-1]->distrib = NULL;

    }
    break;

  case PDM_MESH_NODAL_POLY_2D  :
    {
      if(dmesh_nodal->mesh_dimension != 2){
        PDM_error(__FILE__, __LINE__, 0, "You cannot specify a 2D polygon elements if your meshes not 2D %d\n",
                  dmesh_nodal->mesh_dimension);
        abort();
      }
      // Comment faire pour disserner 2D/Volume s?
      dmesh_nodal->section_type[dmesh_nodal->n_section_tot-1] = PDM_SECTION_TYPE_POLY2D;
      dmesh_nodal->section_idx [dmesh_nodal->n_section_tot-1] = dmesh_nodal->n_section_poly2d_l1;

      dmesh_nodal->n_section++;
      dmesh_nodal->n_section_poly2d_l1++;

      dmesh_nodal->sections_poly2d_l1 = realloc(dmesh_nodal->sections_poly2d_l1, dmesh_nodal->n_section_poly2d_l1 * sizeof(PDM_DMesh_nodal_section_poly2d_t));
      dmesh_nodal->sections_poly2d_l1[dmesh_nodal->n_section_poly2d_l1-1] = malloc( sizeof(PDM_DMesh_nodal_section_poly2d_t) );
      dmesh_nodal->sections_poly2d_l1[dmesh_nodal->n_section_poly2d_l1-1]->n_elt       = -1;
      dmesh_nodal->sections_poly2d_l1[dmesh_nodal->n_section_poly2d_l1-1]->_connec     = NULL;
      dmesh_nodal->sections_poly2d_l1[dmesh_nodal->n_section_poly2d_l1-1]->_connec_idx = NULL;
      dmesh_nodal->sections_poly2d_l1[dmesh_nodal->n_section_poly2d_l1-1]->distrib     = NULL;

    }

    break;

  case PDM_MESH_NODAL_POLY_3D  :
    {
      if(dmesh_nodal->mesh_dimension != 3){
        PDM_error(__FILE__, __LINE__, 0, "You cannot specify a 3D polyhedral elements if your meshes not 3D %d\n",
                  dmesh_nodal->mesh_dimension);
        abort();
      }
      assert(dmesh_nodal->mesh_dimension == 3);

      dmesh_nodal->section_type[dmesh_nodal->n_section_tot-1] = PDM_SECTION_TYPE_POLY3D;
      dmesh_nodal->section_idx [dmesh_nodal->n_section_tot-1] = dmesh_nodal->n_section_poly3d;

      dmesh_nodal->n_section++;
      dmesh_nodal->n_section_poly3d++;

      dmesh_nodal->sections_poly3d = realloc(dmesh_nodal->sections_poly3d, dmesh_nodal->n_section_poly3d * sizeof(PDM_DMesh_nodal_section_poly3d_t));
      dmesh_nodal->sections_poly3d[dmesh_nodal->n_section_poly3d-1]                =  malloc( sizeof(PDM_DMesh_nodal_section_poly3d_t) );
      dmesh_nodal->sections_poly3d[dmesh_nodal->n_section_poly3d-1]->n_elt          = -1;
      dmesh_nodal->sections_poly3d[dmesh_nodal->n_section_poly3d-1]->n_face         = -1;
      dmesh_nodal->sections_poly3d[dmesh_nodal->n_section_poly3d-1]->_face_vtx_idx  = NULL;
      dmesh_nodal->sections_poly3d[dmesh_nodal->n_section_poly3d-1]->_face_vtx      = NULL;
      dmesh_nodal->sections_poly3d[dmesh_nodal->n_section_poly3d-1]->_cell_face_idx = NULL;
      dmesh_nodal->sections_poly3d[dmesh_nodal->n_section_poly3d-1]->_cell_face     = NULL;
      dmesh_nodal->sections_poly3d[dmesh_nodal->n_section_poly3d-1]->distrib        = NULL;

    }

    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown element type\n");
    break;

  }

  // _update_sections_id (mesh);
  return id_section ;
}


/**
 * \brief Define a standard section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect        Connectivity
 *
 */
// void
// PDM_DMesh_nodal_section_std_set
// (
//PDM_dmesh_nodal_t  *dmesh_nodal,
// const int          id_section,
// const int          n_elt,
//       PDM_g_num_t *connec,
//       PDM_g_num_t  index_parent_gnum   /* -1 si implicite - else give a shift apply for all element in block
//                                           the numbering inside block behave implicit */
// );

void
PDM_DMesh_nodal_section_std_set
(
PDM_dmesh_nodal_t  *dmesh_nodal,
const int           id_section,
const int           n_elt,
      PDM_g_num_t  *connec
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  PDM_DMesh_nodal_section_std_t* section;
  if(dmesh_nodal->section_type[id_section] == PDM_SECTION_TYPE_STD3D){
    section = dmesh_nodal->sections_std[dmesh_nodal->section_idx[id_section]];
  } else if(dmesh_nodal->section_type[id_section] == PDM_SECTION_TYPE_STD2D){
    section = dmesh_nodal->sections_std_l1[dmesh_nodal->section_idx[id_section]];
  } else if(dmesh_nodal->section_type[id_section] == PDM_SECTION_TYPE_STD1D){
    section = dmesh_nodal->sections_std_l2[dmesh_nodal->section_idx[id_section]];
  }

  PDM_printf(" PDM_Mesh_nodal_section_std_set - _id_section : %i  \n ", _id_section);
  PDM_printf(" PDM_Mesh_nodal_section_std_set - n_elt       : %i  \n ", n_elt);

  /* Mapping */
  section->n_elt   = n_elt;
  section->_connec = connec;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmesh_nodal->n_rank + 1));

  /* Creation of distribution */
  PDM_g_num_t beg_num_abs;
  PDM_g_num_t _n_elt = n_elt;

  PDM_MPI_Scan (&_n_elt, &beg_num_abs, 1, PDM__PDM_MPI_G_NUM,
                PDM_MPI_SUM, dmesh_nodal->pdm_mpi_comm);
  beg_num_abs -= _n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    dmesh_nodal->pdm_mpi_comm);

  section->distrib[0] = 0;
  for (int i = 1; i < dmesh_nodal->n_rank + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }

}


/**
 * \brief Return standard section description
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return  connect           Connectivity
 *
 */

PDM_g_num_t *
PDM_DMesh_nodal_section_std_get
(
      PDM_dmesh_nodal_t *dmesh_nodal,
const int                id_section
)
{
  abort();
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(id_section);
  return NULL;
  // PDM_dmesh_nodal_t *mesh =
  //         (PDM_dmesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  // if (dmesh_nodal == NULL) {
  //   PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  // }

  // int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  // const PDM_DMesh_nodal_section_std_t *section =
  // (const PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (dmesh_nodal->sections_std, _id_section);

  // if (section == NULL) {
  //   PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  // }

  // return section->_connec;
}

/**
 * \brief Get number of section elements
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return      Number of elements
 *
 */

int
PDM_DMesh_nodal_section_n_elt_get
(
      PDM_dmesh_nodal_t *dmesh_nodal,
const int                id_section
)
{
  abort();
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(id_section);
  return -10000;
  // PDM_dmesh_nodal_t *mesh =
  //         (PDM_dmesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  // if (dmesh_nodal == NULL) {
  //   PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  // }

  // int _id_section;

  // if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

  //   _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

  //   const PDM_DMesh_nodal_section_poly3d_t *section =
  //   (const PDM_DMesh_nodal_section_poly3d_t *)
  //    PDM_Handles_get (dmesh_nodal->sections_poly3d, _id_section);

  //   if (section == NULL) {
  //     PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  //   }

  //   return section->n_elt;
  // }

  // else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

  //   _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  //   const PDM_DMesh_nodal_section_poly2d_t *section =
  //   (const PDM_DMesh_nodal_section_poly2d_t *)
  //    PDM_Handles_get (dmesh_nodal->sections_poly2d, _id_section);

  //   if (section == NULL) {
  //     PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
  //   }

  //   return section->n_elt;
  // }

  // else {

  //   _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  //   const PDM_DMesh_nodal_section_std_t *section =
  //   (const PDM_DMesh_nodal_section_std_t *)
  //    PDM_Handles_get (dmesh_nodal->sections_std, _id_section);

  //   if (section == NULL) {
  //     PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
  //   }

  //   return section->n_elt;
  // }

}


/**
 * \brief Define a polygon section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_DMesh_nodal_section_poly2d_set
(
      PDM_dmesh_nodal_t *dmesh_nodal,
const int                id_section,
const PDM_l_num_t        n_elt,
      PDM_l_num_t       *connec_idx,
      PDM_g_num_t       *connec
)
{
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(id_section);
  PDM_UNUSED(n_elt);
  PDM_UNUSED(connec_idx);
  PDM_UNUSED(connec);
  abort();

  // PDM_dmesh_nodal_t *mesh =
  //         (PDM_dmesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  // if (dmesh_nodal == NULL) {
  //   PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  // }

  // int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  // PDM_DMesh_nodal_section_poly2d_t *section =
  //         (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (dmesh_nodal->sections_poly2d, _id_section);

  // if (section == NULL) {
  //   PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  // }

  // /* Mapping */

  // section->n_elt       = n_elt;
  // section->_connec_idx = connec_idx;
  // section->_connec     = connec;

  // section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmesh_nodal->n_rank + 1));

  // /* Creation of distribution */
  // PDM_g_num_t beg_num_abs;
  // PDM_g_num_t _n_elt = n_elt;

  // PDM_MPI_Scan (&_n_elt, &beg_num_abs, 1, PDM__PDM_MPI_G_NUM,
  //               PDM_MPI_SUM, dmesh_nodal->pdm_mpi_comm);
  // beg_num_abs -= _n_elt;

  // PDM_MPI_Allgather((void *) &_n_elt,
  //                   1,
  //                   PDM__PDM_MPI_G_NUM,
  //                   (void *) (&section->distrib[1]),
  //                   1,
  //                   PDM__PDM_MPI_G_NUM,
  //                   dmesh_nodal->pdm_mpi_comm);

  // section->distrib[0] = 0;
  // for (int i = 1; i < dmesh_nodal->n_rank + 1; i++) {
  //   section->distrib[i] +=  section->distrib[i-1];
  // }

}


/**
 * \brief Return a polygon section description
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_DMesh_nodal_section_poly2d_get
(
      PDM_dmesh_nodal_t  *dmesh_nodal,
const int                 id_section,
      PDM_l_num_t       **connec_idx,
      PDM_g_num_t       **connec
)
{
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(id_section);
  // PDM_dmesh_nodal_t *mesh =
  //         (PDM_dmesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

  abort();
  PDM_DMesh_nodal_section_poly2d_t *section = NULL;

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  *connec_idx = section->_connec_idx;
  *connec     = section->_connec;

}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_section_poly3d_set
(
      PDM_dmesh_nodal_t  *dmesh_nodal,
const int                 id_section,
const PDM_l_num_t         n_elt,
const PDM_l_num_t         n_face,
      PDM_l_num_t        *facvtx_idx,
      PDM_g_num_t        *facvtx,
      PDM_l_num_t        *cellfac_idx,
      PDM_g_num_t        *cellfac
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_UNUSED(id_section);
  // int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_DMesh_nodal_section_poly3d_t *section = NULL;
  abort();

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  section->n_elt          = n_elt;
  section->n_face         = n_face;
  section->_face_vtx_idx  = facvtx_idx;
  section->_face_vtx      = facvtx;
  section->_cell_face_idx = cellfac_idx;
  section->_cell_face     = cellfac;

  section->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmesh_nodal->n_rank + 1));

  /* Creation of distribution */
  PDM_g_num_t beg_num_abs;
  PDM_g_num_t _n_elt = n_elt;

  PDM_MPI_Scan (&_n_elt, &beg_num_abs, 1, PDM__PDM_MPI_G_NUM,
                PDM_MPI_SUM, dmesh_nodal->pdm_mpi_comm);
  beg_num_abs -= _n_elt;

  PDM_MPI_Allgather((void *) &_n_elt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&section->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    dmesh_nodal->pdm_mpi_comm);

  section->distrib[0] = 0;
  for (int i = 1; i < dmesh_nodal->n_rank + 1; i++) {
    section->distrib[i] +=  section->distrib[i-1];
  }


}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out]  n_face         Number of faces used to describe polyhedra
 * \param [out]  facvtx_idx     Index of face vertex connectivity
 * \param [out]  facvtx         Face vertex connectivity
 * \param [out]  cellfac_idx    Index of cell face connectivity
 * \param [out]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_section_poly3d_get
(
      PDM_dmesh_nodal_t  *dmesh_nodal,
const int                 id_section,
      PDM_l_num_t        *n_face,
      PDM_l_num_t       **facvtx_idx,
      PDM_g_num_t       **facvtx,
      PDM_l_num_t       **cellfac_idx,
      PDM_g_num_t       **cellfac
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  // int _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

  PDM_DMesh_nodal_section_poly3d_t *section = NULL;
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(id_section);
  abort();

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  *n_face      = section->n_face;
  *facvtx_idx  = section->_face_vtx_idx;
  *facvtx      = section->_face_vtx;
  *cellfac_idx = section->_cell_face_idx;
  *cellfac     = section->_cell_face;

}


/**
 * \brief  Return total number of elements of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return number elements of a partition
 *
 */

PDM_g_num_t
PDM_dmesh_nodal_total_n_cell_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if(dmesh_nodal->mesh_dimension != 3){
    PDM_error(__FILE__, __LINE__, 0, "You cannot compute cell number if you have mesh dimension != 3 PDM_dmesh_nodal_total_n_cell_get %d\n",
              dmesh_nodal->mesh_dimension);
    abort();
  }
  assert(dmesh_nodal->mesh_dimension == 3);

  PDM_g_num_t total_n_cell = 0;
  for(int i_section = 0; i_section < dmesh_nodal->n_section_std; ++i_section){
    total_n_cell += dmesh_nodal->sections_std[i_section]->distrib[dmesh_nodal->n_rank];
  }

  for(int i_section = 0; i_section < dmesh_nodal->n_section_poly3d; ++i_section){
    total_n_cell += dmesh_nodal->sections_poly3d[i_section]->distrib[dmesh_nodal->n_rank];
  }

  return total_n_cell;
}


/**
 * \brief  Return total number of faces of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return total number of faces
 *
 */

PDM_g_num_t
PDM_dmesh_nodal_total_n_face_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  // if (dmesh_nodal->dcell_face == NULL) {
  //   PDM_error (__FILE__, __LINE__, 0, "Not implemented yet\n");
  // }

  return dmesh_nodal->face_distrib[dmesh_nodal->n_rank];

}


/**
 * \brief  Return total number of vertices of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return total number of vertices
 *
 */

PDM_g_num_t
PDM_dmesh_nodal_total_n_vtx_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  return vtx->distrib[dmesh_nodal->n_rank];
}


/**
*
* \brief PDM_sections_decompose_faces
*
* \param [in]     hdl                Distributed nodal mesh handle
* \param [inout]  n_face_elt_tot     Number of faces
* \param [inout]  n_sum_vtx_face_tot Number of vtx for all faces (cumulative)
*
*/
void
PDM_dmesh_nodal_decompose_faces_get_size
(
PDM_dmesh_nodal_t *dmesh_nodal,
int               *n_face_elt_tot,
int               *n_sum_vtx_face_tot
)
{
  /* Get current structure to treat */


  *n_face_elt_tot     = 0;
  *n_sum_vtx_face_tot = 0;

  for (int i_section = 0; i_section < dmesh_nodal->n_section_std; i_section++) {

    int n_face_elt     = PDM_n_face_elt_per_elmt    (dmesh_nodal->sections_std[i_section]->t_elt);
    int n_sum_vtx_face = PDM_n_sum_vtx_face_per_elmt(dmesh_nodal->sections_std[i_section]->t_elt);

    *n_face_elt_tot     += dmesh_nodal->sections_std[i_section]->n_elt*n_face_elt;
    *n_sum_vtx_face_tot += dmesh_nodal->sections_std[i_section]->n_elt*n_sum_vtx_face;

  }

  printf("n_face_elt_tot     ::%i\n", *n_face_elt_tot   );
  printf("n_sum_vtx_face_tot::%i\n" , *n_sum_vtx_face_tot);
}

/**
 *
 * \brief Setup global distribution of all elements register in current structure
 *
 * \param [inout]  mesh
 *
 * \return         Null
 *
 */
void
PDM_dmesh_nodal_generate_distribution
(
 PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  /* Get current structure to treat */


  assert(dmesh_nodal->section_distribution == NULL);

  /* Creation of element distribution among all sections */
  dmesh_nodal->section_distribution    = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmesh_nodal->n_section + 1));
  dmesh_nodal->section_distribution[0] = 0;

  printf("dmesh_nodal->n_section : %i \n", dmesh_nodal->n_section);

  PDM_g_num_t shift = 0;
  for(int i_section = 0; i_section < dmesh_nodal->n_section_std; ++i_section){
    dmesh_nodal->section_distribution[i_section+1] = dmesh_nodal->sections_std[i_section]->distrib[dmesh_nodal->n_rank];
  }
  shift += dmesh_nodal->n_section_std;

  for(int i_section = 0; i_section < dmesh_nodal->n_section_poly3d; ++i_section){
    dmesh_nodal->section_distribution[shift+i_section+1] = dmesh_nodal->sections_poly3d[i_section]->distrib[dmesh_nodal->n_rank];
  }
  // shift += dmesh_nodal->n_section_poly3d;

  for (int i_section = 1; i_section < dmesh_nodal->n_section + 1; i_section++) {
    dmesh_nodal->section_distribution[i_section] +=  dmesh_nodal->section_distribution[i_section-1];
  }

  /* Creation of element distribution among all sections */
  dmesh_nodal->section_distribution_l1    = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmesh_nodal->n_section_l1 + 1));
  dmesh_nodal->section_distribution_l1[0] = 0;
  shift = 0;
  for(int i_section = 0; i_section < dmesh_nodal->n_section_std_l1; ++i_section){
    dmesh_nodal->section_distribution_l1[i_section+1] = dmesh_nodal->sections_std_l1[i_section]->distrib[dmesh_nodal->n_rank];
  }
  shift += dmesh_nodal->n_section_std_l1;
  for(int i_section = 0; i_section < dmesh_nodal->n_section_poly2d_l1; ++i_section){
    dmesh_nodal->section_distribution_l1[shift+i_section+1] = dmesh_nodal->sections_poly2d_l1[i_section]->distrib[dmesh_nodal->n_rank];
  }
  shift += dmesh_nodal->n_section_poly2d_l1;

  for (int i_section = 1; i_section < dmesh_nodal->n_section_l1 + 1; i_section++) {
    dmesh_nodal->section_distribution_l1[i_section] +=  dmesh_nodal->section_distribution_l1[i_section-1];
  }

  /* Creation of element distribution among all sections */
  dmesh_nodal->section_distribution_l2    = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmesh_nodal->n_section_l2 + 1));
  dmesh_nodal->section_distribution_l2[0] = 0;
  shift = 0;
  for(int i_section = 0; i_section < dmesh_nodal->n_section_std_l2; ++i_section){
    dmesh_nodal->section_distribution_l2[i_section+1] = dmesh_nodal->sections_std_l2[i_section]->distrib[dmesh_nodal->n_rank];
  }
  shift += dmesh_nodal->n_section_std_l2;

  for (int i_section = 1; i_section < dmesh_nodal->n_section_l2 + 1; i_section++) {
    dmesh_nodal->section_distribution_l2[i_section] +=  dmesh_nodal->section_distribution_l2[i_section-1];
  }

  // int shift = 0;
  // if (dmesh_nodal->sections_std != NULL) {
  //   int n_section_std  = PDM_Handles_n_get  (dmesh_nodal->sections_std);
  //   const int *list_ind = PDM_Handles_idx_get(dmesh_nodal->sections_std);
  //   for (int i_section = 0; i_section < n_section_std; i_section++) {
  //     PDM_DMesh_nodal_section_std_t *section_std = (PDM_DMesh_nodal_section_std_t *) PDM_Handles_get (dmesh_nodal->sections_std, list_ind[i_section]);
  //     dmesh_nodal->section_distribution[i_section+1] = section_std->distrib[dmesh_nodal->n_rank];
  //   }
  //   shift += n_section_std;
  // }


  // if (dmesh_nodal->sections_poly2d != NULL) {
  //   int n_section_poly2d = PDM_Handles_n_get  (dmesh_nodal->sections_poly2d);
  //   const int *list_ind   = PDM_Handles_idx_get(dmesh_nodal->sections_poly2d);
  //   for (int i_section = 0; i_section < n_section_poly2d; i_section++) {
  //     PDM_DMesh_nodal_section_poly2d_t *section_poly2d = (PDM_DMesh_nodal_section_poly2d_t *) PDM_Handles_get (dmesh_nodal->sections_std, list_ind[i_section]);
  //     dmesh_nodal->section_distribution[i_section+shift+1] = section_poly2d->distrib[dmesh_nodal->n_rank];
  //   }
  //   shift += n_section_poly2d;
  // }

  // if (dmesh_nodal->sections_poly3d != NULL) {
  //   int n_section_poly3d = PDM_Handles_n_get  (dmesh_nodal->sections_poly3d);
  //   const int *list_ind   = PDM_Handles_idx_get(dmesh_nodal->sections_poly3d);
  //   for (int i_section = 0; i_section < n_section_poly3d; i_section++) {
  //     PDM_DMesh_nodal_section_poly3d_t *section_poly3d = (PDM_DMesh_nodal_section_poly3d_t *) PDM_Handles_get (dmesh_nodal->sections_std, list_ind[i_section]);
  //     dmesh_nodal->section_distribution[i_section+shift+1] = section_poly3d->distrib[dmesh_nodal->n_rank];
  //   }
  //   shift += n_section_poly3d;
  // }

  // for (int i_section = 1; i_section < dmesh_nodal->n_section + 1; i_section++) {
  //   dmesh_nodal->section_distribution[i_section] +=  dmesh_nodal->section_distribution[i_section-1];
  // }

  /* Verbose */
  if(1 == 1)
  {
    PDM_printf(" ------------------------------ \n ");
    for(int i_section=0; i_section < dmesh_nodal->n_section+1; i_section++){
      PDM_printf("%i ", dmesh_nodal->section_distribution[i_section]);
    }
  }

}

/**
*
* \brief PDM_dmesh_nodal_decompose_edges_get_size
*
* \param [in]     hdl                Distributed nodal mesh handle
* \param [inout]  n_edge_elt_tot     Number of edges
* \param [inout]  n_sum_vtx_edge_tot Number of vtx for all edges (cumulative)
*
*/
void
PDM_dmesh_nodal_decompose_edges_get_size
(
PDM_dmesh_nodal_t *dmesh_nodal,
int               *n_edge_elt_tot,
int               *n_sum_vtx_edge_tot
)
{

  *n_edge_elt_tot     = 0;
  *n_sum_vtx_edge_tot = 0;

  for (int i_section = 0; i_section < dmesh_nodal->n_section_std; i_section++) {

    int n_edge_elt     = PDM_n_nedge_elt_per_elmt   (dmesh_nodal->sections_std[i_section]->t_elt);
    int n_sum_vtx_edge = PDM_n_sum_vtx_edge_per_elmt(dmesh_nodal->sections_std[i_section]->t_elt);

    *n_edge_elt_tot     += dmesh_nodal->sections_std[i_section]->n_elt*n_edge_elt;
    *n_sum_vtx_edge_tot += dmesh_nodal->sections_std[i_section]->n_elt*n_sum_vtx_edge;

  }

  printf("n_edge_elt_tot     ::%i\n", *n_edge_elt_tot   );
  printf("n_sum_vtx_edge_tot::%i\n" , *n_sum_vtx_edge_tot);
}



/**
 * \brief  Compute cell->cell connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 *
 */
void
PDM_dmesh_nodal_dual_graph
(
PDM_dmesh_nodal_t  *dmesh_nodal,
PDM_g_num_t       **dual_graph_idx,
PDM_g_num_t       **dual_graph,
int                 dim
)
{
  assert(dim == 3); /* dim if for dispatching 2D/3D */
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(dual_graph_idx);
  PDM_UNUSED(dual_graph);

  // Pour Bérenger :
  //     --> On doit concatener les multiple block pour faire un seul cell_vtx
  //         pdm_multi_block_to_part avec parnN = 1 et ln_to_gn = [1, 2, 3, ...., N]
  //     Et on met tout les blocks de données d'entrés :
  //     Une fois ok --> on a le dcell_vtx + dcell_vtx_idx
  //     On calcul la transpose : PDM_dconnectivity_transpose --> dvtx_cell + dvtx_cell_idx
  //     Puis le dual : PDM_deduce_combine_connectivity_dual
  //     See : pdm_dconnectivity_transform.c


}

/**
 * \brief  Compute cell->face connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 *
 */
void
PDM_DMesh_nodal_cell_face_compute2
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{

  PDM_printf("PDM_DMesh_nodal_cell_face_compute2 \n ");

  /* Get current structure to treat */


  // Count the size of the merge connectivity

  // Fonction du choix on génére la bonne table
  // if(GENERATE_EDGE_GLOBAL_NUMBERING)
  // else (GENERATE_FACE_GLOBAL_NUMBERING)

  // cell_vtx : ??? cell_face / face_edge / edge_vtx

  // Option  : (dface_cell oriente)/(ou pas) / cell_face en option
  // if(GENERATE_EDGE_GLOBAL_NUMBERING)
  //    -->
  //    -->
  //    -->
  //    -->
  // else (GENERATE_FACE_GLOBAL_NUMBERING)
  //    --> (option) face_cell
  //    --> face_vtx
  //    --> face_vtx_idx
  //    --> (option) cell_face


  int n_face_elt_tot     = 0;
  int n_sum_vtx_face_tot = 0;

  PDM_dmesh_nodal_decompose_faces_get_size(dmesh_nodal, &n_face_elt_tot, &n_sum_vtx_face_tot);
  // PDM_dmesh_nodal_decompose_edges_get_size(hdl, &n_face_elt_tot, &n_sum_vtx_face_tot);

  // int n_edge_elt_tot     = 0;
  // int n_sum_vtx_edge_tot = 0;

  printf("n_face_elt_tot     ::%i\n", n_face_elt_tot   );
  printf("n_sum_vtx_face_tot::%i\n", n_sum_vtx_face_tot);

  PDM_g_num_t* delmt_face_cell    = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  int*         dcell_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  PDM_g_num_t* dcell_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  dcell_face_vtx_idx[0] = 0;
  PDM_sections_decompose_faces(dmesh_nodal,
                               dcell_face_vtx_idx,
                               dcell_face_vtx,
                               delmt_face_cell,
                               NULL, NULL);
  // PDM_sections_decompose_edges(mesh,
  //                              dcell_face_vtx_idx,
  //                              dcell_face_vtx,
  //                              delmt_face_cell,
  //                              NULL, NULL);

  // Begin the rotation by the plus petit - connectity
  // PDM_connectivity_order_by_min_and_sens (delmt_face_cell, dcell_face_vtx_idx, dcell_face_vtx) // inplace

  /*
   * We are now all information flatten - we only need to compute hash_keys for each faces
   */
  PDM_g_num_t* ln_to_gn = (PDM_g_num_t*) malloc( n_face_elt_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t key_mod = 4 * dmesh_nodal->n_vtx_abs;
  // printf("key_mod::%i \n", key_mod);

  // Option des clés
  _compute_keys(n_face_elt_tot,
                dcell_face_vtx_idx,
                dcell_face_vtx,
                ln_to_gn,
                key_mod);

  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_face_elt_tot , "ln_to_gn:: ");
    PDM_log_trace_array_int(dcell_face_vtx_idx  , n_face_elt_tot , "dcell_face_vtx_idx:: ");
    PDM_log_trace_array_long(dcell_face_vtx, dcell_face_vtx_idx[n_face_elt_tot] , "dcell_face_vtx:: ");
  }

  /*
   * Prepare exchange by computing stride
   */
  int* dcell_face_vtx_n = (int        *) malloc( n_face_elt_tot * sizeof(int        ));
  int* stride_one       = (int        *) malloc( n_face_elt_tot * sizeof(int        ));
  for(int i_face = 0; i_face < n_face_elt_tot; ++i_face) {
    dcell_face_vtx_n[i_face] = dcell_face_vtx_idx[i_face+1] - dcell_face_vtx_idx[i_face];
    stride_one[i_face]       = 1;
  }
  free(dcell_face_vtx_idx);

  /*
   * Setup part_to_block to filter all keys
   */
  double* weight = (double *) malloc( n_face_elt_tot * sizeof(double));
  for(int i = 0; i < n_face_elt_tot; ++i) {
    weight[i] = 1.;
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &ln_to_gn,
                                                      &weight,
                                                      &n_face_elt_tot,
                                                      1,
                                                      dmesh_nodal->pdm_mpi_comm);
  free(weight);
  /*
   * Exchange data
   */
  int*         blk_tot_face_vtx_n = NULL;
  PDM_g_num_t* blk_tot_face_vtx   = NULL;

  int blk_tot_face_vtx_size = PDM_part_to_block_exch(         ptb,
                                                              sizeof(PDM_g_num_t),
                                                              PDM_STRIDE_VAR,
                                                              -1,
                                                              &dcell_face_vtx_n,
                                                    (void **) &dcell_face_vtx,
                                                              &blk_tot_face_vtx_n,
                                                    (void **) &blk_tot_face_vtx);

  int* blk_n_face_per_key = NULL;
  int* blk_face_vtx_n     = NULL;
  int blk_face_vtx_n_size = PDM_part_to_block_exch(         ptb,
                                                            sizeof(int),
                                                            PDM_STRIDE_VAR,
                                                            -1,
                                                            &stride_one,
                                                  (void **) &dcell_face_vtx_n,
                                                            &blk_n_face_per_key,
                                                  (void **) &blk_face_vtx_n);
  free(dcell_face_vtx_n);

  int*         blk_elmt_face_cell_stri = NULL;
  PDM_g_num_t* blk_elmt_face_cell      = NULL;
  int blk_face_cell_size = PDM_part_to_block_exch(         ptb,
                                                           sizeof(PDM_g_num_t),
                                                           PDM_STRIDE_VAR,
                                                           -1,
                                                           &stride_one,
                                                 (void **) &delmt_face_cell,
                                                           &blk_elmt_face_cell_stri,
                                                 (void **) &blk_elmt_face_cell);

  /*
   *  Get the size of the current process bloc
   */
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);

  if( 0 == 1 ) {
    PDM_log_trace_array_int(blk_tot_face_vtx_n, blk_size             , "blk_tot_face_vtx_n:: ");
    PDM_log_trace_array_long(blk_tot_face_vtx , blk_tot_face_vtx_size, "blk_tot_face_vtx:: "  );

    PDM_log_trace_array_int(blk_n_face_per_key, blk_size         , "blk_n_face_per_key:: ");
    PDM_log_trace_array_int(blk_face_vtx_n    , blk_face_vtx_n_size, "blk_face_vtx_n:: ");

    PDM_log_trace_array_long(blk_elmt_face_cell, blk_face_cell_size, "blk_elmt_face_cell:: ");
  }

  PDM_part_to_block_free(ptb);
  free(stride_one);
  free(blk_elmt_face_cell_stri); // Same as blk_n_face_per_key
  free(blk_tot_face_vtx_n);

  /*
   * Get the max number of vertex of faces
   */
  int* blk_face_vtx_idx  = (int        *) malloc( (blk_face_vtx_n_size+1) * sizeof(int        ));
  int n_max_face_per_key = -1;
  for(int i_face = 0; i_face < blk_size; ++i_face) {
    n_max_face_per_key = PDM_MAX(n_max_face_per_key, blk_n_face_per_key[i_face]);
  }

  int n_max_vtx       = -1;
  blk_face_vtx_idx[0] = 0;
  for(int i_face = 0; i_face < blk_face_vtx_n_size; ++i_face) {
    n_max_vtx          = PDM_MAX(n_max_vtx         , blk_face_vtx_n    [i_face]);
    blk_face_vtx_idx[i_face+1] = blk_face_vtx_idx[i_face] + blk_face_vtx_n[i_face];
  }

  // PDM_log_trace_array_int(blk_face_vtx_idx, blk_face_vtx_n_size, "blk_face_vtx_idx:: ");
  /*
   * We need to identify each uniques faces :
   *      - We have multiple packet to treat
   *      - The connectivity can be sorted in place
   *      - Multiple case can occur :
   *           - Alone face normaly boundary
   *           - Multiple faces
   *           - Same face, we remove and replace by the first
   */
  PDM_g_num_t* loc_face_vtx_1 = (PDM_g_num_t *) malloc( n_max_vtx          * sizeof(PDM_g_num_t) );
  PDM_g_num_t* loc_face_vtx_2 = (PDM_g_num_t *) malloc( n_max_vtx          * sizeof(PDM_g_num_t) );
  int*         already_treat  = (int         *) malloc( n_max_face_per_key * sizeof(int        ) );

  /*
   * Allocate Memory - face_vtx - face_cell
   */
  dmesh_nodal->_dface_vtx     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_tot_face_vtx_size  );
  dmesh_nodal->_dface_vtx_idx = (int         *) malloc( sizeof(int        ) * (blk_face_vtx_n_size+1 ));
  dmesh_nodal->_dface_cell    = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_face_cell_size * 2 );

  printf("blk_face_cell_size::%i\n", blk_face_cell_size);

  /*
   * Init global numbering
   */
  int i_abs_face = 0;
  dmesh_nodal->_dface_vtx_idx[0] = 0;

  int idx = 0;
  int idx_face_vtx = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {
    // printf("Number of conflicting keys :: %i \n", blk_n_face_per_key[i_key]);

    int n_conflict_faces = blk_n_face_per_key[i_key];

    /* Reset */
    for(int j = 0; j < n_conflict_faces; ++j) {
      already_treat[j] = -1;
    }

    /* Loop over all faces in conflict */
    for(int i_face = 0; i_face < n_conflict_faces; ++i_face) {
      // printf("Number of vtx on faces %i :: %i with index [%i] \n", i_face, blk_face_vtx_n[idx+i_face], idx+i_face);

      int n_vtx_face_1 = blk_face_vtx_n[idx+i_face];
      int beg_1        = blk_face_vtx_idx[idx+i_face];
      if(already_treat[i_face] != 1) {

        PDM_g_num_t key_1 = 0;
        for(int j = 0; j < n_vtx_face_1; ++j) {
          loc_face_vtx_1[j] = blk_tot_face_vtx[beg_1+j];
          key_1 += loc_face_vtx_1[j];
        }
        PDM_quick_sort_long(loc_face_vtx_1, 0, n_vtx_face_1-1);

        for(int i_face_next = i_face+1; i_face_next < n_conflict_faces; ++i_face_next) {

          int n_vtx_face_2 = blk_face_vtx_n[idx+i_face_next];

          if(n_vtx_face_1 == n_vtx_face_2 ) {

            int beg_2 = blk_face_vtx_idx[idx+i_face_next];
            PDM_g_num_t key_2 = 0;
            for(int j = 0; j < n_vtx_face_1; ++j) {
              loc_face_vtx_2[j] = blk_tot_face_vtx[beg_2+j];
              key_2 += loc_face_vtx_2[j];
            }
            PDM_quick_sort_long(loc_face_vtx_2, 0, n_vtx_face_2-1);

            assert(key_1 == key_2);

            int is_same_face = 1;

            for(int i_vtx = 0; i_vtx < n_vtx_face_1; ++i_vtx) {
              if(loc_face_vtx_1[i_vtx] != loc_face_vtx_2[i_vtx]) {
                is_same_face = -1;
                break;
              }
            }

            if(is_same_face == 1 ){
              // printf(" It's a match ! \n");

              dmesh_nodal->_dface_cell[2*i_abs_face  ] = blk_elmt_face_cell[idx+i_face     ];
              dmesh_nodal->_dface_cell[2*i_abs_face+1] = blk_elmt_face_cell[idx+i_face_next];

              for(int i_vtx = 0; i_vtx < n_vtx_face_1; ++i_vtx) {
                dmesh_nodal->_dface_vtx[idx_face_vtx++] = blk_tot_face_vtx[beg_1+i_vtx];
              }
              dmesh_nodal->_dface_vtx_idx[i_abs_face+1] = dmesh_nodal->_dface_vtx_idx[i_abs_face] + n_vtx_face_1;
              i_abs_face++;

              already_treat[i_face]      = 1;
              already_treat[i_face_next] = 1;
            }


          } /* End if same number of vertex */
        } /* End loop next face */
      } /* End loop already treated */

      /* If a face is not found a match inside the pool, it's a boundary condition */
      if(already_treat[i_face] != 1) {

        dmesh_nodal->_dface_cell[2*i_abs_face  ] = blk_elmt_face_cell[idx+i_face];
        dmesh_nodal->_dface_cell[2*i_abs_face+1] = 0;

        for(int i_vtx = 0; i_vtx < n_vtx_face_1; ++i_vtx) {
          dmesh_nodal->_dface_vtx[idx_face_vtx++] = blk_tot_face_vtx[beg_1+i_vtx];
        }
        dmesh_nodal->_dface_vtx_idx[i_abs_face+1] = dmesh_nodal->_dface_vtx_idx[i_abs_face] + n_vtx_face_1;
        i_abs_face++;

        already_treat[i_face]      = 1;

      }

    } /* End loop face in conflict */

    idx += n_conflict_faces;

  }

  /*
   * Free all unused structure
   */
  free(loc_face_vtx_1);
  free(loc_face_vtx_2);
  free(already_treat);
  free(delmt_face_cell);
  free(dcell_face_vtx);
  free(ln_to_gn);
  free(blk_tot_face_vtx);
  free(blk_n_face_per_key);
  free(blk_face_vtx_n);
  free(blk_elmt_face_cell);

  if( 0 == 1 ){
    printf("i_abs_face::%i \n", i_abs_face);
    PDM_log_trace_array_int(dmesh_nodal->_dface_vtx_idx, i_abs_face+1                    , "dmesh_nodal->_dface_vtx_idx:: ");
    PDM_log_trace_array_long(dmesh_nodal->_dface_vtx   , dmesh_nodal->_dface_vtx_idx[i_abs_face], "_dface_vtx:: ");
    PDM_log_trace_array_long(dmesh_nodal->_dface_cell  , 2*i_abs_face                    , "_dface_cell:: ");
  }

  /*
   * Fill up structure
   */
  dmesh_nodal->dn_face = i_abs_face;

  /*
   * Realloc
   */
  dmesh_nodal->_dface_vtx_idx = (int *        ) realloc((dmesh_nodal->_dface_vtx_idx), (dmesh_nodal->dn_face + 1) * sizeof(int * ) );
  dmesh_nodal->_dface_vtx     = (PDM_g_num_t *) realloc((dmesh_nodal->_dface_vtx    ), dmesh_nodal->_dface_vtx_idx[dmesh_nodal->dn_face] * sizeof(PDM_g_num_t * ));
  dmesh_nodal->_dface_cell    = (PDM_g_num_t *) realloc((dmesh_nodal->_dface_cell   ), dmesh_nodal->dn_face * 2 * sizeof(PDM_g_num_t * ));

  /*
   * Generate absolute numerotation of faces
   */
  _make_absolute_face_numbering(dmesh_nodal);

  /*
   * Rebuild cell face
   */
  int n_face_cell = 2*dmesh_nodal->dn_face;

  int*         part_stri_face_cell = (int         *) malloc( sizeof(int        ) * n_face_cell );
  PDM_g_num_t* dface_cell_tmp      = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * n_face_cell );
  PDM_g_num_t* ln_to_gn_elem       = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) * n_face_cell );

  int idx_g = 0;
  assert( dmesh_nodal->dn_face == dmesh_nodal->face_distrib[dmesh_nodal->i_rank+1] - dmesh_nodal->face_distrib[dmesh_nodal->i_rank]);
  // for (int i = dmesh_nodal->face_distrib[dmesh_nodal->i_rank]; i < dmesh_nodal->face_distrib[dmesh_nodal->i_rank+1]; i++) {
  for (int i_face = 0; i_face < dmesh_nodal->dn_face; ++i_face) {

    PDM_g_num_t g_num_face = (PDM_g_num_t) i_face + dmesh_nodal->face_distrib[dmesh_nodal->i_rank] + 1;

    dface_cell_tmp[idx_g] = dmesh_nodal->_dface_cell[2*i_face];
    ln_to_gn_elem [idx_g] = g_num_face;
    part_stri_face_cell[idx_g] = 1;
    idx_g++;
    if(dmesh_nodal->_dface_cell[2*i_face+1] != 0){
      dface_cell_tmp[idx_g] = dmesh_nodal->_dface_cell[2*i_face+1];
      ln_to_gn_elem [idx_g] = g_num_face;
      part_stri_face_cell[idx_g] = 1;
      idx_g++;
    } else {
      dface_cell_tmp[idx_g] = dmesh_nodal->_dface_cell[2*i_face];
      ln_to_gn_elem [idx_g] = g_num_face;
      part_stri_face_cell[idx_g] = 0;
      idx_g++;
    }
  }
  n_face_cell = idx_g; // Adapt size

  if(0 == 1 ){
    printf("n_face_cell::%i\n", n_face_cell);
    PDM_log_trace_array_int(part_stri_face_cell, n_face_cell, "part_stri_face_cell:: ");
    PDM_log_trace_array_long(ln_to_gn_elem     , n_face_cell, "ln_to_gn_elem:: ");
    PDM_log_trace_array_long(dface_cell_tmp    , n_face_cell, "dface_cell_tmp:: ");
  }

  /*
   *  Use part_to_block with the cell numbering
   */
  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &dface_cell_tmp,
                                                       NULL,
                                                       &n_face_cell,
                                                       1,
                                                       dmesh_nodal->pdm_mpi_comm);

  int         *blk_cell_face_n = NULL;
  PDM_g_num_t *blk_cell_face   = NULL;

  int blk_cell_face_size = PDM_part_to_block_exch(          ptb2,
                                                            sizeof(PDM_g_num_t),
                                                            PDM_STRIDE_VAR,
                                                            -1,
                                                            &part_stri_face_cell,
                                                  (void **) &ln_to_gn_elem,
                                                            &blk_cell_face_n,
                                                  (void **) &blk_cell_face);
  PDM_UNUSED(blk_cell_face_size);

  /*
   *  Get the size of the current process bloc
   */
  int delmt_tot = PDM_part_to_block_n_elt_block_get(ptb2);
  dmesh_nodal->dn_cell = delmt_tot;

  /*
   * Free
   */
  PDM_part_to_block_free(ptb2);
  free(part_stri_face_cell);
  free(dface_cell_tmp     );
  free(ln_to_gn_elem      );

  /*
   * Allcoate
   */
  assert(dmesh_nodal->dcell_face_idx == NULL);
  dmesh_nodal->dcell_face_idx = (int * ) malloc( (delmt_tot + 1) * sizeof(int) );

  dmesh_nodal->dcell_face_idx[0] = 0;
  for(int i = 0; i < delmt_tot; i++){
    dmesh_nodal->dcell_face_idx[i+1] = dmesh_nodal->dcell_face_idx[i] + blk_cell_face_n[i];
  }

  // PDM_log_trace_array_int (blk_cell_face_n, delmt_tot         , "blk_cell_face_n:: ");
  // PDM_log_trace_array_long(blk_cell_face  , blk_cell_face_size, "blk_cell_face:: ");

  assert(dmesh_nodal->dcell_face == NULL);
  dmesh_nodal->dcell_face = blk_cell_face;


  /* Compress connectivity in place */
  PDM_para_graph_compress_connectivity(dmesh_nodal->dn_cell,
                                       dmesh_nodal->dcell_face_idx,
                                       blk_cell_face_n,
                                       dmesh_nodal->dcell_face);


  if( 0 == 1 ){
    printf("dmesh_nodal->dn_cell ::%i\n", dmesh_nodal->dn_cell );
    PDM_log_trace_array_int(dmesh_nodal->dcell_face_idx, dmesh_nodal->dn_cell+1                    , "dmesh_nodal->dcell_face_idx:: ");
    PDM_log_trace_array_long(dmesh_nodal->dcell_face   , dmesh_nodal->dcell_face_idx[dmesh_nodal->dn_cell], "dcell_face:: ");
  }

  /*
   *  Realloc
   */
  dmesh_nodal->dcell_face = (PDM_g_num_t *) realloc( dmesh_nodal->dcell_face, dmesh_nodal->dcell_face_idx[dmesh_nodal->dn_cell] * sizeof(PDM_g_num_t));


  free(blk_cell_face_n);
}

/**
 * \brief  Compute cell->face connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 *
 */
void
PDM_DMesh_nodal_cell_face_compute
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  /* Get current structure to treat */

  if(dmesh_nodal->section_distribution == NULL) {
    PDM_dmesh_nodal_generate_distribution(dmesh_nodal);
  }
  PDM_DMesh_nodal_cell_face_compute2(dmesh_nodal);

  return;
}

/**
 * \brief  Return cell \rightarrow face connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  dcell_face_idx   Index of distributed cell->face connectivity
 * \param [out]  dcell_face       Distributed cell->face connectivity
 *
 * \return     Number of cells on the current process
 *
 */

int
PDM_DMesh_nodal_cell_face_get
(
PDM_dmesh_nodal_t  *dmesh_nodal,
int               **dcell_face_idx,
PDM_g_num_t       **dcell_face
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  *dcell_face_idx = dmesh_nodal->dcell_face_idx;
  *dcell_face = dmesh_nodal->dcell_face;

  return dmesh_nodal->dn_cell;

}

/**
 * \brief  Return face->cell connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  face_cell       Distributed face->cell connectivity
 *
 * \return     Number of cells on the current process
 *
 */
int
PDM_DMesh_nodal_face_cell_get
(
PDM_dmesh_nodal_t  *dmesh_nodal,
PDM_g_num_t       **dface_cell
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  *dface_cell = dmesh_nodal->_dface_cell;

  return dmesh_nodal->dn_face;

}



/**
 * \brief  Return face \rightarrow vertex connectivity
 *
 * \param [in]   hdl              Distributed nodal mesh handle
 * \param [out]  dcell_face_idx   Index of distributed cell->face connectivity
 * \param [out]  dcell_face       Distributed cell->face connectivity
 *
 * \return     Number of faces on the current process
 *
 */

int
PDM_DMesh_nodal_face_vtx_get
(
PDM_dmesh_nodal_t  *dmesh_nodal,
int               **_dface_vtx_idx,
PDM_g_num_t       **_dface_vtx
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  *_dface_vtx_idx = dmesh_nodal->_dface_vtx_idx;
  *_dface_vtx     = dmesh_nodal->_dface_vtx;

  return dmesh_nodal->dn_face;

}

/**
 * \brief  Return vertices distribution
 *
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_ranks + 1
 *
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_vtx_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{


  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  return vtx->distrib;

}

/**
 * \brief  Return section distribution
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  A array of size \ref n_ranks + 1
 *
 */

// const PDM_g_num_t *
// PDM_DMesh_nodal_distrib_section_get
// (
// PDM_dmesh_nodal_t  *dmesh_nodal,
//  const int   id_section
// )
// {
//   PDM_dmesh_nodal_t *mesh =
//           (PDM_dmesh_nodal_t *) PDM_Handles_get (mesh_handles, hdl);

//   if (dmesh_nodal == NULL) {
//     PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
//   }

//   int _id_section;

//   if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

//     _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

//     const PDM_DMesh_nodal_section_poly3d_t *section =
//     (const PDM_DMesh_nodal_section_poly3d_t *)
//      PDM_Handles_get (dmesh_nodal->sections_poly3d, _id_section);

//     if (section == NULL) {
//       PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
//     }

//     return section->distrib;
//   }

//   else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

//     _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

//     const PDM_DMesh_nodal_section_poly2d_t *section =
//     (const PDM_DMesh_nodal_section_poly2d_t *)
//      PDM_Handles_get (dmesh_nodal->sections_poly2d, _id_section);

//     if (section == NULL) {
//       PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
//     }

//     return section->distrib;
//   }

//   else {

//     _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

//     const PDM_DMesh_nodal_section_std_t *section =
//     (const PDM_DMesh_nodal_section_std_t *)
//      PDM_Handles_get (dmesh_nodal->sections_std, _id_section);

//     if (section == NULL) {
//       PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
//     }

//     return section->distrib;
//   }
// }


/**
 * \brief  Return cell distribution
 *
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_ranks + 1
 *
 */

PDM_g_num_t *
PDM_DMesh_nodal_distrib_cell_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{

  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return dmesh_nodal->cell_distrib;

}


/**
 * \brief  Return face distribution
 *
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_ranks + 1
 *
 */

PDM_g_num_t *
PDM_DMesh_nodal_distrib_face_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{

  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return dmesh_nodal->face_distrib;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
