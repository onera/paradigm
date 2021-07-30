
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
#include "pdm_dmesh_nodal_elmts.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_quick_sort.h"
#include "pdm_geom_elem.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dconnectivity_transform.h"

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

    if(vtx->owner == PDM_OWNERSHIP_KEEP) {
      if(vtx->_coords != NULL) {
        free(vtx->_coords);
        vtx->_coords = NULL;
      }
    }
    free (vtx);
  }
  return NULL;
}

static
PDM_dmesh_nodal_elmts_t*
_get_from_geometry_kind
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
)
{
  PDM_dmesh_nodal_elmts_t* dmne = NULL;
  if(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC){
    assert(dmesh_nodal->mesh_dimension == 3);
    dmne = dmesh_nodal->volumic;
  } else if( geom_kind == PDM_GEOMETRY_KIND_SURFACIC){
    assert(dmesh_nodal->mesh_dimension >= 2);
    dmne = dmesh_nodal->surfacic;
  } else if( geom_kind == PDM_GEOMETRY_KIND_RIDGE){
    assert(dmesh_nodal->mesh_dimension >= 1);
    dmne = dmesh_nodal->ridge;
  } else if( geom_kind == PDM_GEOMETRY_KIND_CORNER){
    dmne = dmesh_nodal->corner;
  } else {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom_kind in _get_from_geometry_kind \n");
  }
  assert(dmne != NULL);
  return dmne;
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

  // dmesh_nodal->n_group_elmt             = 0;
  // dmesh_nodal->dgroup_elmt_idx          = NULL;
  // dmesh_nodal->dgroup_elmt              = NULL;
  // dmesh_nodal->dgroup_elmt_owner        = PDM_OWNERSHIP_BAD_VALUE;

  // dmesh_nodal->n_section_tot            = 0;
  // dmesh_nodal->n_section                = 0;
  // dmesh_nodal->n_section_std            = 0;
  // dmesh_nodal->n_section_poly2d         = 0;
  // dmesh_nodal->n_section_poly3d         = 0;

  // dmesh_nodal->sections_id              = NULL;
  // dmesh_nodal->sections_std             = NULL;
  // dmesh_nodal->sections_poly3d          = NULL;
  // dmesh_nodal->sections_poly2d          = NULL;
  // dmesh_nodal->section_distribution     = NULL;

  dmesh_nodal->dn_cell                  = -1;
  dmesh_nodal->dcell_face               = NULL;
  dmesh_nodal->dcell_face_idx           = NULL;
  dmesh_nodal->cell_distrib             = NULL;

  dmesh_nodal->volumic                  = NULL;
  dmesh_nodal->surfacic                 = NULL;
  dmesh_nodal->ridge                    = NULL;
  dmesh_nodal->corner                   = NULL;

  dmesh_nodal->dn_face                  = -1;
  dmesh_nodal->_dface_vtx               = NULL;
  dmesh_nodal->_dface_vtx_idx           = NULL;
  dmesh_nodal->_dface_cell              = NULL;
  dmesh_nodal->face_distrib             = NULL;

}

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

  return mesh;
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

  if (dmesh_nodal != NULL) {

    _vtx_free(dmesh_nodal->vtx);

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

    if(dmesh_nodal->volumic != NULL){
      PDM_DMesh_nodal_elmts_free(dmesh_nodal->volumic);
    }
    if(dmesh_nodal->surfacic != NULL){
      PDM_DMesh_nodal_elmts_free(dmesh_nodal->surfacic);
    }
    if(dmesh_nodal->ridge != NULL){
      PDM_DMesh_nodal_elmts_free(dmesh_nodal->ridge);
    }
    if(dmesh_nodal->corner != NULL){
      PDM_DMesh_nodal_elmts_free(dmesh_nodal->corner);
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
       PDM_real_t        *coords,
       PDM_ownership_t    owner
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
  vtx->owner   = owner;

  vtx->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmesh_nodal->n_rank + 1));

  PDM_g_num_t _n_vtx = n_vtx;
  PDM_MPI_Allgather((void *) &_n_vtx,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&vtx->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    dmesh_nodal->pdm_mpi_comm);

  vtx->distrib[0] = 0;
  for (int i = 1; i < dmesh_nodal->n_rank + 1; i++) {
    vtx->distrib[i] +=  vtx->distrib[i-1];
  }
}


void
PDM_DMesh_nodal_section_g_dims_get
(
  PDM_dmesh_nodal_t *dmesh_nodal,
  PDM_g_num_t       *n_cell_abs,
  PDM_g_num_t       *n_face_abs,
  PDM_g_num_t       *n_edge_abs,
  PDM_g_num_t       *n_vtx_abs
)
{
  *n_cell_abs = dmesh_nodal->n_cell_abs;
  *n_face_abs = dmesh_nodal->n_face_abs;
  *n_edge_abs = dmesh_nodal->n_edge_abs;
  *n_vtx_abs  = dmesh_nodal->n_vtx_abs;
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

double *
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
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  return dmne->n_section;
}

/**
 * \brief  Return sections identifier
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 *
 * \return  Blocks identifier
 *
 */
int *
PDM_DMesh_nodal_sections_id_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  return dmne->sections_id;
}

/**
 * \brief  Return type of element of section
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  Type of section
 *
 */
PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_section_elt_type_get
(
        PDM_dmesh_nodal_t   *dmesh_nodal,
        PDM_geometry_kind_t  geom_kind,
  const int                  id_section
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  if (id_section <= PDM_BLOCK_ID_BLOCK_POLY2D) { // std
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
    return dmne->sections_std[_id_section]->t_elt;
  }
  assert(0); // only useful for std elements
}

PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_section_type_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);

  PDM_Mesh_nodal_elt_t t_elt = PDM_MESH_NODAL_POLY_3D;

  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

    t_elt = PDM_MESH_NODAL_POLY_3D;
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
    PDM_DMesh_nodal_section_std_t *section = dmne->sections_std[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
    }

    t_elt = section->t_elt;
  }

  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

    t_elt = PDM_MESH_NODAL_POLY_2D;

  }

  else {

    t_elt = PDM_MESH_NODAL_POLY_3D;

  }

  return t_elt;
}

/**
 * \brief  Return distri of section
 *
 * \param [in] dmesh_nodal
 * \param [in] id_section   Block identifier
 *
 * \return  distri
 *
 */
PDM_g_num_t*
PDM_DMesh_nodal_section_distri_std_get
(
        PDM_dmesh_nodal_t   *dmesh_nodal,
        PDM_geometry_kind_t  geom_kind,
  const int                  id_section
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  if (id_section <= PDM_BLOCK_ID_BLOCK_POLY2D) { // std
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
    return dmne->sections_std[_id_section]->distrib;
  }
  assert(0); // only useful for std elements
}

int
PDM_DMesh_nodal_section_add
(
      PDM_dmesh_nodal_t    *dmesh_nodal,
      PDM_geometry_kind_t   geom_kind,
const PDM_Mesh_nodal_elt_t  t_elt
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  return PDM_DMesh_nodal_elmts_section_add(dmne, t_elt);
}


void
PDM_DMesh_nodal_section_std_set
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
const int                  n_elt,
      PDM_g_num_t         *connec,
      PDM_ownership_t      owner
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  PDM_DMesh_nodal_elmts_section_std_set(dmne, id_section, n_elt, connec, owner);
}

void
PDM_DMesh_nodal_section_group_elmt_set
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  n_group_elmt,
      int                 *dgroup_elmt_idx,
      PDM_g_num_t         *dgroup_elmt,
      PDM_ownership_t      owner
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  PDM_DMesh_nodal_elmts_group_set(dmne, n_group_elmt, dgroup_elmt_idx, dgroup_elmt, owner);
}

void
PDM_DMesh_nodal_section_group_elmt_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind,
 int                 *n_group_elmt,
 int                 **dgroup_elmt_idx,
 PDM_g_num_t         **dgroup_elmt
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  *n_group_elmt    = dmne->n_group_elmt;
  *dgroup_elmt_idx = dmne->dgroup_elmt_idx;
  *dgroup_elmt     = dmne->dgroup_elmt;
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
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  return PDM_DMesh_nodal_elmts_section_std_get(dmne, id_section);
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
PDM_DMesh_nodal_elmts_section_std_get
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section
)
{
  if (dmn_elts == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  PDM_DMesh_nodal_section_std_t *section = dmn_elts->sections_std[_id_section];

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  return section->_connec;
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
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  return PDM_DMesh_nodal_elmts_section_n_elt_get(dmne, id_section);
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
PDM_DMesh_nodal_elmts_section_n_elt_get
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section
)
{
  int _id_section;

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_DMesh_nodal_section_poly3d_t *section = dmn_elts->sections_poly3d[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
    }

    return section->n_elt;
  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_DMesh_nodal_section_poly2d_t *section = dmn_elts->sections_poly2d[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
    }

    return section->n_elt;
  }

  else {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

    PDM_DMesh_nodal_section_std_t *section = dmn_elts->sections_std[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
    }

    return section->n_elt;
  }

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
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
const PDM_l_num_t          n_elt,
      PDM_l_num_t         *connec_idx,
      PDM_g_num_t         *connec,
      PDM_ownership_t      owner
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  PDM_DMesh_nodal_elmts_section_poly2d_set(dmne, id_section, n_elt, connec_idx, connec, owner);
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
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
      PDM_l_num_t        **connec_idx,
      PDM_g_num_t        **connec
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  PDM_DMesh_nodal_elmts_section_poly2d_get(dmne, id_section, connec_idx, connec);
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
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
const PDM_l_num_t          n_elt,
const PDM_l_num_t          n_face,
      PDM_l_num_t         *facvtx_idx,
      PDM_g_num_t         *facvtx,
      PDM_l_num_t         *cellfac_idx,
      PDM_g_num_t         *cellfac,
      PDM_ownership_t      owner
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  PDM_DMesh_nodal_elmts_section_poly3d_set(dmne, id_section, n_elt, n_face, facvtx_idx, facvtx, cellfac_idx, cellfac, owner);
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
      PDM_dmesh_nodal_t    *dmesh_nodal,
      PDM_geometry_kind_t   geom_kind,
const int                   id_section,
      PDM_l_num_t          *n_face,
      PDM_l_num_t         **facvtx_idx,
      PDM_g_num_t         **facvtx,
      PDM_l_num_t         **cellfac_idx,
      PDM_g_num_t         **cellfac
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  PDM_DMesh_nodal_elmts_section_poly3d_get(dmne, id_section, n_face, facvtx_idx, facvtx, cellfac_idx, cellfac);
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
PDM_dmesh_nodal_total_n_elmt_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  return PDM_DMesh_nodal_elmts_total_n_elmt_get(dmne);
}

/**
 * \brief  Return vtx distribution of a distributed mesh
 *
 * \param [in]  dmesh_nodal
 *
 * \return  Return vtx distribution
 *
 */
PDM_g_num_t*
PDM_dmesh_nodal_vtx_distrib_get
(
  PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  PDM_dmesh_nodal_t* mesh = (PDM_dmesh_nodal_t *) dmesh_nodal;
  return mesh->vtx->distrib;
}

/**
 * \brief  Return vtx distribution of a distributed mesh
 *
 * \param [in]  dmesh_nodal
 *
 * \return  Return vtx distribution
 *
 */
PDM_g_num_t*
PDM_dmesh_nodal_vtx_distrib_copy_get
(
  PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  PDM_dmesh_nodal_t* mesh = (PDM_dmesh_nodal_t *) dmesh_nodal;

  PDM_g_num_t* distrib = (PDM_g_num_t *) malloc( (dmesh_nodal->n_rank + 1) * sizeof(PDM_g_num_t));
  for(int i = 0; i < dmesh_nodal->n_rank + 1; ++i) {
    distrib[i] = mesh->vtx->distrib[i];
  }

  return distrib;
}



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
  if(dmesh_nodal->volumic != NULL) {
     PDM_dmesh_nodal_elmts_generate_distribution(dmesh_nodal->volumic );
  }
  if(dmesh_nodal->surfacic != NULL){
     PDM_dmesh_nodal_elmts_generate_distribution(dmesh_nodal->surfacic);
  }
  if(dmesh_nodal->ridge != NULL){
     PDM_dmesh_nodal_elmts_generate_distribution(dmesh_nodal->ridge);
  }
  if(dmesh_nodal->corner != NULL){
     PDM_dmesh_nodal_elmts_generate_distribution(dmesh_nodal->corner);
  }
}



/**
*
* \brief Concatenates all element sections blocks
*
* \param [in]   dmesh_nodal
* \param [out]  section_idx        index of element section
* \param [out]  cat_delt_vtx_idx   index of element
* \param [out]  cat_delt_vtx       element vtx
*
 * \return     Number sections
*/
int PDM_concat_elt_sections
(
  PDM_dmesh_nodal_t  *dmesh_nodal,
  int               **section_idx,
  int               **cat_delt_vtx_idx,
  PDM_g_num_t       **cat_delt_vtx
)
{
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(section_idx);
  PDM_UNUSED(cat_delt_vtx_idx);
  PDM_UNUSED(cat_delt_vtx);
  abort(); // TODO Remove
  // 0. sizes
  int n_section = 0;
  // int n_section = PDM_DMesh_nodal_n_section_get(dmesh_nodal);
  // int dn_elt_vtx = 0;
  // *section_idx = (int*) malloc((n_section+1) * sizeof(int));
  // int* _section_idx = *section_idx;
  // int* n_vtx_by_elt_by_section = (int*) malloc(n_section * sizeof(int));
  // _section_idx[0] = 0;
  // int* n_elt_vtx_by_section = (int*) malloc(n_section * sizeof(int));
  // for (int i=0; i<n_section; ++i) {
  //   int n_elt_by_section = PDM_DMesh_nodal_section_n_elt_get(dmesh_nodal,i);
  //   _section_idx[i+1] = _section_idx[i] + n_elt_by_section;
  //   PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i);
  //   n_vtx_by_elt_by_section[i] = PDM_Mesh_nodal_n_vertices_element(type,1); // 1: elements of order 1
  //   n_elt_vtx_by_section[i] = n_elt_by_section*n_vtx_by_elt_by_section[i];
  //   dn_elt_vtx += n_elt_vtx_by_section[i];
  // }

  // // 1. cat_delt_vtx_idx
  // int dn_elt = _section_idx[n_section];
  // *cat_delt_vtx_idx = (int*) malloc((dn_elt+1)* sizeof(int));
  // int* _cat_delt_vtx_idx = *cat_delt_vtx_idx;
  // _cat_delt_vtx_idx[0] = 0;
  // int pos_idx = 1;
  // for (int i=0; i<n_section; ++i) {
  //   int n_elt_by_section = _section_idx[i+1] - _section_idx[i];
  //   for (int j=0; j<n_elt_by_section; ++j) {
  //     _cat_delt_vtx_idx[pos_idx+j] = n_vtx_by_elt_by_section[i];
  //   }
  //   pos_idx += n_elt_by_section;
  // }
  // for (int i=1; i<dn_elt+1; ++i) {
  //   _cat_delt_vtx_idx[i] += _cat_delt_vtx_idx[i-1];
  // }

  // // 2. cat_delt_vtx
  // *cat_delt_vtx = (PDM_g_num_t *) malloc(dn_elt_vtx * sizeof(PDM_g_num_t));
  // PDM_g_num_t* _cat_delt_vtx = *cat_delt_vtx;
  // int pos = 0;
  // for (int i=0; i<n_section; ++i) {
  //   PDM_g_num_t* delt_vtx = PDM_DMesh_nodal_section_std_get(dmesh_nodal,i);
  //   for (int j=0; j<n_elt_vtx_by_section[i]; ++j) {
  //     _cat_delt_vtx[pos+j] = delt_vtx[j];
  //   }
  //   pos += n_elt_vtx_by_section[i];
  // }

  // // 3. free
  // free(n_elt_vtx_by_section);
  // free(n_vtx_by_elt_by_section);

  return n_section;
}


/**
*
* \brief Concatenates 3D element sections blocks
*
* \param [in]   dmesh_nodal
* \param [out]  section_idx        index of element section
* \param [out]  cat_delt_vtx_idx   index of element
* \param [out]  cat_delt_vtx       element vtx
*
 * \return     Number sections
*/
int PDM_concat_cell_sections
(
  PDM_dmesh_nodal_t  *dmesh_nodal,
  int               **cell_section_idx,
  int               **cat_dcell_vtx_idx,
  PDM_g_num_t       **cat_dcell_vtx
)
{
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(cell_section_idx);
  PDM_UNUSED(cat_dcell_vtx_idx);
  PDM_UNUSED(cat_dcell_vtx);
  int n_cell_section = -1;
  abort();
  // 0. sizes
  // int n_section = PDM_DMesh_nodal_n_section_get(dmesh_nodal);
  // int n_cell_section = 0;
  // for (int i=0; i<n_section; ++i) {
  //   PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i);
  //   if (PDM_Mesh_nodal_is_3D_element(type)) {
  //     ++n_cell_section;
  //   }
  // }

  // int dn_cell_vtx = 0;
  // *cell_section_idx = (int*) malloc((n_cell_section+1) * sizeof(int));
  // int* _cell_section_idx = *cell_section_idx;
  // int* n_vtx_by_cell_by_section = (int*) malloc(n_cell_section * sizeof(int));
  // _cell_section_idx[0] = 0;
  // int* n_cell_vtx_by_section = (int*) malloc(n_section * sizeof(int));
  // int cnt = 0;
  // for (int i=0; i<n_section; ++i) {
  //   PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i);
  //   if (PDM_Mesh_nodal_is_3D_element(type)) {
  //     int n_cell_by_section = PDM_DMesh_nodal_section_n_elt_get(dmesh_nodal,i);
  //     _cell_section_idx[cnt+1] = _cell_section_idx[cnt] + n_cell_by_section;
  //     n_vtx_by_cell_by_section[cnt] = PDM_Mesh_nodal_n_vertices_element(type,1); // 1: elements of order 1
  //     n_cell_vtx_by_section[i] = n_cell_by_section*n_vtx_by_cell_by_section[cnt];
  //     dn_cell_vtx += n_cell_vtx_by_section[i];
  //     ++cnt;
  //   }
  // }

  // // 1. cat_delt_vtx_idx
  // int dn_cell = _cell_section_idx[n_cell_section];
  // *cat_dcell_vtx_idx = (int*) malloc((dn_cell+1)* sizeof(int));
  // int* _cat_dcell_vtx_idx = *cat_dcell_vtx_idx;
  // _cat_dcell_vtx_idx[0] = 0;
  // int pos_idx = 1;
  // for (int i=0; i<n_cell_section; ++i) {
  //   int n_cell_by_section = _cell_section_idx[i+1] - _cell_section_idx[i];
  //   for (int j=0; j<n_cell_by_section; ++j) {
  //     _cat_dcell_vtx_idx[pos_idx+j] = n_vtx_by_cell_by_section[i];
  //   }
  //   pos_idx += n_cell_by_section;
  // }
  // for (int i=1; i<dn_cell+1; ++i) {
  //   _cat_dcell_vtx_idx[i] += _cat_dcell_vtx_idx[i-1];
  // }

  // // 2. cat_delt_vtx
  // *cat_dcell_vtx = (PDM_g_num_t *) malloc(dn_cell_vtx * sizeof(PDM_g_num_t));
  // PDM_g_num_t* _cat_dcell_vtx = *cat_dcell_vtx;
  // int pos = 0;
  // for (int i=0; i<n_section; ++i) {
  //   PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i);
  //   if (PDM_Mesh_nodal_is_3D_element(type)) {
  //     PDM_g_num_t* dcell_vtx = PDM_DMesh_nodal_section_std_get(dmesh_nodal,i);
  //     for (int j=0; j<n_cell_vtx_by_section[i]; ++j) {
  //       _cat_dcell_vtx[pos+j] = dcell_vtx[j];
  //     }
  //     pos += n_cell_vtx_by_section[i];
  //   }
  // }

  // // 3. free
  // free(n_cell_vtx_by_section);
  // free(n_vtx_by_cell_by_section);

  return n_cell_section;
}


/**
 * \brief  Compute elt->elt connectivity
 *
 * \param [out]  dual_graph_idx
 * \param [out]  dual_graph
 * \param [in]   dim Distributed nodal mesh handle
 *
 */
void
PDM_dmesh_nodal_dual_graph
(
  PDM_g_num_t   *vtx_dist,
  PDM_g_num_t   *elt_dist,
  int           *delt_vtx_idx,
  PDM_g_num_t   *delt_vtx,
  PDM_g_num_t  **delt_elt_idx,
  PDM_g_num_t  **delt_elt,
  PDM_MPI_Comm   comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // 0. transpose
  int* dvtx_elt_idx;
  PDM_g_num_t* dvtx_elt;

  PDM_dconnectivity_transpose(comm,
                              elt_dist, vtx_dist,
                              delt_vtx_idx,delt_vtx,
                              0, // not signed
                              &dvtx_elt_idx,&dvtx_elt);

  // 1. dual
  PDM_deduce_combine_connectivity_dual(comm,
                                       elt_dist, vtx_dist,
                                       delt_vtx_idx,delt_vtx,
                                       dvtx_elt_idx,dvtx_elt,
                                       0, // not signed
                                       delt_elt_idx,
                                       delt_elt);

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
  *dcell_face     = dmesh_nodal->dcell_face;

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

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_section_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  return PDM_DMesh_nodal_elmts_distrib_section_get(dmne, id_section);
}


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


void
PDM_Mesh_nodal_add_dmesh_nodal_elmts
(
 PDM_dmesh_nodal_t       *dmesh_nodal,
 PDM_dmesh_nodal_elmts_t *dmn_elts
)
{
  if(dmn_elts->mesh_dimension == 3) {
    dmesh_nodal->volumic = dmn_elts;
  } else if(dmn_elts->mesh_dimension == 2){
    dmesh_nodal->surfacic = dmn_elts;
  } else if(dmn_elts->mesh_dimension == 1){
    dmesh_nodal->ridge = dmn_elts;
  } else if(dmn_elts->mesh_dimension == 0){
    dmesh_nodal->corner = dmn_elts;
  } else {
    PDM_error (__FILE__, __LINE__, 0, "PDM_Mesh_nodal_add_dmesh_nodal_elmts bad mesh_dimension\n");
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

// /**
// *
// * \brief PDM_dmesh_nodal_decompose_edges_get_size
// *
// * \param [in]     hdl                Distributed nodal mesh handle
// * \param [inout]  n_edge_elt_tot     Number of edges
// * \param [inout]  n_sum_vtx_edge_tot Number of vtx for all edges (cumulative)
// *
// */
// void
// PDM_dmesh_nodal_decompose_edges_get_size
// (
// PDM_dmesh_nodal_t *dmesh_nodal,
// int               *n_edge_elt_tot,
// int               *n_sum_vtx_edge_tot
// )
// {
//   *n_edge_elt_tot     = 0;
//   *n_sum_vtx_edge_tot = 0;

//   for (int i_section = 0; i_section < dmesh_nodal->n_section_std; i_section++) {

//     int n_edge_elt     = PDM_n_nedge_elt_per_elmt   (dmesh_nodal->sections_std[i_section]->t_elt);
//     int n_sum_vtx_edge = PDM_n_sum_vtx_edge_per_elmt(dmesh_nodal->sections_std[i_section]->t_elt);

//     *n_edge_elt_tot     += dmesh_nodal->sections_std[i_section]->n_elt*n_edge_elt;
//     *n_sum_vtx_edge_tot += dmesh_nodal->sections_std[i_section]->n_elt*n_sum_vtx_edge;

//   }

//   for (int i = 0; i < dmesh_nodal->n_section_poly3d; i++) {
//     int _n_face = dmesh_nodal->sections_poly3d[i]->n_face;
//     *n_edge_elt_tot     +=     dmesh_nodal->sections_poly3d[i]->_face_vtx[dmesh_nodal->sections_poly3d[i]->_face_vtx_idx[_n_face]];
//     *n_sum_vtx_edge_tot += 2 * dmesh_nodal->sections_poly3d[i]->_face_vtx[dmesh_nodal->sections_poly3d[i]->_face_vtx_idx[_n_face]];
//   }

//   for (int i = 0; i < dmesh_nodal->n_section_poly2d; i++) {
//     *n_edge_elt_tot     +=     dmesh_nodal->sections_poly2d[i]->_connec_idx[dmesh_nodal->sections_poly2d[i]->n_elt];
//     *n_sum_vtx_edge_tot += 2 * dmesh_nodal->sections_poly2d[i]->_connec_idx[dmesh_nodal->sections_poly2d[i]->n_elt];
//   }

//   // assert(dmesh_nodal->n_section_poly3d == 0); // Not implemented
//   // assert(dmesh_nodal->n_section_poly2d == 0); // Not implemented

//   // printf("n_edge_elt_tot     ::%i\n", *n_edge_elt_tot   );
//   // printf("n_sum_vtx_edge_tot::%i\n" , *n_sum_vtx_edge_tot);
// }
