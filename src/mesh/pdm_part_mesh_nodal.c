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

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_geom_elem.h"
#include "pdm_array.h"
#include "pdm_vtk.h"

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
 * Maximum number of blocks depending of block type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_vtx_free
(
 PDM_Mesh_nodal_vtx_t *vtx
)
{
  if (vtx != NULL) {
    if (vtx->parent != NULL) {
      _vtx_free (vtx->parent);
    }

    if (vtx->coords != NULL) {
      free (vtx->coords);
      vtx->coords = NULL;
    }

  }
}

static
PDM_part_mesh_nodal_elmts_t*
_get_from_geometry_kind
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = NULL;
  if(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC){
    assert(pmn->mesh_dimension == 3);
    pmne = pmn->volumic;
  } else if( geom_kind == PDM_GEOMETRY_KIND_SURFACIC){
    assert(pmn->mesh_dimension >= 2);
    pmne = pmn->surfacic;
  } else if( geom_kind == PDM_GEOMETRY_KIND_RIDGE){
    assert(pmn->mesh_dimension >= 1);
    pmne = pmn->ridge;
  } else if( geom_kind == PDM_GEOMETRY_KIND_CORNER){
    pmne = pmn->corner;
  } else {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom_kind in _get_from_geometry_kind \n");
  }
  return pmne;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Create a \ref PDM_part_mesh_nodal_t structure
 *
 * \param [in]   mesh_dimension   Mesh dimension
 * \param [in]   n_part           Number of partition on the current process
 * \param [in]   comm             MPI communicator
 *
 * \return       Pointer to \ref PDM_part_mesh_nodal_t object
 *
 */

PDM_part_mesh_nodal_t*
PDM_part_mesh_nodal_create
(
 const int          mesh_dimension,
 const int          n_part,
 const PDM_MPI_Comm comm
)
{
  PDM_part_mesh_nodal_t *pmn = (PDM_part_mesh_nodal_t *) malloc (sizeof(PDM_part_mesh_nodal_t));

  pmn->comm           = comm;
  pmn->mesh_dimension = mesh_dimension;
  pmn->n_part         = n_part;

  pmn->vtx      = malloc(n_part * sizeof(PDM_Mesh_nodal_vtx_t *));
  for (int i = 0; i < n_part; i++) {
    pmn->vtx[i] = malloc(sizeof(PDM_Mesh_nodal_vtx_t));
    pmn->vtx[i]->_coords    = NULL;
    pmn->vtx[i]->_numabs    = NULL;
    pmn->vtx[i]->_numparent = NULL;
    pmn->vtx[i]->n_vtx      = 0;
    pmn->vtx[i]->parent     = NULL;
    pmn->vtx[i]->coords     = NULL;
    pmn->vtx[i]->owner      = PDM_OWNERSHIP_KEEP;
  }

  pmn->n_vol    = PDM_array_zeros_int(n_part);
  pmn->n_surf   = PDM_array_zeros_int(n_part);
  pmn->n_ridge  = PDM_array_zeros_int(n_part);
  pmn->n_corner = PDM_array_zeros_int(n_part);

  pmn->volumic  = NULL;
  pmn->surfacic = NULL;
  pmn->ridge    = NULL;
  pmn->corner   = NULL;

  pmn->is_owner_volumic  = PDM_OWNERSHIP_USER;
  pmn->is_owner_surfacic = PDM_OWNERSHIP_USER;
  pmn->is_owner_ridge    = PDM_OWNERSHIP_USER;
  pmn->is_owner_corner   = PDM_OWNERSHIP_USER;

  return pmn;
}


/**
 * \brief Define partition vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part   Partition identifier
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 * \param [in]  numabs    Global numbering
 * \param [in]  owner      Vertices ownship
 *
 */

void
PDM_part_mesh_nodal_coord_set
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part,
 const int                    n_vtx,
 const double                *coords,
 const PDM_g_num_t           *numabs,
       PDM_ownership_t        owner
)
{

  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];

  if ((vtx->_coords != NULL) ||
      (vtx->_numabs != NULL)) {
    PDM_error(__FILE__, __LINE__, 0, "these partition vertices are already defined\n");
  }

  /* Mapping memoire */
  vtx->n_vtx   = n_vtx;
  vtx->_coords = (double *) coords;
  vtx->_numabs = (PDM_g_num_t*) numabs;
  vtx->owner   = owner;

}

/**
 * \brief Define partition vertices from parents
 *
 * \param [in]  pmn           Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part       Partition identifier
 * \param [in]  n_vtx         Number of vertices
 * \param [in]  n_vtx_parent  Number of parent vertices
 * \param [in]  numabs        Global numbering (size = \ref n_vtx)
 * \param [in]  numabs        Global numbering of parent vertices (size = \ref n_vtx_parent)
 * \param [in]  coords        Interlaced coordinates (size = 3 * \ref n_vtx)
 * \param [in]  coords        Interlaced coordinates of parent vertices (size = 3 * \ref n_vtx_parent)
 * \param [in]  owner         Vertices ownship
 *
 */

void
PDM_part_mesh_nodal_coord_from_parent_set
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part,
 const int                    n_vtx,
 const int                    n_vtx_parent,
 const PDM_g_num_t           *numabs,
 const int                   *num_parent,
 const PDM_real_t            *coords_parent,
 const PDM_g_num_t           *numabs_parent,
 const PDM_ownership_t        ownership
)
{

  if (pmn == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];

  if ((vtx->_coords != NULL) ||
      (vtx->_numabs != NULL)) {
    PDM_error(__FILE__, __LINE__, 0, "Vertices are already defined\n");
  }

  vtx->parent = (PDM_Mesh_nodal_vtx_t *) malloc (sizeof (PDM_Mesh_nodal_vtx_t));
  PDM_Mesh_nodal_vtx_t *_parent = vtx->parent;
  _parent->parent = NULL;
  _parent->n_vtx = n_vtx_parent;
  _parent->coords = NULL;
  _parent->_coords = (double *) coords_parent ;
  _parent->_numabs = (PDM_g_num_t *) numabs_parent;
  _parent->_numparent = NULL;
  _parent->owner = PDM_OWNERSHIP_USER;
  _parent->is_coords_get = 0;


  vtx->n_vtx      = n_vtx;
  vtx->coords     = malloc (sizeof(double) * 3 * n_vtx);
  vtx->_coords    = (double *) vtx->coords;
  vtx->_numabs    = (PDM_g_num_t *) numabs;
  vtx->_numparent = (int *) num_parent;
  vtx->owner = ownership;
  vtx->is_coords_get = 0;

  for (int i = 0; i < n_vtx; i++) {
    int i_parent = num_parent[i] - 1;
    for (int j = 0; j < 3; j++) {
      vtx->coords[3*i+j] = _parent->_coords[3*i_parent+j];
    }
  }
  pmn->is_vtx_def_from_parent = 1;

}


/**
 * \brief Add a \ref PDM_part_mesh_nodal_elmts_t to a \ref PDM_part_mesh_nodal_t
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  pmne         Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in]  owner        Ownership
 *
 */

void
PDM_part_mesh_nodal_add_part_mesh_nodal_elmts
(
 PDM_part_mesh_nodal_t       *pmn,
 PDM_part_mesh_nodal_elmts_t *pmne,
 PDM_ownership_t              owner
)
{
  assert(pmn->n_part == pmne->n_part);
  assert(pmn->mesh_dimension >= pmne->mesh_dimension);
  if(pmne->mesh_dimension == 3) {
    pmn->volumic          = pmne;
    pmn->is_owner_volumic = owner;
  } else if(pmne->mesh_dimension == 2){
    pmn->surfacic          = pmne;
    pmn->is_owner_surfacic = owner;
  } else if(pmne->mesh_dimension == 1){
    pmn->ridge          = pmne;
    pmn->is_owner_ridge = owner;
  } else if(pmne->mesh_dimension == 0){
    pmn->corner          = pmne;
    pmn->is_owner_corner = owner;
  } else {
    PDM_error (__FILE__, __LINE__, 0, "PDM_Mesh_nodal_add_dmesh_nodal_elmts bad mesh_dimension\n");
  }
}


/**
 * \brief  Return number of partitions
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return  Number of partitions
 *
 */

int
PDM_part_mesh_nodal_n_part_get
(
       PDM_part_mesh_nodal_t *pmn
)
{
  return pmn->n_part;
}


/**
 * \brief  Return number of vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Number of vertices
 *
 */

int
PDM_part_mesh_nodal_n_vtx_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part
)
{
  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];
  return vtx->n_vtx;
}


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Coordinates of vertices
 *
 */

double*
PDM_part_mesh_nodal_vtx_coord_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part
)
{
  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];
  return (double *) vtx->_coords;
}


/**
 * \brief  Return global ids of vertices
 *
 * \param [in]  pmn       Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  id_part   Partition identifier
 *
 * \return  Golbal ids of vertices
 *
 */

PDM_g_num_t*
PDM_part_mesh_nodal_vtx_g_num_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part
)
{
  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];
  return (PDM_g_num_t*) vtx->_numabs;
}


/**
 * \brief  Return number of sections
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 *
 * \return  Number of sections
 *
 */

int
PDM_part_mesh_nodal_n_section_get
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind
)
{
  if (pmn == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part mesh nodal identifier\n");
  }
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  if(pmne){
    return pmne->n_section;
  } else {
    return 0;
  }
}


/**
 * \brief  Return ids of sections
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 *
 * \return  Ids of sections
 *
 */

int *
PDM_part_mesh_nodal_sections_id_get
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind
)
{
  if (pmn == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part mesh nodal identifier\n");
  }
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  return pmne->sections_id;
}


/**
 * \brief  Return type of section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section   Section identifier
 *
 * \return  Type of section
 *
 */

PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_section_elt_type_get
(
        PDM_part_mesh_nodal_t *pmn,
        PDM_geometry_kind_t    geom_kind,
  const int                    id_section
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_block_type_get(pmne, id_section);
}


/**
 * \brief  Add a new section to the current mesh
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  t_elt        Section type
 *
 * \return Section identifier
 *
 */

int
PDM_part_mesh_nodal_section_add
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const PDM_Mesh_nodal_elt_t   t_elt
)
{
  if( _get_from_geometry_kind(pmn, geom_kind) == NULL) {
    if(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC) {
      pmn->volumic = PDM_part_mesh_nodal_elmts_create(pmn->mesh_dimension, pmn->n_part, pmn->comm);
    } else if( geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {
      pmn->surfacic = PDM_part_mesh_nodal_elmts_create(pmn->mesh_dimension, pmn->n_part, pmn->comm);
    } else if( geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
      pmn->ridge = PDM_part_mesh_nodal_elmts_create(pmn->mesh_dimension, pmn->n_part, pmn->comm);
    } else if( geom_kind == PDM_GEOMETRY_KIND_CORNER) {
      pmn->corner = PDM_part_mesh_nodal_elmts_create(pmn->mesh_dimension, pmn->n_part, pmn->comm);
    }
  }

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_add(pmne, t_elt);
}


/**
 * \brief Define a standard section
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind               Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section              Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [in]  n_elt                   Number of elements
 * \param [in]  connec                  Connectivity
 * \param [in]  numabs                  Global numbering
 * \param [in]  parent_num              Parent numbering or NULL
 * \param [in]  parent_entity_g_num     Parent global numbering or NULL
 * \param [in]  owner                   Ownership
 *
 */

void
PDM_part_mesh_nodal_section_std_set
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const int                    id_section,
const int                    id_part,
const int                    n_elt,
const int                   *connec,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
const PDM_g_num_t           *parent_entity_g_num,
      PDM_ownership_t        owner
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_std_set(pmne, id_section, id_part, n_elt, connec, numabs, parent_num, parent_entity_g_num, owner);
}


/**
 * \brief Define a standard high-order section
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind               Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section              Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [in]  n_elt                   Number of elements
 * \param [in]  connec                  Connectivity
 * \param [in]  numabs                  Global numbering
 * \param [in]  parent_num              Parent numbering or NULL
 * \param [in]  parent_entity_g_num     Parent global numbering or NULL
 * \param [in]  order                   Element order
 * \param [in]  ho_ordering             HO ordering
 * \param [in]  owner                   Ownership
 *
 */

void
PDM_part_mesh_nodal_section_std_ho_set
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const int                    id_section,
const int                    id_part,
const int                    n_elt,
const int                   *connec,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
const PDM_g_num_t           *parent_entity_g_num,
const int                    order,
const char                  *ho_ordering,
      PDM_ownership_t        owner
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_std_ho_set(pmne,
                                       id_section,
                                       id_part,
                                       n_elt,
                                       connec,
                                       numabs,
                                       parent_num,
                                       parent_entity_g_num,
                                       order,
                                       ho_ordering,
                                       owner);
}


/**
 * \brief Get number of section elements
 *
 * \param [in]  pmn        Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind  Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_block   Section identifier
 * \param [in]  id_part    Partition identifier
 *
 * \return      Number of elements
 *
 */

int
PDM_part_mesh_nodal_block_n_elt_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_section,
const int                     id_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_block_n_elt_get(pmne, id_section, id_part);
}


/**
 * \brief Return standard section description
 *
 *  - PDM_MESH_NODAL_POINT :
 *
 *   1 x
 *
 *  - PDM_MESH_NODAL_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_TRIA3 :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_MESH_NODAL_QUAD4 :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_MESH_NODAL_TETRA4 :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_MESH_NODAL_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_MESH_NODAL_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_MESH_NODAL_HEXA8 :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind               Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_block                Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [out] connec                  Connectivity
 * \param [out] numabs                  Global numbering
 * \param [out] parent_num              Parent numbering or NULL
 * \param [out] parent_entity_g_num     Parent global numbering or NULL
 *
 */

void
PDM_part_mesh_nodal_block_std_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_section,
const int                     id_part,
      int                   **connec,
      PDM_g_num_t           **numabs,
      int                   **parent_num,
      PDM_g_num_t           **parent_entity_g_num
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_block_std_get(pmne, id_section, id_part, connec, numabs, parent_num, parent_entity_g_num);
}


/**
 * \brief Return standard high-order section description
 *
 * \param [in]  pmn                     Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind               Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section              Section identifier
 * \param [in]  id_part                 Partition identifier
 * \param [out] connec                  Connectivity
 * \param [out] numabs                  Global numbering
 * \param [out] parent_num              Parent numbering or NULL
 * \param [out] parent_entity_g_num     Parent global numbering or NULL
 * \param [out] order                   Element order
 * \param [out] ho_ordering             HO ordering
 *
 */

void
PDM_part_mesh_nodal_section_std_ho_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_section,
const int                     id_part,
      int                   **connec,
      PDM_g_num_t           **numabs,
      int                   **parent_num,
      PDM_g_num_t           **parent_entity_g_num,
      int                    *order,
const char                  **ho_ordering
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_block_std_ho_get(pmne,
                                             id_section,
                                             id_part,
                                             connec,
                                             numabs,
                                             parent_num,
                                             parent_entity_g_num,
                                             order,
                                             ho_ordering);
}


/**
 * \brief Get parent numbering of block elements
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_block     Block identifier
 * \param [in]  id_part      Partition identifier
 *
 * \return      Return parent numbering of block elements
 *
 */

int *
PDM_part_mesh_nodal_block_parent_num_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_section,
const int                     id_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_parent_num_get(pmne, id_section, id_part);
}


/**
 * \brief Get global element numbering of block elements
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_block     Block identifier
 * \param [in]  id_part      Partition identifier
 *
 * \return      Return global element numbering of block elements
 *
 */
 // ATTENTION != PDM_Mesh_nodal_block_g_num_get

PDM_g_num_t *
PDM_part_mesh_nodal_block_g_num_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_section,
const int                     id_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_g_num_get(pmne, id_section, id_part);
}


/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 *
 */

void
PDM_part_mesh_nodal_free
(
 PDM_part_mesh_nodal_t* pmn
)
{

  if(pmn->is_owner_volumic == PDM_OWNERSHIP_KEEP) {
    PDM_part_mesh_nodal_elmts_free(pmn->volumic);
  }
  if(pmn->is_owner_surfacic == PDM_OWNERSHIP_KEEP) {
    PDM_part_mesh_nodal_elmts_free(pmn->surfacic);
  }
  if(pmn->is_owner_ridge == PDM_OWNERSHIP_KEEP) {
    PDM_part_mesh_nodal_elmts_free(pmn->ridge);
  }
  if(pmn->is_owner_corner == PDM_OWNERSHIP_KEEP) {
    PDM_part_mesh_nodal_elmts_free(pmn->corner);
  }

  if (pmn->vtx != NULL) {
    for (int i = 0; i < pmn->n_part; i++) {
      if(pmn->vtx[i]->owner == PDM_OWNERSHIP_KEEP){
        _vtx_free (pmn->vtx[i]);
      }
      free(pmn->vtx[i]);
    }

    free(pmn->vtx);
    pmn->vtx = NULL;
  }

  free(pmn->n_vol   );
  free(pmn->n_surf  );
  free(pmn->n_ridge );
  free(pmn->n_corner);

  free(pmn);
}


void
PDM_part_mesh_nodal_dump_vtk
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind,
 const char            *filename_patter
)
{
  int i_rank = -1;
  PDM_MPI_Comm_rank(pmn->comm, &i_rank);

  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int pn_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    double      *pvtx_coord    = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);
    PDM_g_num_t *pvtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part);

    int  n_section  = PDM_part_mesh_nodal_n_section_get  (pmn, geom_kind);
    int *section_id = PDM_part_mesh_nodal_sections_id_get(pmn, geom_kind);

    // printf("pn_vtx = %i\n", pn_vtx);

    for(int i_section = 0; i_section < n_section; ++i_section) {

      int id_section = section_id  [i_section];
      int                  n_elt = PDM_part_mesh_nodal_elmts_block_n_elt_get(pmne, id_section, i_part);
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_block_type_get (pmne, id_section);

      char filename[999];
      sprintf(filename, "%s_section_%2.2d_%2.2d_%2.2d.vtk", filename_patter, i_section, i_part, i_rank);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        int *connec_idx = NULL;
        int *connec     = NULL;
        PDM_part_mesh_nodal_elmts_block_poly2d_get(pmne,
                                                   section_id[i_section],
                                                   i_part,
                                                   &connec_idx,
                                                   &connec);


        PDM_g_num_t *pelmt_ln_to_gn = PDM_part_mesh_nodal_elmts_g_num_get(pmne,
                                                                          section_id[i_section],
                                                                          i_part);

        PDM_vtk_write_polydata(filename,
                               pn_vtx,
                               pvtx_coord,
                               pvtx_ln_to_gn,
                               n_elt,
                               connec_idx,
                               connec,
                               pelmt_ln_to_gn,
                               NULL);
      }

      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        abort();
      }

      else {
        int is_ho = PDM_Mesh_nodal_elmt_is_ho(t_elt);
        if(is_ho) {
          int order;
          const char *ho_ordering      = NULL;
          int         *pcell_vtx       = NULL;
          PDM_g_num_t *pelmt_ln_to_gn  = NULL;
          int         *parent_num      = NULL;
          PDM_g_num_t *parent_elmt_num = NULL;
          PDM_part_mesh_nodal_elmts_block_std_ho_get(pmne,
                                                     section_id[i_section],
                                                     i_part,
                                                     &pcell_vtx,
                                                     &pelmt_ln_to_gn,
                                                     &parent_num,
                                                     &parent_elmt_num,
                                                     &order,
                                                     &ho_ordering);

          int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
          int *pcell_vtx_out = malloc(n_vtx_per_elmt * n_elt * sizeof(int));
          for(int i = 0; i < n_vtx_per_elmt * n_elt; ++i) {
            pcell_vtx_out[i] = pcell_vtx[i];
          }

          PDM_Mesh_nodal_reorder_elt_vtx(t_elt,
                                         order,
                                         ho_ordering,
                                         "PDM_HO_ORDERING_VTK",
                                         n_elt,
                                         pcell_vtx,
                                         pcell_vtx_out);


          PDM_vtk_write_std_elements_ho(filename,
                                        order,
                                        pn_vtx,
                                        pvtx_coord,
                                        pvtx_ln_to_gn,
                                        t_elt,
                                        n_elt,
                                        pcell_vtx_out,
                                        pelmt_ln_to_gn,
                                        0,
                                        NULL,
                                        NULL);
          free(pcell_vtx_out);
        } else {

          int         *pcell_vtx       = NULL;
          PDM_g_num_t *pelmt_ln_to_gn  = NULL;
          int         *parent_num      = NULL;
          PDM_g_num_t *parent_elmt_num = NULL;
          PDM_part_mesh_nodal_elmts_block_std_get(pmne,
                                                  section_id[i_section],
                                                  i_part,
                                                  &pcell_vtx,
                                                  &pelmt_ln_to_gn,
                                                  &parent_num,
                                                  &parent_elmt_num);

          PDM_vtk_write_std_elements(filename,
                                     pn_vtx,
                                     pvtx_coord,
                                     pvtx_ln_to_gn,
                                     t_elt,
                                     n_elt,
                                     pcell_vtx,
                                     pelmt_ln_to_gn,
                                     0,
                                     NULL,
                                     NULL);
        }
      }
    }
  }
}


/**
 * \brief Compute element extents of a part of a section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind      Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  tolerance      Expansion tolerance for bounding boxes
 * \param [out] extents        Extents of mesh elements in current part of current block
 *
 */

void
PDM_part_mesh_nodal_block_elt_extents_compute
(
       PDM_part_mesh_nodal_t *pmn,
       PDM_geometry_kind_t    geom_kind,
 const int                    id_section,
 const int                    i_part,
 const double                 tolerance,
       double                *extents
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);

  PDM_part_mesh_nodal_elmts_elt_extents_compute(pmne,
                                                id_section,
                                                i_part,
                                                tolerance,
                                                vtx_coord,
                                                extents);
}


/**
 * \brief Compute cell centers of a part of section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind      Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  ownership      Ownership
 *
 */

void
PDM_part_mesh_nodal_section_elt_center_compute
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const int                    id_section,
const int                    i_part,
const PDM_ownership_t        ownership
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);

  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);

  PDM_part_mesh_nodal_elmts_elt_center_compute(pmne,
                                               id_section,
                                               i_part,
                                               n_vtx,
                                               vtx_coord,
                                               ownership);
}


/**
 * \brief  Return cell centers
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind      Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 *
 * \return  Return cell centers
 *
 */

const double *
PDM_part_mesh_nodal_section_elt_center_get
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const int                    id_section,
const int                    i_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  return PDM_part_mesh_nodal_elmts_elt_center_get(pmne, id_section, i_part);
}


/**
 * \brief Reset cell center computation
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind      Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 *
 */

void
PDM_part_mesh_nodal_section_elt_center_reset
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const int                    id_section,
const int                    i_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  PDM_part_mesh_nodal_elmts_elt_center_reset(pmne, id_section, i_part);
}

/**
 * \brief Define a polygon section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind      Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connec_idx     Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connec         Connectivity
 * \param [in]  numabs         Global numbering
 * \param [in]  parent_num     Parent numbering or NULL
 * \param [in]  owner          Ownership
 *
 */

void
PDM_part_mesh_nodal_section_poly2d_set
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const int                    id_section,
const int                    id_part,
const int                    n_elt,
const int                   *connec_idx,
const int                   *connec,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
      PDM_ownership_t        owner
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_block_poly2d_set(pmne,
                                             id_section,
                                             id_part,
                                             n_elt,
                                             connec_idx,
                                             connec,
                                             numabs,
                                             parent_num,
                                             owner);
}


/**
 * \brief Return a polygon section description
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind      Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] connec_idx     Connectivity index
 * \param [out] connec         Connectivity
 *
 */

void
PDM_part_mesh_nodal_block_poly2d_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_section,
const int                     id_part,
      int                   **connec_idx,
      int                   **connec
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_block_poly2d_get(pmne,
                                             id_section,
                                             id_part,
                                             connec_idx,
                                             connec);
}


/**
 * \brief Define a polyhedron section
 *
 * \param [in]  pmn            Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind      Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section     Section identifier
 * \param [in]  id_part        Partition identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  n_face         Number of faces
 * \param [in]  facvtx_idx     Face->vertex connectivity index (size = \ref n_face + 1)
 * \param [in]  facvtx         Face->vertex connectivity
 * \param [in]  face_ln_to_gn  Face global ids
 * \param [in]  cellfac_idx    Cell->face connectivity index (size = \ref n_cell + 1)
 * \param [in]  cellfac        Cell->face connectivity
 * \param [in]  numabs         Cell global ids
 * \param [in]  parent_num     Parent numbering or NULL
 * \param [in]  owner          Ownership
 *
 */

void
PDM_part_mesh_nodal_section_poly3d_set
(
      PDM_part_mesh_nodal_t *pmn,
      PDM_geometry_kind_t    geom_kind,
const int                    id_section,
const int                    id_part,
const int                    n_elt,
const int                    n_face,
const int                   *facvtx_idx,
const int                   *facvtx,
const PDM_g_num_t           *face_ln_to_gn,
const int                   *cellfac_idx,
const int                   *cellfac,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
      PDM_ownership_t        owner
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_block_poly3d_set(pmne,
                                             id_section,
                                             id_part,
                                             n_elt,
                                             n_face,
                                             facvtx_idx,
                                             facvtx,
                                             face_ln_to_gn,
                                             cellfac_idx,
                                             cellfac,
                                             numabs,
                                             parent_num,
                                             owner);
}


/**
 * \brief Return a polyhedron section
 *
 * \param [in]  pmn                  Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind            Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section           Section identifier
 * \param [in]  id_part              Partition identifier
 * \param [out] n_face               Number of faces
 * \param [out] face_ln_to_gn        Face global ids
 * \param [out] facvtx_idx           Face->vertex connectivity index (size = \ref n_face + 1)
 * \param [out] facvtx               Face->vertex connectivity
 * \param [out] numabs               Cell global ids
 * \param [out] cell_face_idx        Cell->face connectivity index (size = \ref n_cell + 1)
 * \param [out] cell_face            Cell->face connectivity
 * \param [out] parent_num           Parent numbering or NULL
 * \param [out] parent_entity_g_num  Parent global ids or NULL
 *
 */

void
PDM_part_mesh_nodal_section_poly3d_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_section,
const int                     id_part,
      int                    *n_face,
      PDM_g_num_t           **face_ln_to_gn,
      int                   **face_vtx_idx,
      int                   **face_vtx,
      PDM_g_num_t           **numabs,
      int                   **cell_face_idx,
      int                   **cell_face,
      int                   **parent_num,
      PDM_g_num_t           **parent_entity_g_num
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_block_poly3d_get(pmne,
                                             id_section,
                                             id_part,
                                             n_face,
                                             face_ln_to_gn,
                                             face_vtx_idx,
                                             face_vtx,
                                             numabs,
                                             cell_face_idx,
                                             cell_face,
                                             parent_num,
                                             parent_entity_g_num);
}


/**
 * \brief Get the cell-vertex connectivity of a polyhedron section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section   Section identifier
 * \param [in]  id_part      Partition identifier
 * \param [out] cellvtx_idx  Index of cell vertex connectivity
 * \param [out] cellvtx      Cell vertex connectivity
 *
 */

void
PDM_part_mesh_nodal_section_poly3d_cell_vtx_connect_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_section,
const int                     id_part,
      int                   **cellvtx_idx,
      int                   **cellvtx
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_block_poly3d_cell_vtx_connect_get(pmne,
                                                              id_section,
                                                              id_part,
                                                              cellvtx_idx,
                                                              cellvtx);
}


/**
 * \brief Reset a nodal mesh structure
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 *
 * \return      NULL
 *
 */

void
PDM_part_mesh_nodal_reset
(
 PDM_part_mesh_nodal_t *pmn
)
{
  for (PDM_geometry_kind_t geom_kind = 0; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {
    PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);

    if (pmne != NULL) {
      PDM_part_mesh_nodal_elmts_reset(pmne);
    }

  }

  if (pmn->vtx != NULL) {
    for (int i = 0; i < pmn->n_part; i++) {
      pmn->vtx[i]->_coords = NULL;
      pmn->vtx[i]->_numabs = NULL;
      pmn->vtx[i]->_numparent = NULL;
      pmn->vtx[i]->n_vtx   = 0;
      if (pmn->vtx[i]->parent != NULL) {
        _vtx_free (pmn->vtx[i]->parent);
        pmn->vtx[i]->parent = NULL;
      }
      if (pmn->vtx[i]->coords != NULL) {
        free (pmn->vtx[i]->coords);
        pmn->vtx[i]->coords = NULL;
      }
    }
  }
}


/**
 * \brief  Compute a global numbering in a section
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section   Section identifier
 * \param [in]  ownership    Ownership
 *
 */
//---> PDM_Mesh_nodal_g_num_in_block_compute

void
PDM_part_mesh_nodal_g_num_in_section_compute
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_section,
const PDM_ownership_t         ownership
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  PDM_part_mesh_nodal_elmts_g_num_in_section_compute(pmne,
                                                     id_section,
                                                     ownership);
}


/**
 * \brief  Return number elements of a partition
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_part      Partition identifier
 *
 * \return  Return number elements of a partition
 *
 */
// ---> PDM_Mesh_nodal_n_cell_get

int
PDM_part_mesh_nodal_n_elmts_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  return PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, id_part);
}


/**
 * \brief Get the element global numbering taking into account parent_num
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_part      Partition identifier
 *
 * \return  Global ids of element in current partition
 *
 */
// ---> PDM_Mesh_nodal_g_num_get_from_part

PDM_g_num_t *
PDM_part_mesh_nodal_g_num_get_from_part
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  return PDM_part_mesh_nodal_elmts_g_num_get_from_part(pmne, id_part);
}


/**
 * \brief Free partially a part_mesh_nodal structure
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return      NULL
 *
 */
//---> PDM_Mesh_nodal_partial_free

void
PDM_part_mesh_nodal_partial_free
(
 PDM_part_mesh_nodal_t *pmn
)
{
  for (PDM_geometry_kind_t geom_kind = 0; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {
    PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);

    if (pmne != NULL) {
      PDM_part_mesh_nodal_elmts_partial_free(pmne);
    }

  }
}


/**
 * \brief Extract vertices from parent vertices
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 *
 * \return true if the vertices are defined from parents
 */
//---> PDM_Mesh_nodal_is_set_coord_from_parent

int
PDM_part_mesh_nodal_is_set_coord_from_parent
(
 PDM_part_mesh_nodal_t *pmn
)
{
  if (pmn == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  return pmn->is_vtx_def_from_parent;
}


/**
 * \brief Get global element numbering of block elements inside the block
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_section   Section identifier
 * \param [in]  id_part      Partition identifier
 *
 * \return      Return global numbering of block elements inside the block
 *
 */
//---> PDM_Mesh_nodal_block_g_num_get

PDM_g_num_t *
PDM_part_mesh_nodal_section_g_num_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_section,
const int                     id_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  return PDM_part_mesh_nodal_elmts_section_g_num_get(pmne, id_section, id_part);
}


/**
 * \brief  Return parent element number to local number
 *
 * \param [in]  pmn          Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind    Geometry kind (corner, ridge, surface or volume)
 * \param [in]  id_part      Partition identifier
 *
 * \return  Parent element number to local number
 *
 */
//---> PDM_Mesh_nodal_num_cell_parent_to_local_get

int *
PDM_part_mesh_nodal_num_elmt_parent_to_local_get
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  return PDM_part_mesh_nodal_elmts_num_elmt_parent_to_local_get(pmne, id_part);
}
