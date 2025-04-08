/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm.h"
#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_mem_tool.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal_elmts_priv.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_unique.h"
#include "pdm_logging.h"
#include "pdm_gnum.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_order.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_part_comm_graph.h"

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
      _vtx_free(vtx->parent);
      vtx->parent = NULL;
    }

    if (vtx->_coords != NULL && vtx->owner_coords == PDM_OWNERSHIP_KEEP) {
      PDM_free(vtx->_coords);
    }

    if (vtx->_numabs != NULL && vtx->owner_numabs == PDM_OWNERSHIP_KEEP) {
      PDM_free(vtx->_numabs);
    }

    if (vtx->_numparent != NULL && vtx->owner_numparent == PDM_OWNERSHIP_KEEP) {
      PDM_free(vtx->_numparent);
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

PDM_part_mesh_nodal_t*
PDM_part_mesh_nodal_create
(
 const int          mesh_dimension,
 const int          n_part,
 const PDM_MPI_Comm comm
)
{
  PDM_part_mesh_nodal_t *pmn;
  PDM_malloc(pmn, 1, PDM_part_mesh_nodal_t);
  memset(pmn, 0, sizeof(PDM_part_mesh_nodal_t));

  pmn->comm           = comm;
  pmn->mesh_dimension = mesh_dimension;
  pmn->n_part         = n_part;

  PDM_malloc(pmn->vtx, n_part, PDM_Mesh_nodal_vtx_t *);
  for (int i = 0; i < n_part; i++) {
    PDM_malloc(pmn->vtx[i], 1, PDM_Mesh_nodal_vtx_t);
    pmn->vtx[i]->_coords         = NULL;
    pmn->vtx[i]->_numabs         = NULL;
    pmn->vtx[i]->_numparent      = NULL;
    pmn->vtx[i]->n_vtx           = 0;
    pmn->vtx[i]->parent          = NULL;
    pmn->vtx[i]->coords          = NULL;
    pmn->vtx[i]->owner_coords    = PDM_OWNERSHIP_KEEP;
    pmn->vtx[i]->owner_numabs    = PDM_OWNERSHIP_KEEP;
    pmn->vtx[i]->owner_numparent = PDM_OWNERSHIP_KEEP;
  }

  pmn->s_section = 10;
  pmn->n_section = 0;
  PDM_malloc(pmn->section_kind, pmn->s_section, PDM_geometry_kind_t);
  PDM_malloc(pmn->section_id,   pmn->s_section, int                );

  return pmn;
}

void
PDM_part_mesh_nodal_coord_set
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part,
 const int                    n_vtx,
 const double                *coords,
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
  vtx->n_vtx        = n_vtx;
  vtx->_coords      = (double *) coords;
  vtx->owner_coords = owner;

}


void
PDM_part_mesh_nodal_vtx_gnum_set
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part,
 const PDM_g_num_t           *numabs,
       PDM_ownership_t        owner
)
{

  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];

  if ((vtx->_numabs != NULL)) {
    PDM_error(__FILE__, __LINE__, 0, "these partition vertices are already defined\n");
  }

  vtx->_numabs      = (PDM_g_num_t*) numabs;
  vtx->owner_numabs = owner;
}




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

  PDM_malloc(vtx->parent,1,PDM_Mesh_nodal_vtx_t);
  PDM_Mesh_nodal_vtx_t *_parent = vtx->parent;
  _parent->parent          = NULL;
  _parent->n_vtx           = n_vtx_parent;
  _parent->coords          = NULL;
  _parent->_coords         = (double      *) coords_parent;
  _parent->_numabs         = (PDM_g_num_t *) numabs_parent;
  _parent->_numparent      = NULL;
  _parent->owner_coords    = PDM_OWNERSHIP_USER;
  _parent->owner_numabs    = PDM_OWNERSHIP_USER;
  _parent->owner_numparent = PDM_OWNERSHIP_USER;


  vtx->n_vtx      = n_vtx;
  PDM_malloc(vtx->coords, 3 * n_vtx, double);
  vtx->_coords         = (double *) vtx->coords;
  vtx->_numabs         = (PDM_g_num_t *) numabs;
  vtx->_numparent      = (int *) num_parent;
  vtx->owner_coords    = ownership;
  vtx->owner_numabs    = ownership;
  vtx->owner_numparent = ownership;

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
 PDM_part_mesh_nodal_elmts_t *pmne
)
{
  if (pmne == NULL) {
    return;
  }

  assert(pmn->n_part == pmne->n_part);
  assert(pmn->mesh_dimension >= pmne->mesh_dimension);
  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_MAX;
  if(pmne->mesh_dimension == 3) {
    pmn->volumic          = pmne;
    geom_kind             = PDM_GEOMETRY_KIND_VOLUMIC;
  } else if(pmne->mesh_dimension == 2){
    pmn->surfacic          = pmne;
    geom_kind             = PDM_GEOMETRY_KIND_SURFACIC;
  } else if(pmne->mesh_dimension == 1){
    pmn->ridge          = pmne;
    geom_kind             = PDM_GEOMETRY_KIND_RIDGE;
  } else if(pmne->mesh_dimension == 0){
    pmn->corner          = pmne;
    geom_kind             = PDM_GEOMETRY_KIND_CORNER;
  } else {
    PDM_error (__FILE__, __LINE__, 0, "PDM_Mesh_nodal_add_dmesh_nodal_elmts bad mesh_dimension\n");
  }

  // update pmn->n_section, pmn->section_kind, pmn->section_id
  int n_section = PDM_part_mesh_nodal_elmts_n_section_get(pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  if (pmn->n_section + n_section >= pmn->s_section) {
    pmn->s_section = PDM_MAX(pmn->s_section, pmn->n_section + n_section);
    PDM_realloc(pmn->section_kind ,pmn->section_kind , pmn->s_section,PDM_geometry_kind_t);
    PDM_realloc(pmn->section_id   ,pmn->section_id   , pmn->s_section,int                );
  }


  for (int i = 0; i < n_section; i++) {
    int _id_section = pmn->n_section++;
    pmn->section_kind[_id_section] = geom_kind;
    pmn->section_id  [_id_section] = sections_id[i];
  }
}

void
PDM_part_mesh_nodal_part_comm_graph_set
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_part_comm_graph_t *pcg,
  PDM_geometry_kind_t    geom_kind,
  PDM_ownership_t        ownership
)
{
  pmn->pcg[geom_kind] = pcg;
  if (ownership!=PDM_OWNERSHIP_BAD_VALUE) {
    pmn->pcg_ownership[geom_kind] = ownership;
  }
}

void
PDM_part_mesh_nodal_part_comm_graph_get
(
  PDM_part_mesh_nodal_t  *pmn,
  PDM_geometry_kind_t     geom_kind,
  PDM_part_comm_graph_t **pcg,
  PDM_ownership_t         ownership
)
{
  *pcg = pmn->pcg[geom_kind];
  if (ownership!=PDM_OWNERSHIP_BAD_VALUE) {
    pmn->pcg_ownership[geom_kind] = ownership;
  }
}

void
PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind
)
{
  if (pmn->pcg[geom_kind]!=NULL) {
    return;
  }

  int          *n_entity    = NULL;
  PDM_g_num_t **entity_gnum = NULL;
  PDM_malloc(n_entity   , pmn->n_part, int);
  PDM_malloc(entity_gnum, pmn->n_part, PDM_g_num_t*);

  if (geom_kind==PDM_GEOMETRY_KIND_CORNER) {
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      n_entity   [i_part] = PDM_part_mesh_nodal_n_vtx_get    (pmn, i_part);
      entity_gnum[i_part] = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);
    }
  }
  else if (geom_kind==PDM_GEOMETRY_KIND_RIDGE || geom_kind==PDM_GEOMETRY_KIND_SURFACIC) {

    PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);

    int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
    int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

    for (int i_part = 0; i_part < pmn->n_part; i_part++) {

      int n_elt_tot = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);
      PDM_malloc(entity_gnum[i_part], n_elt_tot, PDM_g_num_t);
      n_entity[i_part] = 0;

      for (int i_section = 0; i_section < n_section; i_section++) {

        int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                                sections_id[i_section],
                                                                i_part);

        PDM_g_num_t *section_gnum = PDM_part_mesh_nodal_elmts_g_num_get(pmne, i_section, i_part, PDM_OWNERSHIP_BAD_VALUE);
        memcpy(&entity_gnum[i_part][n_entity[i_part]], section_gnum, n_elt*sizeof(PDM_g_num_t));
        n_entity[i_part] += n_elt;
      }
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum: invalid geom_kind (=%d)\n", geom_kind);
  }


  // Retrieve partition boundary entities using global IDs
  int n_rank;
  PDM_MPI_Comm_size(pmn->comm, &n_rank);

  PDM_g_num_t *part_distribution = PDM_compute_entity_distribution(pmn->comm, pmn->n_part);

  int  *n_entity_part_bound        = NULL;
  int **entity_proc_bound_idx      = NULL;
  int **entity_part_bound_idx      = NULL;
  int **entity_part_bound          = NULL;
  int **entity_part_bound_priority = NULL;
  PDM_part_generate_entity_graph_comm(pmn->comm,
                                      part_distribution,
                                      NULL,
                                      pmn->n_part,
                                      n_entity,
               (const PDM_g_num_t **) entity_gnum,
                                      NULL,
                                      &entity_proc_bound_idx,
                                      &entity_part_bound_idx,
                                      &entity_part_bound,
                                      &entity_part_bound_priority);
  PDM_malloc(n_entity_part_bound, pmn->n_part, int);
  for (int i_part = 0; i_part < pmn->n_part; i_part++) {
    n_entity_part_bound[i_part] = entity_part_bound_idx[i_part][part_distribution[n_rank]];

    PDM_log_trace_array_int(entity_part_bound[i_part], 4*n_entity_part_bound[i_part], "entity_part_bound[i_part] :: ");
    PDM_free(entity_proc_bound_idx     [i_part]);
    PDM_free(entity_part_bound_idx     [i_part]);
    PDM_free(entity_part_bound_priority[i_part]);
  }
  PDM_free(entity_proc_bound_idx);
  PDM_free(entity_part_bound_idx);
  PDM_free(entity_part_bound_priority);
  PDM_free(part_distribution);

  PDM_free(n_entity);
  if (geom_kind==PDM_GEOMETRY_KIND_RIDGE || geom_kind==PDM_GEOMETRY_KIND_SURFACIC) {
    for (int i_part = 0; i_part < pmn->n_part; i_part++) {
      PDM_free(entity_gnum[i_part]);
    }
  }
  PDM_free(entity_gnum);

  // Build part comm graph
  pmn->pcg[geom_kind] = PDM_part_comm_graph_create(pmn->n_part,
                                                   n_entity_part_bound,
                                                   entity_part_bound,
                                                   PDM_OWNERSHIP_KEEP,
                                                   pmn->comm);
  pmn->pcg_ownership[geom_kind] = PDM_OWNERSHIP_KEEP;

  PDM_free(n_entity_part_bound);
}

int
PDM_part_mesh_nodal_mesh_dimension_get
(
  PDM_part_mesh_nodal_t *pmn
)
{
  return pmn->mesh_dimension;
}

int
PDM_part_mesh_nodal_n_part_get
(
  PDM_part_mesh_nodal_t *pmn
)
{
  return pmn->n_part;
}

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

double*
PDM_part_mesh_nodal_vtx_coord_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part,
       PDM_ownership_t        ownership
)
{
  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    vtx->owner_coords = ownership;
  }

  return (double *) vtx->_coords;
}

PDM_g_num_t*
PDM_part_mesh_nodal_vtx_g_num_get
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    id_part,
       PDM_ownership_t        ownership
)
{
  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];

  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    vtx->owner_numabs = ownership;
  }
  return (PDM_g_num_t*) vtx->_numabs;
}

int
PDM_part_mesh_nodal_n_section_in_geom_kind_get
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


int *
PDM_part_mesh_nodal_sections_id_in_geom_kind_get
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
    return pmne->sections_id;
  } else {
    return NULL;
  }
}

PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_section_elt_type_get
(
        PDM_part_mesh_nodal_t *pmn,
  const int                    i_section
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  return PDM_part_mesh_nodal_section_in_geom_kind_elt_type_get(pmn,
                                                               geom_kind,
                                                               id_section);
}

PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_section_in_geom_kind_elt_type_get
(
        PDM_part_mesh_nodal_t *pmn,
        PDM_geometry_kind_t    geom_kind,
  const int                    id_section
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);
}


int
PDM_part_mesh_nodal_section_add
(
      PDM_part_mesh_nodal_t *pmn,
const PDM_Mesh_nodal_elt_t   t_elt
)
{
  PDM_geometry_kind_t geom_kind = PDM_Mesh_nodal_geom_kind_from_elt_type(t_elt);

  if( _get_from_geometry_kind(pmn, geom_kind) == NULL) {
    if(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC) {
      pmn->volumic = PDM_part_mesh_nodal_elmts_create(3,//pmn->mesh_dimension,
                                                      pmn->n_part,
                                                      pmn->comm);
    } else if( geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {
      pmn->surfacic = PDM_part_mesh_nodal_elmts_create(2,//pmn->mesh_dimension,
                                                       pmn->n_part,
                                                       pmn->comm);
    } else if( geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
      pmn->ridge = PDM_part_mesh_nodal_elmts_create(1,//pmn->mesh_dimension,
                                                    pmn->n_part,
                                                    pmn->comm);
    } else if( geom_kind == PDM_GEOMETRY_KIND_CORNER) {
      pmn->corner = PDM_part_mesh_nodal_elmts_create(0,//pmn->mesh_dimension,
                                                     pmn->n_part,
                                                     pmn->comm);
    }
  }


  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  int id_section = PDM_part_mesh_nodal_elmts_add(pmne, t_elt);

  if (pmn->n_section >= pmn->s_section) {
    pmn->s_section *= 2;
    PDM_realloc(pmn->section_kind ,pmn->section_kind , pmn->s_section,PDM_geometry_kind_t);
    PDM_realloc(pmn->section_id   ,pmn->section_id   , pmn->s_section,int                );
  }

  int _id_section = pmn->n_section++;
  pmn->section_kind[_id_section] = geom_kind;
  pmn->section_id  [_id_section] = id_section;

  return _id_section;
}

void
PDM_part_mesh_nodal_section_std_set
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    id_part,
const int                    n_elt,
const int                   *connec,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
const PDM_g_num_t           *parent_entity_g_num,
      PDM_ownership_t        owner
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_std_set(pmne, id_section, id_part, n_elt, connec, numabs, parent_num, parent_entity_g_num, owner);
}

void
PDM_part_mesh_nodal_section_std_ho_set
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
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
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

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

int
PDM_part_mesh_nodal_section_n_elt_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, id_part);
}

void
PDM_part_mesh_nodal_section_std_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      int                   **connec,
      PDM_g_num_t           **numabs,
      int                   **parent_num,
      PDM_g_num_t           **parent_entity_g_num,
      PDM_ownership_t         ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_section_std_get(pmne, id_section, id_part, connec, numabs, parent_num, parent_entity_g_num, ownership);
}

void
PDM_part_mesh_nodal_section_std_ho_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      int                   **connec,
      PDM_g_num_t           **numabs,
      int                   **parent_num,
      PDM_g_num_t           **parent_entity_g_num,
      int                    *order,
const char                  **ho_ordering,
      PDM_ownership_t         ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                               id_section,
                                               id_part,
                                               connec,
                                               numabs,
                                               parent_num,
                                               parent_entity_g_num,
                                               order,
                                               ho_ordering,
                                               ownership);
}


int *
PDM_part_mesh_nodal_section_parent_num_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      PDM_ownership_t         ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_parent_num_get(pmne, id_section, id_part, ownership);
}

PDM_g_num_t *
PDM_part_mesh_nodal_g_num_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      PDM_ownership_t         ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_g_num_get(pmne, id_section, id_part, ownership);
}

void
PDM_part_mesh_nodal_free
(
 PDM_part_mesh_nodal_t* pmn
)
{
  // volumic
  PDM_part_mesh_nodal_elmts_free(pmn->volumic);

  // surfacic
  PDM_part_mesh_nodal_elmts_free(pmn->surfacic);

  // ridge
  PDM_part_mesh_nodal_elmts_free(pmn->ridge);

  // corner
  PDM_part_mesh_nodal_elmts_free(pmn->corner);

  if (pmn->vtx != NULL) {
    for (int i_part = 0; i_part < pmn->n_part; i_part++) {
      _vtx_free (pmn->vtx[i_part]);
      PDM_free(pmn->vtx[i_part]);
    }

    PDM_free(pmn->vtx);
  }

  for (int geom_kind=0; geom_kind<PDM_GEOMETRY_KIND_MAX; ++geom_kind) {
    if (pmn->pcg_ownership[geom_kind]==PDM_OWNERSHIP_KEEP) {
      PDM_part_comm_graph_free(pmn->pcg[geom_kind]);
    }
  }

  PDM_free(pmn->section_kind);
  PDM_free(pmn->section_id);

  PDM_free(pmn);
}


void
PDM_part_mesh_nodal_dump_vtk
(
 PDM_part_mesh_nodal_t *pmn,
 PDM_geometry_kind_t    geom_kind,
 const char            *filename_pattern
)
{
  int i_rank = -1;
  PDM_MPI_Comm_rank(pmn->comm, &i_rank);

  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  if (pmne == NULL) {
    printf("Warning : PDM_part_mesh_nodal_dump_vtk : NULL pmne\n");
    return;
  }
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int pn_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    double      *pvtx_coord    = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);
    PDM_g_num_t *pvtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);

    int  n_section  = PDM_part_mesh_nodal_n_section_in_geom_kind_get  (pmn, geom_kind);
    int *section_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(pmn, geom_kind);

    /* Export group also */
    int n_group = pmne->n_group;
    int n_elt_tot = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);
    double *elt_group = PDM_array_const_double(n_elt_tot, -1);
    for(int i_group = 0; i_group < n_group; ++i_group) {
      for(int i = 0; i < pmne->n_group_elmt[i_part][i_group]; ++i) {
        int i_elt = pmne->group_elmt[i_part][i_group][i]-1;
        elt_group[i_elt] = i_group;
      }
    }

    int *elt_vtx_idx = NULL;
    int *elt_vtx     = NULL;
    PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne,
                                                   i_part,
                                                   &elt_vtx_idx,
                                                   &elt_vtx);

    PDM_g_num_t          *elt_g_num   = NULL;
    PDM_Mesh_nodal_elt_t *elt_type    = NULL;
    double               *elt_section = NULL;
    double               *elt_entity  = NULL;
    PDM_malloc(elt_g_num,   n_elt_tot, PDM_g_num_t         );
    PDM_malloc(elt_type,    n_elt_tot, PDM_Mesh_nodal_elt_t);
    PDM_malloc(elt_section, n_elt_tot, double              );
    PDM_malloc(elt_entity,  n_elt_tot, double              );

    int n_field = 3;

    int idx = 0;
    for (int i_section = 0; i_section < n_section; ++i_section) {
      int id_section = section_id[i_section];
      int                  n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get (pmne, id_section);

      PDM_g_num_t *g_num = PDM_part_mesh_nodal_elmts_g_num_get(pmne, id_section, i_part, PDM_OWNERSHIP_BAD_VALUE);

      int *_elt_to_entity = PDM_part_mesh_nodal_elmts_section_elmt_to_entity_get(pmne, id_section, i_part, PDM_OWNERSHIP_BAD_VALUE);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne, id_section, i_part, PDM_OWNERSHIP_BAD_VALUE);

      for (int i_elt = 0; i_elt < n_elt; i_elt++) {

        int i_parent = idx;
        if (parent_num != NULL) {
          i_parent = parent_num[i_elt];
        }

        if (g_num != NULL) {
          elt_g_num[i_parent] = g_num[i_elt];
        } else {
          elt_g_num[i_parent] = -1;
        }

        elt_type   [i_parent] = t_elt;
        elt_section[i_parent] = i_section;
        if (_elt_to_entity != NULL) {
          elt_entity[i_parent] = _elt_to_entity[i_elt];
        }
        else {
          n_field = 2;
        }
        idx++;
      }
    }

    const char   *field_name[] = {"groud_id", "section_id", "elt_to_entity"};
    const double *field_val [] = {elt_group, elt_section, elt_entity};

    char filename[999];
    sprintf(filename, "%s_%d_%d.vtk", filename_pattern, i_part, i_rank);
    PDM_vtk_write_unstructured_grid(filename,
                                    pn_vtx,
                                    pvtx_coord,
                                    pvtx_ln_to_gn,
                                    n_elt_tot,
                                    elt_type,
                                    elt_vtx_idx,
                                    elt_vtx,
                                    elt_g_num,
                                    n_field,
                                    field_name,
                                    field_val,
                                    0,
                                    NULL,
                                    NULL);

    PDM_free(elt_group  );
    PDM_free(elt_g_num  );
    PDM_free(elt_type   );
    PDM_free(elt_section);
    PDM_free(elt_entity );
    PDM_free(elt_vtx_idx);
    PDM_free(elt_vtx    );
  }
}

void
PDM_part_mesh_nodal_section_elt_extents_compute
(
       PDM_part_mesh_nodal_t *pmn,
 const int                    i_section,
 const int                    i_part,
 const double                 tolerance,
       double                *extents
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);

  PDM_part_mesh_nodal_elmts_elt_extents_compute(pmne,
                                                id_section,
                                                i_part,
                                                tolerance,
                                                vtx_coord,
                                                extents);
}


void
PDM_part_mesh_nodal_section_elt_center_compute
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    i_part,
const PDM_ownership_t        ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);

  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);

  PDM_part_mesh_nodal_elmts_elt_center_compute(pmne,
                                               id_section,
                                               i_part,
                                               n_vtx,
                                               vtx_coord,
                                               ownership);
}


const double *
PDM_part_mesh_nodal_section_elt_center_get
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    i_part,
      PDM_ownership_t        ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  return PDM_part_mesh_nodal_elmts_elt_center_get(pmne, id_section, i_part, ownership);
}

void
PDM_part_mesh_nodal_section_elt_center_reset
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    i_part
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  PDM_part_mesh_nodal_elmts_elt_center_reset(pmne, id_section, i_part);
}


void
PDM_part_mesh_nodal_section_poly2d_set
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    id_part,
const int                    n_elt,
const int                   *connec_idx,
const int                   *connec,
const PDM_g_num_t           *numabs,
const int                   *parent_num,
      PDM_ownership_t        owner
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);
  assert(geom_kind == PDM_GEOMETRY_KIND_SURFACIC);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_section_poly2d_set(pmne,
                                               id_section,
                                               id_part,
                                               n_elt,
                                               connec_idx,
                                               connec,
                                               numabs,
                                               parent_num,
                                               owner);
}


void
PDM_part_mesh_nodal_section_poly2d_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      int                   **connec_idx,
      int                   **connec,
      PDM_ownership_t         ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);
assert(geom_kind == PDM_GEOMETRY_KIND_SURFACIC);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                               id_section,
                                               id_part,
                                               connec_idx,
                                               connec,
                                               ownership);
}

void
PDM_part_mesh_nodal_section_poly3d_set
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
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
const PDM_g_num_t           *parent_entity_g_num,
      PDM_ownership_t        owner
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);
  assert(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_section_poly3d_set(pmne,
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
                                               parent_entity_g_num,
                                               owner);
}

void
PDM_part_mesh_nodal_section_poly3d_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      int                    *n_face,
      PDM_g_num_t           **face_ln_to_gn,
      int                   **face_vtx_idx,
      int                   **face_vtx,
      PDM_g_num_t           **numabs,
      int                   **cell_face_idx,
      int                   **cell_face,
      int                   **parent_num,
      PDM_g_num_t           **parent_entity_g_num,
      PDM_ownership_t         ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);
  assert(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_section_poly3d_get(pmne,
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
                                               parent_entity_g_num,
                                               ownership);
}

void
PDM_part_mesh_nodal_section_poly3d_cell_vtx_connect_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      int                   **cellvtx_idx,
      int                   **cellvtx,
      PDM_ownership_t         ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);
  assert(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_section_poly3d_cell_vtx_connect_get(pmne,
                                                                id_section,
                                                                id_part,
                                                                cellvtx_idx,
                                                                cellvtx,
                                                                ownership);
}

void
PDM_part_mesh_nodal_reset
(
 PDM_part_mesh_nodal_t *pmn
)
{
  for (PDM_geometry_kind_t geom_kind = (PDM_geometry_kind_t) 0; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {
    PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);

    if (pmne != NULL) {
      PDM_part_mesh_nodal_elmts_reset(pmne);
    }

  }

  pmn->n_section = 0;

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
        PDM_free(pmn->vtx[i]->coords);
        pmn->vtx[i]->coords = NULL;
      }
    }
  }
}

void
PDM_part_mesh_nodal_g_num_in_section_compute
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const PDM_ownership_t         ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  PDM_part_mesh_nodal_elmts_g_num_in_section_compute(pmne,
                                                     id_section,
                                                     ownership);
}

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


PDM_g_num_t *
PDM_part_mesh_nodal_g_num_get_from_part
(
      PDM_part_mesh_nodal_t  *pmn,
      PDM_geometry_kind_t     geom_kind,
const int                     id_part,
      PDM_ownership_t         ownership
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  return PDM_part_mesh_nodal_elmts_g_num_get_from_part(pmne, id_part, ownership);
}

void
PDM_part_mesh_nodal_partial_free
(
 PDM_part_mesh_nodal_t *pmn
)
{
  for (PDM_geometry_kind_t geom_kind = (PDM_geometry_kind_t)  0; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {
    PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);

    if (pmne != NULL) {
      PDM_part_mesh_nodal_elmts_partial_free(pmne);
    }

  }
}

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

PDM_g_num_t *
PDM_part_mesh_nodal_section_g_num_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     i_section,
const int                     id_part,
      PDM_ownership_t         ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                  i_section,
                                                  &geom_kind,
                                                  &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  return PDM_part_mesh_nodal_elmts_section_g_num_get(pmne, id_section, id_part, ownership);
}

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

int *
PDM_part_mesh_nodal_section_elmt_to_entity_get
(
      PDM_part_mesh_nodal_t *pmn,
const int                    i_section,
const int                    id_part,
      PDM_ownership_t        ownership
)
{
  PDM_geometry_kind_t geom_kind;
  int                 id_section;
  PDM_part_mesh_nodal_section_id_and_geom_kind_get(pmn,
                                                   i_section,
                                                   &geom_kind,
                                                   &id_section);

  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);

  return PDM_part_mesh_nodal_elmts_section_elmt_to_entity_get(pmne,
                                                              id_section,
                                                              id_part,
                                                              ownership);
}

void
PDM_part_mesh_nodal_group_get
(
       PDM_part_mesh_nodal_t  *pmn,
       PDM_geometry_kind_t     geom_kind,
 const int                     i_part,
 const int                     i_group,
       int                    *n_group_elmt,
       int                   **group_elmt,
       PDM_g_num_t           **group_ln_to_gn,
       PDM_ownership_t         ownership
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_group_get(pmne,
                                      i_part,
                                      i_group,
                                      n_group_elmt,
                                      group_elmt,
                                      group_ln_to_gn,
                                      ownership);
}


int*
PDM_part_mesh_nodal_compute_sections_idx
(
 PDM_part_mesh_nodal_t  *pmn,
 PDM_geometry_kind_t     geom_kind,
 const int               id_part
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_compute_sections_idx(pmne, id_part);
}


const int *
PDM_part_mesh_nodal_vertices_parent_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part
 )
{
  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];

  return vtx->_numparent;
}


const PDM_g_num_t *
PDM_part_mesh_nodal_vertices_g_num_parent_get
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part
 )
{
  if (pmn == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_part >= pmn->n_part) {
    PDM_error (__FILE__, __LINE__, 0, "Bad part identifier\n");
  }

  PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[id_part];

  assert(vtx->parent != NULL);

  return vtx->parent->_numabs;
}


void
PDM_part_mesh_nodal_cell3d_cellface_add
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part,
const int                     n_cell,
const int                     n_face,
const int                    *face_vtx_idx,
const int                    *face_vtx,
const PDM_g_num_t            *face_ln_to_gn,
const int                    *cell_face_idx,
const int                    *cell_face,
const PDM_g_num_t            *cell_ln_to_gn,
const PDM_ownership_t         ownership
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, PDM_GEOMETRY_KIND_VOLUMIC);

  if (pmne == NULL) {
    pmne = PDM_part_mesh_nodal_elmts_create(3, pmn->n_part, pmn->comm);
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmne);
  }

  int n_section_before = PDM_part_mesh_nodal_n_section_in_geom_kind_get(pmn,
                                                                        PDM_GEOMETRY_KIND_VOLUMIC);

  PDM_part_mesh_elmts_nodal_cell3d_cellface_add(pmne,
                                                id_part,
                                                n_cell,
                                                n_face,
                                                face_vtx_idx,
                                                face_vtx,
                                                face_ln_to_gn,
                                                cell_face_idx,
                                                cell_face,
                                                cell_ln_to_gn,
                                                pmn->vtx,
                                                ownership);

  // update pmn->n_section, pmn->section_kind, pmn->section_id
  int n_section_after = PDM_part_mesh_nodal_n_section_in_geom_kind_get(pmn,
                                                                       PDM_GEOMETRY_KIND_VOLUMIC);
  int *sections_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(pmn,
                                                                      PDM_GEOMETRY_KIND_VOLUMIC);
  if (pmn->n_section + n_section_after - n_section_before >= pmn->s_section) {
    pmn->s_section = PDM_MAX(pmn->s_section, pmn->n_section + n_section_after - n_section_before);
    PDM_realloc(pmn->section_kind ,pmn->section_kind , pmn->s_section,PDM_geometry_kind_t);
    PDM_realloc(pmn->section_id   ,pmn->section_id   , pmn->s_section,int                );
  }

  for (int i = n_section_before; i < n_section_after; i++) {
    int _id_section = pmn->n_section++;
    pmn->section_kind[_id_section] = PDM_GEOMETRY_KIND_VOLUMIC;
    pmn->section_id  [_id_section] = sections_id[i];
  }
}

void
PDM_part_mesh_nodal_face2d_faceedge_add
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part,
const int                     n_face,
const int                     n_edge,
const int                    *edge_vtx,
const int                    *face_edge_idx,
const int                    *face_edge,
const PDM_g_num_t            *face_ln_to_gn,
const PDM_ownership_t         ownership
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, PDM_GEOMETRY_KIND_SURFACIC);

  if (pmne == NULL) {
    pmne = PDM_part_mesh_nodal_elmts_create(2, pmn->n_part, pmn->comm);
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmne);
  }

  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, id_part);

  int n_section_before = PDM_part_mesh_nodal_n_section_in_geom_kind_get(pmn,
                                                                        PDM_GEOMETRY_KIND_SURFACIC);

  PDM_part_mesh_nodal_elmts_face2d_faceedge_add(pmne,
                                                id_part,
                                                n_face,
                                                n_edge,
                                                edge_vtx,
                                                face_edge_idx,
                                                face_edge,
                                                face_ln_to_gn,
                                                n_vtx,
                                                ownership);

  // update pmn->n_section, pmn->section_kind, pmn->section_id
  int n_section_after = PDM_part_mesh_nodal_n_section_in_geom_kind_get(pmn,
                                                                       PDM_GEOMETRY_KIND_SURFACIC);
  int *sections_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(pmn,
                                                                      PDM_GEOMETRY_KIND_SURFACIC);

  if (pmn->n_section + n_section_after - n_section_before >= pmn->s_section) {
    pmn->s_section = PDM_MAX(pmn->s_section, pmn->n_section + n_section_after - n_section_before);
    PDM_realloc(pmn->section_kind ,pmn->section_kind , pmn->s_section,PDM_geometry_kind_t);
    PDM_realloc(pmn->section_id   ,pmn->section_id   , pmn->s_section,int                );
  }

  for (int i = n_section_before; i < n_section_after; i++) {
    int _id_section = pmn->n_section++;
    pmn->section_kind[_id_section] = PDM_GEOMETRY_KIND_SURFACIC;
    pmn->section_id  [_id_section] = sections_id[i];
  }
}

void
PDM_part_mesh_nodal_cells_cellvtx_add
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part,
const int                     n_cell,
const int                    *cell_vtx_idx,
const int                    *cell_vtx,
const PDM_g_num_t            *numabs,
const PDM_ownership_t         ownership
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, PDM_GEOMETRY_KIND_VOLUMIC);

  if (pmne == NULL) {
    pmne = PDM_part_mesh_nodal_elmts_create(3, pmn->n_part, pmn->comm);
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmne);
  }

  int n_section_before = PDM_part_mesh_nodal_n_section_in_geom_kind_get(pmn,
                                                                        PDM_GEOMETRY_KIND_VOLUMIC);

  PDM_part_mesh_nodal_elmts_cells_cellvtx_add(pmne,
                                              id_part,
                                              n_cell,
                                              cell_vtx_idx,
                                              cell_vtx,
                                              numabs,
                                              ownership);

  // update pmn->n_section, pmn->section_kind, pmn->section_id
  int n_section_after = PDM_part_mesh_nodal_n_section_in_geom_kind_get(pmn,
                                                                       PDM_GEOMETRY_KIND_VOLUMIC);
  int *sections_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(pmn,
                                                                      PDM_GEOMETRY_KIND_VOLUMIC);

  if (pmn->n_section + n_section_after - n_section_before >= pmn->s_section) {
    pmn->s_section = PDM_MAX(pmn->s_section, pmn->n_section + n_section_after - n_section_before);
    PDM_realloc(pmn->section_kind ,pmn->section_kind , pmn->s_section,PDM_geometry_kind_t);
    PDM_realloc(pmn->section_id   ,pmn->section_id   , pmn->s_section,int                );
  }

  for (int i = n_section_before; i < n_section_after; i++) {
    int _id_section = pmn->n_section++;
    pmn->section_kind[_id_section] = PDM_GEOMETRY_KIND_VOLUMIC;
    pmn->section_id  [_id_section] = sections_id[i];
  }
}


void
PDM_part_mesh_nodal_faces_facevtx_add
(
      PDM_part_mesh_nodal_t  *pmn,
const int                     id_part,
const int                     n_face,
const int                    *face_vtx_idx,
const int                    *face_vtx,
const PDM_g_num_t            *numabs,
const PDM_ownership_t         ownership
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, PDM_GEOMETRY_KIND_SURFACIC);

  if (pmne == NULL) {
    pmne = PDM_part_mesh_nodal_elmts_create(2, pmn->n_part, pmn->comm);
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmne);
  }

  int n_section_before = PDM_part_mesh_nodal_n_section_in_geom_kind_get(pmn,
                                                                        PDM_GEOMETRY_KIND_SURFACIC);

  PDM_part_mesh_nodal_elmts_faces_facevtx_add(pmne,
                                              id_part,
                                              n_face,
                                              face_vtx_idx,
                                              face_vtx,
                                              numabs,
                                              ownership);

  // update pmn->n_section, pmn->section_kind, pmn->section_id
  int n_section_after = PDM_part_mesh_nodal_n_section_in_geom_kind_get(pmn,
                                                                       PDM_GEOMETRY_KIND_SURFACIC);
  int *sections_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(pmn,
                                                                      PDM_GEOMETRY_KIND_SURFACIC);

  if (pmn->n_section + n_section_after - n_section_before >= pmn->s_section) {
    pmn->s_section = PDM_MAX(pmn->s_section, pmn->n_section + n_section_after - n_section_before);
    PDM_realloc(pmn->section_kind ,pmn->section_kind , pmn->s_section,PDM_geometry_kind_t);
    PDM_realloc(pmn->section_id   ,pmn->section_id   , pmn->s_section,int                );
  }

  for (int i = n_section_before; i < n_section_after; i++) {
    int _id_section = pmn->n_section++;
    pmn->section_kind[_id_section] = PDM_GEOMETRY_KIND_SURFACIC;
    pmn->section_id  [_id_section] = sections_id[i];
  }
}

void
PDM_part_mesh_nodal_section_id_and_geom_kind_get
(
       PDM_part_mesh_nodal_t  *pmn,
 const int                     i_section,
       PDM_geometry_kind_t    *geom_kind,
       int                    *id_section_in_geom_kind
 )
{
  if (i_section >= pmn->n_section) {
    PDM_error(__FILE__, __LINE__, 0, "i_section (%d) > n_section (%d)\n", i_section, pmn->n_section);
  }

  *geom_kind               = pmn->section_kind[i_section];
  *id_section_in_geom_kind = pmn->section_id  [i_section];
}

int
PDM_part_mesh_nodal_section_id_from_geom_kind_get
(
       PDM_part_mesh_nodal_t  *pmn,
 const PDM_geometry_kind_t     geom_kind,
 const int                     id_section_in_geom_kind
 )
{
  int i_section = 0;

  for (i_section = 0; i_section < pmn->n_section; i_section++) {
    if (pmn->section_id  [i_section] == id_section_in_geom_kind &&
        pmn->section_kind[i_section] == geom_kind) {
      return i_section;
    }
  }

  return -1;
}

int
PDM_part_mesh_nodal_n_section_get
(
 PDM_part_mesh_nodal_t *pmn
)
{
  assert(pmn != NULL);

  return pmn->n_section;
}

int *
PDM_part_mesh_nodal_sections_id_get
(
 PDM_part_mesh_nodal_t *pmn
)
{
  assert(pmn != NULL);

  return pmn->section_id;
}


void
PDM_part_mesh_nodal_n_group_set
(
       PDM_part_mesh_nodal_t  *pmn,
       PDM_geometry_kind_t     geom_kind,
 const int                     n_group
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_n_group_set(pmne, n_group);
}

void
PDM_part_mesh_nodal_group_set
(
       PDM_part_mesh_nodal_t  *pmn,
       PDM_geometry_kind_t     geom_kind,
 const int                     i_part,
 const int                     i_group,
       int                     n_group_elmt,
       int                    *group_elmt,
       PDM_g_num_t            *group_ln_to_gn,
       PDM_ownership_t         ownership
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  PDM_part_mesh_nodal_elmts_group_set(pmne,
                                      i_part,
                                      i_group,
                                      n_group_elmt,
                                      group_elmt,
                                      group_ln_to_gn,
                                      ownership);
}

int
PDM_part_mesh_nodal_n_group_get
(
       PDM_part_mesh_nodal_t  *pmn,
       PDM_geometry_kind_t     geom_kind
)
{
  PDM_part_mesh_nodal_elmts_t* pmne = _get_from_geometry_kind(pmn, geom_kind);
  assert(pmne != NULL);
  return PDM_part_mesh_nodal_elmts_n_group_get(pmne);
}

PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_part_mesh_nodal_elmts_get
(
 PDM_part_mesh_nodal_t  *pmn,
 PDM_geometry_kind_t     geom_kind
)
{
  return _get_from_geometry_kind(pmn, geom_kind);
}

PDM_geometry_kind_t
PDM_part_mesh_nodal_principal_geom_kind_get
(
 PDM_part_mesh_nodal_t  *pmn
 )
{
  switch (pmn->mesh_dimension) {
  case 3:
    return PDM_GEOMETRY_KIND_VOLUMIC;
    break;
  case 2:
    return PDM_GEOMETRY_KIND_SURFACIC;
    break;
  case 1:
    return PDM_GEOMETRY_KIND_RIDGE;
    break;
  case 0:
    return PDM_GEOMETRY_KIND_CORNER;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid mesh_dimension %d\n", pmn->mesh_dimension);
  }

  return PDM_GEOMETRY_KIND_MAX;
}

int
PDM_part_mesh_nodal_cell_vtx_connect_get
(
        PDM_part_mesh_nodal_t  *pmn,
        PDM_geometry_kind_t     geom_kind,
  const int                     i_part,
        int                   **cell_vtx_idx,
        int                   **cell_vtx
)
{
  if (pmn == NULL) {
    return 0;
  }

  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);
  if (i_part >= n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Invalid i_part (%d / %d)\n", i_part, n_part);
  }

  PDM_part_mesh_nodal_elmts_t *pmne = _get_from_geometry_kind(pmn,
                                                              geom_kind);

  return PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne,
                                                        i_part,
                                                        cell_vtx_idx,
                                                        cell_vtx);
}


static void
_transpose_group_information
(
  PDM_MPI_Comm    comm,
  int             n_part,
  int            *n_entity,
  int           **entity_group_idx,
  PDM_g_num_t   **entity_group,
  int            *out_n_group,
  int          ***out_group_entity_idx,
  int          ***out_group_entity
)
{
  int debug_verbose = 1;

  int **_entity_group     = NULL;
  int **_group_entity_idx = NULL;
  int **_group_entity     = NULL;
  PDM_malloc(_entity_group    , n_part, int *);
  PDM_malloc(_group_entity_idx, n_part, int *);
  PDM_malloc(_group_entity    , n_part, int *);

  int ln_group = 0;
  int gn_group = 0;
  for (int i_part=0; i_part<n_part; ++i_part) {

    PDM_malloc(_entity_group[i_part], entity_group_idx[i_part][n_entity[i_part]], int);

    for (int i_group=0; i_group<entity_group_idx[i_part][n_entity[i_part]]; ++i_group) {
      ln_group = PDM_MAX(entity_group[i_part][i_group], ln_group);
      _entity_group[i_part][i_group] = (int) entity_group[i_part][i_group];
    }

    if (debug_verbose==1) {
      int entity_group_size = entity_group_idx[i_part][n_entity[i_part]];
      PDM_log_trace_array_int(entity_group_idx[i_part], n_entity[i_part]+1       , "entity_group_idx :: ");
      PDM_log_trace_array_int(_entity_group   [i_part], entity_group_size, "_entity_group     :: ");
    }
  }

  PDM_MPI_Allreduce(&ln_group, &gn_group, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  for (int i_part=0; i_part<n_part; ++i_part) {
    PDM_connectivity_transpose(n_entity[i_part],
                               gn_group,
                               entity_group_idx [i_part],
                               _entity_group    [i_part],
                              &_group_entity_idx[i_part],
                              &_group_entity    [i_part]);

    if (debug_verbose==1) {
      int group_entity_size = _group_entity_idx[i_part][gn_group];
      PDM_log_trace_array_int(_group_entity_idx[i_part], gn_group+1   , "group_entity_idx :: ");
      PDM_log_trace_array_int(_group_entity    [i_part], group_entity_size, "group_entity     :: ");
    }

    PDM_free(_entity_group[i_part]);
  }
  PDM_free(_entity_group);

  *out_n_group = gn_group;
  *out_group_entity_idx = _group_entity_idx;
  *out_group_entity     = _group_entity;

}


static void
_generate_group_entity_gid
(
  PDM_MPI_Comm    comm,
  int             n_part,
  PDM_g_num_t   **entity_gid,
  int             n_group,
  int           **group_entity_idx,
  int           **group_entity,
  PDM_g_num_t  ***out_group_entity_gid
)
{
  PDM_g_num_t **group_entity_gid = NULL;
  PDM_malloc(group_entity_gid, n_part, PDM_g_num_t *);

  for (int i_part=0; i_part<n_part; ++i_part) {

    PDM_malloc(group_entity_gid[i_part], group_entity_idx[i_part][n_group], PDM_g_num_t);

    for (int i_entity=0; i_entity<group_entity_idx[i_part][n_group]; ++i_entity) {
      group_entity_gid[i_part][i_entity] = entity_gid[i_part][group_entity[i_part][i_entity]-1];
    }
  }

  for (int i_group=0; i_group<n_group; ++i_group) {

    PDM_gen_gnum_t* gen_group_gid = PDM_gnum_create(3,
                                                     n_part,
                                                     PDM_TRUE,
                                                     1e-4,
                                                     comm,
                                                     PDM_OWNERSHIP_USER);
    PDM_gnum_set_parents_nuplet(gen_group_gid, 1);
    for (int i_part=0; i_part<n_part; ++i_part) {
      int n_group_entity = group_entity_idx[i_part][i_group+1]-group_entity_idx[i_part][i_group];
      PDM_gnum_set_from_parents(gen_group_gid, i_part, n_group_entity, &group_entity_gid[i_part][group_entity_idx[i_part][i_group]]);
    }

    PDM_gnum_compute(gen_group_gid);

    for (int i_part=0; i_part<n_part; ++i_part) {
      int n_group_entity = group_entity_idx[i_part][i_group+1]-group_entity_idx[i_part][i_group];
      PDM_g_num_t *group_gid = PDM_gnum_get(gen_group_gid, i_part);
      for (int i_entity=0; i_entity<n_group_entity; ++i_entity) {
        group_entity_gid[i_part][group_entity_idx[i_part][i_group]+i_entity] = group_gid[i_entity];
      }
      PDM_free(group_gid);
    }

    PDM_gnum_free(gen_group_gid);
  }

  *out_group_entity_gid = group_entity_gid;
}


static void
_generate_group_gid
(
  PDM_MPI_Comm   comm,
  int            n_part,
  int           *n_entity,
  int          **entity_tag_idx,
  int          **entity_tag,
  int         ***out_entity_group_idx,
  PDM_g_num_t ***out_entity_group
)
{

  int           n_parent_max            = 0;
  int         **group_parent_n          = NULL;
  int         **entity_group_idx        = NULL;
  // int         **entity_parent_group_idx = NULL;
  // int         **entity_parent_group     = NULL;
  PDM_g_num_t **group_parent_nplt       = NULL;
  PDM_malloc(group_parent_n         , n_part, int         *);
  PDM_malloc(entity_group_idx       , n_part, int         *);
  // PDM_malloc(entity_parent_group_idx, n_part, int         *);
  // PDM_malloc(entity_parent_group    , n_part, int         *);
  PDM_malloc(group_parent_nplt      , n_part, PDM_g_num_t *);

  for (int i_part=0; i_part<n_part; ++i_part) {
    /**
     * All parent group are known by entities and coherent over all procs, unique parents and : 
     *  - count n_parent_max for gnum gnum
     *  - count n_parent for each entity
     *  - count straddling entity in entity_group_idx
     */
    // int entity_parent_group_size = 0;

    PDM_malloc(group_parent_n  [i_part], n_entity[i_part]  , int);
    PDM_malloc(entity_group_idx[i_part], n_entity[i_part]+1, int); entity_group_idx[i_part][0] = 0;
    for(int i_entity = 0; i_entity < n_entity[i_part]; ++i_entity) {
      group_parent_n[i_part][i_entity  ] = 0;
      entity_group_idx [i_part][i_entity+1] = entity_group_idx[i_part][i_entity];
      int beg = entity_tag_idx[i_part][i_entity  ];
      int end = entity_tag_idx[i_part][i_entity+1];
      if(end - beg > 0){
        int n_unique = PDM_inplace_unique(entity_tag[i_part], beg, end-1);
        n_parent_max = PDM_MAX(n_parent_max, n_unique);
        if(n_unique > 1) {
          // entity_parent_group_size  += n_unique;
          group_parent_n  [i_part][i_entity  ] = n_unique;
          entity_group_idx[i_part][i_entity+1] = entity_group_idx[i_part][i_entity]+1;
        }
      }
    }


    /**
     * Prepare nuplet from group gid computation
     *
     * Here we could store parent group at same time
     */
    // int i_write_parent = 0;
    int i_write_nuplet = 0;
    // PDM_malloc(entity_parent_group_idx[i_part], n_entity[i_part]+1                                     , int        ); entity_parent_group_idx[i_part][0] = 0;
    // PDM_malloc(entity_parent_group    [i_part], entity_parent_group_size                               , int        );
    PDM_malloc(group_parent_nplt      [i_part], n_parent_max*entity_group_idx[i_part][n_entity[i_part]], PDM_g_num_t);
    for(int i_entity = 0; i_entity < n_entity[i_part]; ++i_entity) {
      // entity_parent_group_idx[i_part][i_entity+1] = entity_parent_group_idx[i_part][i_entity];
      int beg = entity_tag_idx[i_part][i_entity  ];
      if (group_parent_n[i_part][i_entity]>0) {
        // entity_parent_group_idx[i_part][i_entity+1] += group_parent_n[i_part][i_entity];
        for (int i_read=beg; i_read<beg+group_parent_n[i_part][i_entity]; ++i_read) {
          group_parent_nplt[i_part][i_write_nuplet++] = (PDM_g_num_t) entity_tag[i_part][i_read];
          // entity_parent_group [i_part][i_write_parent++] = entity_tag[i_part][i_read];
        }
        for (int i_read=beg+group_parent_n[i_part][i_entity]; i_read<beg+n_parent_max; ++i_read) {
          group_parent_nplt[i_part][i_write_nuplet++] = 0;
        }
      }
    }
    PDM_free(group_parent_n[i_part]);
  }
  PDM_free(group_parent_n);


  /**
   * Generate ids for groups
   */
  PDM_gen_gnum_t* gen_group_id = PDM_gnum_create(3,
                                                 n_part,
                                                 PDM_TRUE,
                                                 1e-4,
                                                 comm,
                                                 PDM_OWNERSHIP_USER);
  PDM_gnum_set_parents_nuplet(gen_group_id, n_parent_max);
  for (int i_part=0; i_part<n_part; ++i_part) {
    PDM_gnum_set_from_parents(gen_group_id, i_part, entity_group_idx[i_part][n_entity[i_part]], group_parent_nplt[i_part]);
  }
  PDM_gnum_compute(gen_group_id);

  PDM_g_num_t **group_id  = NULL;
  // int         **entity_group = NULL;
  PDM_malloc(group_id , n_part, PDM_g_num_t *);
  // PDM_malloc(entity_group, n_part, int         *);
  for (int i_part=0; i_part<n_part; ++i_part) {
    group_id[i_part] = PDM_gnum_get(gen_group_id, i_part);
    // PDM_malloc(entity_group[i_part], entity_group_idx[i_part][n_entity[i_part]], int);
    PDM_log_trace_connectivity_long(entity_group_idx[i_part], group_id[i_part], n_entity[i_part], "group_gid :: ");
    PDM_free(group_parent_nplt[i_part]);

  }

  PDM_gnum_free(gen_group_id);
  PDM_free(group_parent_nplt);

  *out_entity_group_idx = entity_group_idx;
  *out_entity_group     = group_id;
}


void
PDM_part_mesh_nodal_compute_straddling_entities
(
  PDM_part_mesh_nodal_t  *pmn,
  PDM_geometry_kind_t     geom_kind,
  PDM_geometry_kind_t     geom_kind_tgt
)
{
  int debug_verbose = 1;
  int debug_visu    = 1;

  int i_rank = -1;
  PDM_MPI_Comm_rank(pmn->comm, &i_rank);


  if (!((geom_kind==PDM_GEOMETRY_KIND_SURFACIC && (geom_kind_tgt==PDM_GEOMETRY_KIND_RIDGE ||
                                                   geom_kind_tgt==PDM_GEOMETRY_KIND_CORNER  ) ) ||
        (geom_kind==PDM_GEOMETRY_KIND_RIDGE    &&  geom_kind_tgt==PDM_GEOMETRY_KIND_CORNER    )  )) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_compute_topo_corners cannot build entities of geom_kind %d from entities of geom_kind %d\n", geom_kind_tgt, geom_kind);
  }


  PDM_part_mesh_nodal_elmts_t* pmne = NULL;
  if (geom_kind==PDM_GEOMETRY_KIND_SURFACIC) {
    pmne = pmn->surfacic;
  }
  else if (geom_kind==PDM_GEOMETRY_KIND_RIDGE) {
    pmne = pmn->ridge;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_compute_topo_corners not implemented for geom_kind %d\n", geom_kind);
  }

  if (pmne==NULL) {
    if (i_rank==0) {
      printf("Warning: no entity of geom_kind %d found\n", geom_kind);
    }
    return;
  }

  /**
   * Detect locally which vertices are on multiple entity1 groups
   */
  int  *n_vtx           = NULL;
  int  *n_entity1       = NULL;
  int **entity1_vtx_idx = NULL;
  int **entity1_vtx     = NULL;
  int **tag_vtx_n       = NULL;
  int **tag_vtx_idx     = NULL;
  int **tag_vtx         = NULL;
  PDM_malloc(n_vtx          , pmn->n_part, int  );
  PDM_malloc(n_entity1      , pmn->n_part, int  );
  PDM_malloc(entity1_vtx_idx, pmn->n_part, int *);
  PDM_malloc(entity1_vtx    , pmn->n_part, int *);
  PDM_malloc(tag_vtx_n      , pmn->n_part, int *);
  PDM_malloc(tag_vtx_idx    , pmn->n_part, int *);
  PDM_malloc(tag_vtx        , pmn->n_part, int *);

  int n_group_entity1 = PDM_part_mesh_nodal_elmts_n_group_get(pmne);

  for (int i_part=0; i_part<pmn->n_part; ++i_part) {
    n_vtx    [i_part] = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    n_entity1[i_part] = PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne, i_part,
                                                                      &entity1_vtx_idx[i_part],
                                                                      &entity1_vtx[i_part]);

    tag_vtx_n[i_part]  = PDM_array_const_int(n_vtx[i_part], 0);

    for(int i_group = 0; i_group < n_group_entity1; ++i_group) {
      int          n_entity1_group        = 0;
      int         *entity1_group          = NULL;
      PDM_g_num_t *entity1_group_ln_to_gn = NULL;
      PDM_part_mesh_nodal_elmts_group_get(pmne,
                                          i_part,
                                          i_group,
                                         &n_entity1_group,
                                         &entity1_group,
                                         &entity1_group_ln_to_gn,
                                          PDM_OWNERSHIP_BAD_VALUE);

      for(int idx_elmt = 0; idx_elmt < n_entity1_group; ++idx_elmt) {
        int i_elmt = entity1_group[idx_elmt]-1;
        for (int i_elmt_vtx=entity1_vtx_idx[i_part][i_elmt];
                 i_elmt_vtx<entity1_vtx_idx[i_part][i_elmt+1]; ++i_elmt_vtx) {
          int i_vtx = entity1_vtx[i_part][i_elmt_vtx]-1;
          tag_vtx_n[i_part][i_vtx]++;
        }
      }
    }

    PDM_malloc(tag_vtx_idx[i_part], n_vtx[i_part]+1, int); tag_vtx_idx[i_part][0] = 0;
    for(int i_vtx = 0; i_vtx < n_vtx[i_part]; ++i_vtx) {
      tag_vtx_idx[i_part][i_vtx+1] = tag_vtx_idx[i_part][i_vtx] + tag_vtx_n[i_part][i_vtx];
      tag_vtx_n  [i_part][i_vtx  ] = 0;
    }

    PDM_malloc(tag_vtx[i_part], tag_vtx_idx[i_part][n_vtx[i_part]], int);
    for(int i_group = 0; i_group < n_group_entity1; ++i_group) {
      int          n_entity1_group        = 0;
      int         *entity1_group          = NULL;
      PDM_g_num_t *entity1_group_ln_to_gn = NULL;
      PDM_part_mesh_nodal_elmts_group_get(pmne,
                                          i_part,
                                          i_group,
                                         &n_entity1_group,
                                         &entity1_group,
                                         &entity1_group_ln_to_gn,
                                          PDM_OWNERSHIP_BAD_VALUE);

      for(int idx_elmt = 0; idx_elmt < n_entity1_group; ++idx_elmt) {
        int i_elmt = entity1_group[idx_elmt]-1;
        for (int i_elmt_vtx=entity1_vtx_idx[i_part][i_elmt];
                 i_elmt_vtx<entity1_vtx_idx[i_part][i_elmt+1]; ++i_elmt_vtx) {
          int i_vtx = entity1_vtx[i_part][i_elmt_vtx]-1;
          int idx_write = tag_vtx_idx[i_part][i_vtx] + tag_vtx_n[i_part][i_vtx]++;
          tag_vtx[i_part][idx_write] = i_group+1;
        }
      }
    }

    PDM_free(entity1_vtx_idx[i_part]);
    PDM_free(entity1_vtx    [i_part]);
  }
  PDM_free(entity1_vtx_idx);
  PDM_free(entity1_vtx);


  /**
   * Exchange local information to reduce it globally
   */
  if (pmn->pcg[PDM_GEOMETRY_KIND_CORNER]==NULL) {
    PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum(pmn, PDM_GEOMETRY_KIND_CORNER);
  }


  if (debug_verbose==1) {
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_log_trace_connectivity_int(tag_vtx_idx[i_part], tag_vtx[i_part], n_vtx[i_part], "local tag_vtx :: ");
    }
  }
  
  int **_tag_vtx_idx = NULL;
  int **_tag_vtx     = NULL;
  PDM_part_comm_graph_gather_strided_int_data(pmn->pcg[PDM_GEOMETRY_KIND_CORNER],
                                              n_vtx,
                                              tag_vtx_idx,
                                              tag_vtx,
                                            &_tag_vtx_idx,
                                            &_tag_vtx);
  for (int i_part=0; i_part<pmn->n_part; ++i_part) {
    PDM_free(tag_vtx_idx[i_part]);
    PDM_free(tag_vtx    [i_part]);
  }
  PDM_free(tag_vtx_idx);
  PDM_free(tag_vtx);
  tag_vtx_idx = _tag_vtx_idx;
  tag_vtx     = _tag_vtx;


  if (debug_verbose==1) {
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_log_trace_connectivity_int(tag_vtx_idx[i_part], tag_vtx[i_part], n_vtx[i_part], "global tag_vtx :: ");
    }
  }



  int         **vtx_group_idx = NULL;
  PDM_g_num_t **vtx_group_gid = NULL;
  PDM_l_num_t **ridge_vtx_candidate   = NULL;
  if (geom_kind_tgt==PDM_GEOMETRY_KIND_RIDGE) {
    /**
     * If geom_kind_tgt is ridge, generate vtx group gid is not necessary,
     * we only need to get straddling vertices to generate edges containing them
     */
    PDM_malloc(ridge_vtx_candidate, pmn->n_part, PDM_l_num_t *);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {

      PDM_malloc(ridge_vtx_candidate[i_part], n_vtx[i_part], PDM_l_num_t);

      for(int i_vtx = 0; i_vtx < n_vtx[i_part]; ++i_vtx) {
        int beg = tag_vtx_idx[i_part][i_vtx  ];
        int end = tag_vtx_idx[i_part][i_vtx+1];
        if(end - beg > 0){
          int n_unique = PDM_inplace_unique(tag_vtx[i_part], beg, end-1);
          ridge_vtx_candidate[i_part][i_vtx] = n_unique > 1;
        }
      }
    }
  }
  else {
    _generate_group_gid(pmn->comm,
                        pmn->n_part,
                        n_vtx,
                        tag_vtx_idx,
                        tag_vtx,
                       &vtx_group_idx,
                       &vtx_group_gid);
  }
  for (int i_part=0; i_part<pmn->n_part; ++i_part) {
    PDM_free(tag_vtx_idx[i_part]);
    PDM_free(tag_vtx    [i_part]);
  }
  PDM_free(tag_vtx_idx);
  PDM_free(tag_vtx);


  printf("y'a peut etre un bout a pas faire en fonction de tgt=corner/ridge\n");



  if (geom_kind_tgt==PDM_GEOMETRY_KIND_CORNER) {

    if (pmn->corner!=NULL) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_compute_topo_corners : part_mesh_nodal already has corner section\n");
    }

    /**
     * Count global n_group while casting group id into integer
     * and transpose vtx_group for pmne storage
     */
    int **group_vtx_idx = NULL;
    int **group_vtx     = NULL;
    int  g_n_group = 0;
    _transpose_group_information(pmn->comm,
                                 pmn->n_part,
                                 n_vtx,
                                 vtx_group_idx,
                                 vtx_group_gid,
                                &g_n_group,
                                &group_vtx_idx,
                                &group_vtx);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(vtx_group_idx[i_part]);
      PDM_free(vtx_group_gid[i_part]);
    }
    PDM_free(vtx_group_idx);
    PDM_free(vtx_group_gid);



    /**
     * Generate global id for group entities, then store corners in pmne
     */
    PDM_g_num_t **vtx_gnum       = NULL;
    PDM_g_num_t **group_vtx_gnum = NULL;
    PDM_malloc(vtx_gnum, pmn->n_part, PDM_g_num_t*);
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      vtx_gnum[i_part] = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);
    }

    _generate_group_entity_gid(pmn->comm,
                               pmn->n_part,
                               // n_vtx,
                               vtx_gnum,
                               g_n_group,
                               group_vtx_idx,
                               group_vtx,
                              &group_vtx_gnum);



    pmn->corner = PDM_part_mesh_nodal_elmts_create(0, pmn->n_part, pmn->comm);
    int section_corner = PDM_part_mesh_nodal_elmts_add(pmn->corner, PDM_MESH_NODAL_POINT);
    PDM_part_mesh_nodal_elmts_n_group_set(pmn->corner, g_n_group);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_part_mesh_nodal_elmts_std_set(pmn->corner,
                                        section_corner,
                                        i_part,
                                        group_vtx_idx[i_part][g_n_group],
                                        group_vtx[i_part],
                                        NULL,
                                        NULL,
                                        NULL,
                                        PDM_OWNERSHIP_KEEP);
      for (int i_group=0; i_group<g_n_group; ++i_group) {
        int n_group_vtx = group_vtx_idx[i_part][i_group+1]-group_vtx_idx[i_part][i_group];
        int         *_group_vtx      = NULL;
        PDM_g_num_t *_group_vtx_gnum = NULL;
        PDM_malloc(_group_vtx     , n_group_vtx, int);
        PDM_malloc(_group_vtx_gnum, n_group_vtx, PDM_g_num_t);
        memcpy(_group_vtx     , &group_vtx     [i_part][group_vtx_idx[i_part][i_group]], n_group_vtx*sizeof(int        ));
        memcpy(_group_vtx_gnum, &group_vtx_gnum[i_part][group_vtx_idx[i_part][i_group]], n_group_vtx*sizeof(PDM_g_num_t));
        PDM_log_trace_array_int (_group_vtx     , n_group_vtx, "_group_vtx     ");
        PDM_log_trace_array_long(_group_vtx_gnum, n_group_vtx, "_group_vtx_gnum");
        PDM_part_mesh_nodal_elmts_group_set(pmn->corner, i_part,
                                            i_group, n_group_vtx,
                                            _group_vtx, _group_vtx_gnum,
                                            PDM_OWNERSHIP_KEEP);
      }
      PDM_free(group_vtx_idx [i_part]);
      PDM_free(group_vtx_gnum[i_part]);
    }
    PDM_free(group_vtx_gnum);
    PDM_free(group_vtx_idx);
    PDM_free(group_vtx);

  }
  else { //edge

    if (pmn->ridge!=NULL) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_compute_topo_corners : part_mesh_nodal already has ridge section\n");
    }

    /**
     * 
     * Go through edge element to see if their vertices are tagged by same group
     * 
     */

    int  *n_edge = NULL;
    int **entity1_edge_idx     = NULL;
    int **entity1_edge         = NULL;
    int **entity1_edge_vtx_idx = NULL;
    int **entity1_edge_vtx     = NULL;
    int **edge_parent          = NULL;
    int **edge_parent_pos      = NULL;
    PDM_part_mesh_nodal_elmts_sections_local_decompose_edges(pmne,
                                                             ridge_vtx_candidate,
                                                            &n_edge,
                                                            &entity1_edge_idx,
                                                            &entity1_edge_vtx_idx,
                                                            &entity1_edge_vtx,
                                                            &edge_parent,
                                                            &edge_parent_pos);

    if (debug_verbose==1) {
      for (int i_part=0; i_part<pmn->n_part; ++i_part) {
        log_trace("n_edge = %d\n", n_edge[i_part]);
        PDM_log_trace_array_int(entity1_edge_idx    [i_part], n_entity1[i_part], "entity1_edge_idx::");
        PDM_log_trace_array_int(entity1_edge_vtx_idx[i_part], entity1_edge_idx    [i_part][                         n_entity1[i_part]] , "entity1_edge_vtx_idx::");
        PDM_log_trace_array_int(entity1_edge_vtx    [i_part], entity1_edge_vtx_idx[i_part][entity1_edge_idx[i_part][n_entity1[i_part]]], "entity1_edge_vtx::");
        PDM_log_trace_array_int(edge_parent         [i_part], entity1_edge_idx    [i_part][                         n_entity1[i_part]] , "entity1_parent::");
        PDM_log_trace_array_int(edge_parent_pos     [i_part], entity1_edge_idx    [i_part][                         n_entity1[i_part]] , "entity1_parent_pos::");
      }
    }


    int  *n_unique_edge           = NULL;
    int **unique_edge_vtx_idx     = NULL;
    int **unique_edge_vtx         = NULL;
    int **unique_edge_to_edge_idx = NULL;
    int **unique_edge_to_edge     = NULL;
    PDM_malloc(n_unique_edge          , pmn->n_part, int  );
    PDM_malloc(unique_edge_vtx_idx    , pmn->n_part, int *);
    PDM_malloc(unique_edge_vtx        , pmn->n_part, int *);
    PDM_malloc(unique_edge_to_edge_idx, pmn->n_part, int *);
    PDM_malloc(unique_edge_to_edge    , pmn->n_part, int *);
    PDM_malloc(entity1_edge           , pmn->n_part, int *);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {

      entity1_edge[i_part] = PDM_array_new_idx_from_const_stride_int(1, n_edge[i_part]-1);
      PDM_malloc(unique_edge_vtx_idx[i_part],                              n_edge[i_part]+ 1, int);
      PDM_malloc(unique_edge_vtx    [i_part], entity1_edge_vtx_idx[i_part][n_edge[i_part]]  , int);
      unique_edge_vtx_idx[i_part][0] = 0;

      // Mysterious but still working
      PDM_part_mesh_nodal_elmts_generate_entity_connectivity(n_edge              [i_part],
                                                             entity1_edge_vtx_idx[i_part],
                                                             entity1_edge_vtx    [i_part],
                                                             edge_parent         [i_part],
                                                             0,
                                                             NULL,
                                                             NULL,
                                                            &unique_edge_to_edge_idx[i_part],
                                                            &unique_edge_to_edge    [i_part],
                                                             n_entity1              [i_part],
                                                             entity1_edge_idx       [i_part],
                                                             entity1_edge           [i_part],
                                                            &n_unique_edge          [i_part],
                                                             unique_edge_vtx_idx    [i_part],
                                                             unique_edge_vtx        [i_part],
                                                             PDM_TRUE);
    }

    PDM_g_num_t **edge_gid     = NULL;
    if (pmn->pcg[PDM_GEOMETRY_KIND_RIDGE]==NULL) {

      /**
       * Generate global ids for edges from their vertices in order to initialize a pcg
       * to gather group parent info over all procs
       */
      PDM_g_num_t **edge_vtx_gid = NULL;
      PDM_malloc(edge_gid    , pmn->n_part, PDM_g_num_t*);
      PDM_malloc(edge_vtx_gid, pmn->n_part, PDM_g_num_t*);
      for (int i_part=0; i_part<pmn->n_part; ++i_part) {

        PDM_malloc(edge_vtx_gid[i_part], 2*n_unique_edge[i_part], PDM_g_num_t);

        PDM_g_num_t *vtx_gid = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);

        for (int i_edge=0; i_edge<n_unique_edge[i_part]; ++i_edge) {
          int vtx1 = unique_edge_vtx[i_part][2*i_edge  ];
          int vtx2 = unique_edge_vtx[i_part][2*i_edge+1];
          edge_vtx_gid[i_part][2*i_edge  ] = vtx_gid[vtx1];
          edge_vtx_gid[i_part][2*i_edge+1] = vtx_gid[vtx2];
        }
      }


      PDM_gen_gnum_t* gen_edge_gid = PDM_gnum_create(3,
                                                      pmn->n_part,
                                                      PDM_TRUE,
                                                      1e-4,
                                                      pmn->comm,
                                                      PDM_OWNERSHIP_USER);
      PDM_gnum_set_parents_nuplet(gen_edge_gid, 2);
      for (int i_part=0; i_part<pmn->n_part; ++i_part) {
        PDM_gnum_set_from_parents(gen_edge_gid, i_part, n_unique_edge[i_part], edge_vtx_gid[i_part]);
      }

      PDM_gnum_compute(gen_edge_gid);

      for (int i_part=0; i_part<pmn->n_part; ++i_part) {
        edge_gid[i_part] = PDM_gnum_get(gen_edge_gid, i_part);
      }

      PDM_gnum_free(gen_edge_gid);


      pmn->ridge = PDM_part_mesh_nodal_elmts_create(1, pmn->n_part, pmn->comm);
      int section_ridge = PDM_part_mesh_nodal_elmts_add(pmn->ridge, PDM_MESH_NODAL_BAR2);

      for (int i_part=0; i_part<pmn->n_part; ++i_part) {
        PDM_part_mesh_nodal_elmts_std_set(pmn->ridge,
                                          section_ridge,
                                          i_part,
                                          n_unique_edge[i_part],
                                          unique_edge_vtx[i_part],
                                          edge_gid[i_part],
                                          NULL,
                                          NULL,
                                          PDM_OWNERSHIP_KEEP);
      }

      PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum(pmn, PDM_GEOMETRY_KIND_RIDGE);
    }



    if (debug_visu==1) {
      sprintf(filename, "edge_detected");
      PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_RIDGE, filename);
    }


    int **tag_edge_n       = NULL;
    int **tag_edge_idx     = NULL;
    int **tag_edge         = NULL;
    PDM_malloc(tag_edge_n  , pmn->n_part, int *);
    PDM_malloc(tag_edge_idx, pmn->n_part, int *);
    PDM_malloc(tag_edge    , pmn->n_part, int *);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {

      PDM_calloc(tag_edge_n[i_part], n_unique_edge[i_part], int);

      for(int i_group = 0; i_group < n_group_entity1; ++i_group) {
        int          n_entity1_group        = 0;
        int         *entity1_group          = NULL;
        PDM_g_num_t *entity1_group_ln_to_gn = NULL;
        PDM_part_mesh_nodal_elmts_group_get(pmne,
                                            i_part,
                                            i_group,
                                           &n_entity1_group,
                                           &entity1_group,
                                           &entity1_group_ln_to_gn,
                                            PDM_OWNERSHIP_BAD_VALUE);

        for (int i_entity_group=0; i_entity_group<n_entity1_group; ++i_entity_group) {
          int entity1 = entity1_group[i_entity_group]-1;
          for (int i_edge=entity1_edge_idx[i_part][entity1]; i_edge<entity1_edge_idx[i_part][entity1+1]; ++i_edge) {
            int edge = PDM_ABS(entity1_edge[i_part][i_edge])-1;
            tag_edge_n[i_part][edge]++;
          }
        }
      }

      PDM_malloc(tag_edge_idx[i_part], n_unique_edge[i_part]+1, int); tag_edge_idx[i_part][0] = 0;
      for(int i_vtx = 0; i_vtx < n_unique_edge[i_part]; ++i_vtx) {
        tag_edge_idx[i_part][i_vtx+1] = tag_edge_idx[i_part][i_vtx] + tag_edge_n[i_part][i_vtx];
        tag_edge_n  [i_part][i_vtx  ] = 0;
      }

      PDM_malloc(tag_edge[i_part], tag_edge_idx[i_part][n_unique_edge[i_part]], int);
      for(int i_group = 0; i_group < n_group_entity1; ++i_group) {
        int          n_entity1_group        = 0;
        int         *entity1_group          = NULL;
        PDM_g_num_t *entity1_group_ln_to_gn = NULL;
        PDM_part_mesh_nodal_elmts_group_get(pmne,
                                            i_part,
                                            i_group,
                                           &n_entity1_group,
                                           &entity1_group,
                                           &entity1_group_ln_to_gn,
                                            PDM_OWNERSHIP_BAD_VALUE);

        for (int i_entity_group=0; i_entity_group<n_entity1_group; ++i_entity_group) {
          int entity1 = entity1_group[i_entity_group]-1;
          for (int i_edge=entity1_edge_idx[i_part][entity1]; i_edge<entity1_edge_idx[i_part][entity1+1]; ++i_edge) {
            int edge = PDM_ABS(entity1_edge[i_part][i_edge])-1;
            int idx_write = tag_edge_idx[i_part][edge] + tag_edge_n[i_part][edge]++;
            tag_edge[i_part][idx_write] = i_group+1;
          }
        }
      }

      PDM_log_trace_connectivity_int(tag_edge_idx[i_part], tag_edge[i_part], n_unique_edge[i_part], "tag_edge (local) :: ");

    }

    int **_tag_edge_idx = NULL;
    int **_tag_edge     = NULL;
    PDM_part_comm_graph_gather_strided_int_data(pmn->pcg[PDM_GEOMETRY_KIND_RIDGE],
                                                n_unique_edge,
                                                tag_edge_idx,
                                                tag_edge,
                                              &_tag_edge_idx,
                                              &_tag_edge);
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(tag_edge_idx[i_part]);
      PDM_free(tag_edge    [i_part]);
    }
    PDM_free(tag_edge_idx);
    PDM_free(tag_edge);
    tag_edge_idx = _tag_edge_idx;
    tag_edge     = _tag_edge;

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_log_trace_connectivity_int(tag_edge_idx[i_part], tag_edge[i_part], n_unique_edge[i_part], "tag_edge (glbal) :: ");
    }

    int         **edge_group_idx = NULL;
    PDM_g_num_t **edge_group_gid = NULL;
    _generate_group_gid(pmn->comm,
                        pmn->n_part,
                        n_unique_edge,
                        tag_edge_idx,
                        tag_edge,
                       &edge_group_idx,
                       &edge_group_gid);
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_log_trace_connectivity_long(edge_group_idx[i_part], edge_group_gid[i_part], n_unique_edge[i_part], "edge_group_gid :: ");
      PDM_free(tag_edge_idx[i_part]);
      PDM_free(tag_edge    [i_part]);
    }
    PDM_free(tag_edge_idx);
    PDM_free(tag_edge);




    /**
     * All edge created arent ridges, so we need to extract them and regenerate a gnum
     */
    int          *n_ridge         = NULL;
    int         **ridge_vtx       = NULL;
    PDM_g_num_t **ridge_gid       = NULL;
    int         **ridge_group_idx = NULL;
    PDM_g_num_t **ridge_group_gid = NULL;
    PDM_malloc(n_ridge        , pmn->n_part, int          );
    PDM_malloc(ridge_vtx      , pmn->n_part, int         *);
    PDM_malloc(ridge_gid      , pmn->n_part, PDM_g_num_t *);
    PDM_malloc(ridge_group_idx, pmn->n_part, int         *);
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      n_ridge[i_part] = 0;
      for (int i_edge=0; i_edge<n_unique_edge[i_part]; ++i_edge) {
        if (edge_group_idx[i_part][i_edge+1]-edge_group_idx[i_part][i_edge]>0) {
          n_ridge[i_part]++;
        }
      }

      PDM_malloc(ridge_vtx      [i_part], 2*n_ridge[i_part]  , int        );
      PDM_malloc(ridge_gid      [i_part],   n_ridge[i_part]  , PDM_g_num_t);
      PDM_malloc(ridge_group_idx[i_part],   n_ridge[i_part]+1, int        ); ridge_group_idx[i_part][0] = 0;
      int i_write = 0;
      for (int i_edge=0; i_edge<n_unique_edge[i_part]; ++i_edge) {
        if (edge_group_idx[i_part][i_edge+1]-edge_group_idx[i_part][i_edge]>0) {
          ridge_vtx      [i_part][2*i_write+0] = unique_edge_vtx[i_part][2*i_edge+0];
          ridge_vtx      [i_part][2*i_write+1] = unique_edge_vtx[i_part][2*i_edge+1];
          ridge_gid      [i_part][  i_write  ] = edge_gid       [i_part][  i_edge  ];
          ridge_group_idx[i_part][  i_write+1] = ridge_group_idx[i_part][i_write]+1;
          i_write++;
        }
      }
      PDM_log_trace_array_long(ridge_gid[i_part], n_ridge[i_part], "ridge_gid_parent :: ");
    }
    PDM_gen_gnum_t* gen_ridge_gid = PDM_gnum_create(3,
                                                    pmn->n_part,
                                                    PDM_TRUE,
                                                    1e-4,
                                                    pmn->comm,
                                                    PDM_OWNERSHIP_USER);
    PDM_gnum_set_parents_nuplet(gen_ridge_gid, 1);
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_gnum_set_from_parents(gen_ridge_gid, i_part, n_ridge[i_part], ridge_gid[i_part]);
    }

    PDM_gnum_compute(gen_ridge_gid);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(ridge_gid[i_part]);
      ridge_gid[i_part] = PDM_gnum_get(gen_ridge_gid, i_part);
      PDM_log_trace_array_long(ridge_gid[i_part], n_ridge[i_part], "ridge_gid :: ");

    }

    PDM_gnum_free(gen_ridge_gid);
    ridge_group_gid = edge_group_gid;


    /**
     * Count global n_group while casting group id into integer
     * and transpose vtx_group for pmne storage
     */
    int **group_ridge_idx = NULL;
    int **group_ridge     = NULL;
    int  g_n_group = 0;
    _transpose_group_information(pmn->comm,
                                 pmn->n_part,
                                 n_ridge,
                                 ridge_group_idx,
                                 ridge_group_gid,
                                &g_n_group,
                                &group_ridge_idx,
                                &group_ridge);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(ridge_group_idx[i_part]);
      PDM_free(ridge_group_gid[i_part]);
    }
    PDM_free(ridge_group_idx);
    PDM_free(ridge_group_gid);


    /**
     * Generate global id for group entities, then store corners in pmne
     */

    PDM_g_num_t **group_ridge_gid = NULL;
    _generate_group_entity_gid(pmn->comm,
                               pmn->n_part,
                               ridge_gid,
                               g_n_group,
                               group_ridge_idx,
                               group_ridge,
                              &group_ridge_gid);




    PDM_part_mesh_nodal_elmts_free(pmn->ridge);

    pmn->ridge = PDM_part_mesh_nodal_elmts_create(1, pmn->n_part, pmn->comm);
    int section_ridge = PDM_part_mesh_nodal_elmts_add(pmn->ridge, PDM_MESH_NODAL_BAR2);
    PDM_part_mesh_nodal_elmts_n_group_set(pmn->ridge, g_n_group);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_part_mesh_nodal_elmts_std_set(pmn->ridge,
                                        section_ridge,
                                        i_part,
                                        n_ridge[i_part],
                                        ridge_vtx[i_part],
                                        ridge_gid[i_part],
                                        NULL,
                                        NULL,
                                        PDM_OWNERSHIP_KEEP);
      for (int i_group=0; i_group<g_n_group; ++i_group) {
        int n_group_ridge = group_ridge_idx[i_part][i_group+1]-group_ridge_idx[i_part][i_group];
        int         *_group_ridge      = NULL;
        PDM_g_num_t *_group_ridge_gnum = NULL;
        PDM_malloc(_group_ridge     , n_group_ridge, int);
        PDM_malloc(_group_ridge_gnum, n_group_ridge, PDM_g_num_t);
        memcpy(_group_ridge     , &group_ridge     [i_part][group_ridge_idx[i_part][i_group]], n_group_ridge*sizeof(int        ));
        memcpy(_group_ridge_gnum, &group_ridge_gid[i_part][group_ridge_idx[i_part][i_group]], n_group_ridge*sizeof(PDM_g_num_t));
        PDM_log_trace_array_int (_group_ridge     , n_group_ridge, "_group_ridge     ");
        PDM_log_trace_array_long(_group_ridge_gnum, n_group_ridge, "_group_ridge_gnum");
        PDM_part_mesh_nodal_elmts_group_set(pmn->ridge, i_part,
                                            i_group, n_group_ridge,
                                            _group_ridge, _group_ridge_gnum,
                                            PDM_OWNERSHIP_KEEP);
      }
      PDM_free(group_ridge_idx [i_part]);
      PDM_free(group_ridge_gid[i_part]);

      sprintf(filename, "ridge_final");
      PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_RIDGE, filename);
    }
    PDM_free(group_ridge_gid);
    PDM_free(group_ridge_idx);
    PDM_free(group_ridge);


  }



  // for (int i_part=0; i_part<pmn->n_part; ++i_part) {
  //   PDM_free(vtx_parent_group_idx[i_part]);
  //   PDM_free(vtx_parent_group    [i_part]);
  // }
  // PDM_free(vtx_parent_group_idx);
  // PDM_free(vtx_parent_group);

  // for (int i_part=0; i_part<pmn->n_part; ++i_part) {
  //   PDM_free(vtx_group_gid[i_part]);
  // }
  // PDM_free(vtx_group_gid);

  for (int i_part=0; i_part<pmn->n_part; ++i_part) {
    PDM_free(tag_vtx_n[i_part]);
  }
  PDM_free(tag_vtx_n);

  if (geom_kind_tgt==PDM_GEOMETRY_KIND_RIDGE) {
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(ridge_vtx_candidate[i_part]);
    }
    PDM_free(ridge_vtx_candidate);
  }

  PDM_free(n_entity1);
  PDM_free(n_vtx);

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
