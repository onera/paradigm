#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_printf.h"
#include "pdm_gnum.h"
#include "pdm_dcube_gen.h"
#include "pdm_part.h"
#include "pdm_extract_part.h"
#include "pdm_sort.h"
#include "pdm_part_to_part.h"
#include "pdm_array.h"
#include "pdm_closest_points.h"
#include "pdm_priv.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_binary_search.h"
#include "pdm_para_octree.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"

#include "pdm_mpi.h"
#include "pdm_config.h"

static void
_compute_gnum_in_place
(
  const int           n_part,
  const int          *n_entity,
        PDM_g_num_t **entity_ln_to_gn,
        PDM_MPI_Comm  comm
)
{

  PDM_gen_gnum_t *gnum = PDM_gnum_create(3,
                                         n_part,
                                         PDM_FALSE,
                                         1.,
                                         comm,
                                         PDM_OWNERSHIP_USER);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_gnum_set_from_parents(gnum,
                              i_part,
                              n_entity[i_part],
                              entity_ln_to_gn[i_part]);
  }

  PDM_gnum_compute(gnum);

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(entity_ln_to_gn[i_part]);
    entity_ln_to_gn[i_part] = PDM_gnum_get(gnum, i_part);
  }

  PDM_gnum_free(gnum);

}

static void
_leaf_get
(
  const int   n_explicit_nodes,
  const int  *leaf_id,
        int  *n_leaf,
        int **leaf_nodes_id
)
{

  int _n_leaf = 0;

  int *_leaf_nodes_id = malloc(n_explicit_nodes * sizeof(int));
  for(int i_node = 0; i_node < n_explicit_nodes; i_node++) {
    if(leaf_id[i_node] != -1){
      _leaf_nodes_id[_n_leaf++] = i_node;
    }
  }
  _leaf_nodes_id = realloc(_leaf_nodes_id, _n_leaf * sizeof(int));

  *n_leaf        = _n_leaf;
  *leaf_nodes_id = _leaf_nodes_id;

}

static void
_mean_array_double_per_leaf
(
  const int     n_leaf,
  const int    *leaf_nodes_id,
  const int    *n_points,
  const int    *range,
  const int     stride,
  const double *pts_array_double,
        double *leaf_array_double
)
{

  for (int i_leaf = 0; i_leaf < n_leaf; i_leaf++) {
    for (int i = 0; i < stride; i++) {
      leaf_array_double[stride*i_leaf+i] = 0.0;
    }
    for (int i_pts = range[leaf_nodes_id[i_leaf]]; i_pts < range[leaf_nodes_id[i_leaf]] + n_points[leaf_nodes_id[i_leaf]]; i_pts++) {
      for (int i = 0; i < stride; i++) {
        leaf_array_double[stride*i_leaf+i] = leaf_array_double[stride*i_leaf+i] + pts_array_double[stride*i_pts+i];
      }
    }
    for (int i = 0; i < stride; i++) {
      leaf_array_double[stride*i_leaf+i] /= n_points[leaf_nodes_id[i_leaf]];
    }
  }

}

typedef struct _PDM_surf_deform_t {

  int                  n_part;
  int                 *n_face;
  int                 *n_vtx;
  double             **coords;
  PDM_g_num_t        **vtx_ln_to_gn;
  PDM_g_num_t        **face_ln_to_gn;

  int                **face_vtx_idx;
  int                **face_vtx;

  double             **dcoords;
  double             **aux_geom;

  PDM_para_octree_t  **octree_layer;
  PDM_part_to_part_t **ptp_layer;

} PDM_surf_deform_t;

typedef struct _PDM_cloud_deform_t {

  int                  n_part;
  int                 *n_points;
  double             **coords;
  double             **dcoords;
  PDM_g_num_t        **gnum;

  PDM_part_to_block_t *ptb_blk;
  PDM_part_to_part_t  *ptp_blk;
  int                  blk_n_points;
  double              *blk_coords;
  double              *blk_dcoords;
  int                 *blk_buffer_from_surf_idx;
  int                 *blk_buffer_from_surf;
  double             **blk_coords_from_surf;
  double             **blk_dcoords_from_surf;
  double             **blk_aux_geom_from_surf;

} PDM_cloud_deform_t;

typedef struct _PDM_mesh_deform_t {

  PDM_MPI_Comm        comm;                    /*!< MPI communicator */
  int                 is_cpt;                  /*!< Construction flag */
  int                 is_blk;                  /*!< Construction flag */
  int                 is_res;                  /*!< Construction flag */
  int                 n_aux_geom;
  int                 n_layer;
  int                *n_leaf_per_layer;
  double              min_dist_d;
  int                 min_dist_n_vtx;

  PDM_surf_deform_t  *surf_deform;
  PDM_cloud_deform_t *cloud_deform;

} PDM_mesh_deform_t;

static PDM_mesh_deform_t*
_PDM_mesh_deform_create
(
  const int           n_part_surf,
  const int           n_part_cloud,
  const int           n_aux_geom,
  const int           n_layer,
  const int          *n_leaf_per_layer,
  const double        min_dist_d,
  const int           min_dist_n_vtx,
  const PDM_MPI_Comm  comm
)
{

  PDM_mesh_deform_t *def = (PDM_mesh_deform_t *) malloc(sizeof(PDM_mesh_deform_t));

  def->comm             = comm;
  def->is_cpt           = 0;
  def->is_blk           = 0;
  def->is_res           = 0;
  def->n_aux_geom       = n_aux_geom;
  def->n_layer          = n_layer;
  def->n_leaf_per_layer = (int*) n_leaf_per_layer;
  def->min_dist_d       = min_dist_d;
  def->min_dist_n_vtx   = min_dist_n_vtx;

  def->surf_deform  = (PDM_surf_deform_t  *) malloc(sizeof(PDM_surf_deform_t ));
  def->cloud_deform = (PDM_cloud_deform_t *) malloc(sizeof(PDM_cloud_deform_t));

  PDM_surf_deform_t  *_surf_deform  = def->surf_deform;
  PDM_cloud_deform_t *_cloud_deform = def->cloud_deform;

  _surf_deform->n_part        = n_part_surf;
  _surf_deform->n_face        = malloc(sizeof(int          ) * (n_part_surf+n_layer));
  _surf_deform->n_vtx         = malloc(sizeof(int          ) * (n_part_surf+n_layer));
  _surf_deform->coords        = malloc(sizeof(double      *) * (n_part_surf+n_layer));
  _surf_deform->vtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t *) * (n_part_surf+n_layer));
  _surf_deform->face_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * (n_part_surf+n_layer));
  _surf_deform->face_vtx_idx  = malloc(sizeof(int         *) * (n_part_surf+n_layer));
  _surf_deform->face_vtx      = malloc(sizeof(int         *) * (n_part_surf+n_layer));
  _surf_deform->dcoords       = malloc(sizeof(double      *) * (n_part_surf+n_layer));
  _surf_deform->aux_geom      = malloc(sizeof(double      *) * (n_part_surf+n_layer));

  for (int i = 0; i < n_part_surf+n_layer; i++) {
    _surf_deform->n_face[i]        = 0;
    _surf_deform->n_vtx[i]         = 0;
    _surf_deform->coords[i]        = NULL;
    _surf_deform->vtx_ln_to_gn[i]  = NULL;
    _surf_deform->face_ln_to_gn[i] = NULL;
    _surf_deform->face_vtx_idx[i]  = NULL;
    _surf_deform->face_vtx[i]      = NULL;
    _surf_deform->dcoords[i]       = NULL;
    _surf_deform->aux_geom[i]      = NULL;
  }

  _cloud_deform->n_part   = n_part_cloud;
  _cloud_deform->n_points = malloc(sizeof(int          ) * n_part_cloud);
  _cloud_deform->coords   = malloc(sizeof(double      *) * n_part_cloud);
  _cloud_deform->gnum     = malloc(sizeof(PDM_g_num_t *) * n_part_cloud);

  for (int i = 0; i < n_part_cloud; i++) {
    _cloud_deform->n_points[i] = 0;
    _cloud_deform->coords[i]   = NULL;
    _cloud_deform->gnum[i]     = NULL;
  }

  _cloud_deform->blk_n_points             = 0;
  _cloud_deform->blk_coords               = NULL;
  _cloud_deform->blk_dcoords              = NULL;
  _cloud_deform->blk_coords_from_surf     = NULL;
  _cloud_deform->blk_dcoords_from_surf    = NULL;
  _cloud_deform->blk_aux_geom_from_surf   = NULL;
  _cloud_deform->blk_buffer_from_surf_idx = NULL;
  _cloud_deform->blk_buffer_from_surf     = NULL;

  return def;

}

static void
_PDM_mesh_deform_surf_part_set
(
  PDM_mesh_deform_t *def,
  int                i_part,
  int                n_face,
  int                n_vtx,
  int               *face_vtx_idx,
  int               *face_vtx,
  PDM_g_num_t       *face_ln_to_gn,
  PDM_g_num_t       *vtx_ln_to_gn,
  double            *vtx_coord,
  double            *vtx_dcoord,
  double            *vtx_aux_geom
)
{

  PDM_surf_deform_t  *_surf_deform  = def->surf_deform;

  _surf_deform->n_face[i_part]        = n_face;
  _surf_deform->n_vtx[i_part]         = n_vtx;
  _surf_deform->coords[i_part]        = vtx_coord;
  _surf_deform->vtx_ln_to_gn[i_part]  = vtx_ln_to_gn;
  _surf_deform->face_ln_to_gn[i_part] = face_ln_to_gn;
  _surf_deform->face_vtx_idx[i_part]  = face_vtx_idx;
  _surf_deform->face_vtx[i_part]      = face_vtx;
  _surf_deform->dcoords[i_part]       = vtx_dcoord;
  _surf_deform->aux_geom[i_part]      = vtx_aux_geom;

}

static void
_PDM_mesh_deform_cloud_part_set
(
  PDM_mesh_deform_t *def,
  int                i_part,
  int                n_points,
  PDM_g_num_t       *gnum,
  double            *coord
)
{

  PDM_cloud_deform_t *_cloud_deform = def->cloud_deform;

  _cloud_deform->n_points[i_part] = n_points;
  _cloud_deform->coords[i_part]   = coord;
  _cloud_deform->gnum[i_part]     = gnum;

}

static void
_PDM_mesh_deform_partial_free
(
  PDM_mesh_deform_t *def
)
{

  if (def->is_cpt == 0) {
    return;
  }

  PDM_surf_deform_t  *_surf_deform  = def->surf_deform;
  PDM_cloud_deform_t *_cloud_deform = def->cloud_deform;

  for (int i_layer = 0; i_layer < def->n_layer; i_layer++) {
    free(_surf_deform->vtx_ln_to_gn[_surf_deform->n_part+i_layer]);
    free(_surf_deform->coords      [_surf_deform->n_part+i_layer]);
    free(_surf_deform->dcoords     [_surf_deform->n_part+i_layer]);
    if (def->n_aux_geom > 0) {
      free(_surf_deform->aux_geom[_surf_deform->n_part+i_layer]);
    }
    PDM_part_to_part_free(_surf_deform->ptp_layer[i_layer]);
    PDM_para_octree_free(_surf_deform->octree_layer[i_layer]);
  }

  free(_surf_deform->ptp_layer   );
  free(_surf_deform->octree_layer);

  free(_cloud_deform->blk_buffer_from_surf_idx);
  free(_cloud_deform->blk_buffer_from_surf    );
  PDM_part_to_block_free(_cloud_deform->ptb_blk);
  PDM_part_to_part_free(_cloud_deform->ptp_blk);

  def->is_cpt = 0;

}

static void
_PDM_mesh_deform_free
(
  PDM_mesh_deform_t *def
)
{

  _PDM_mesh_deform_partial_free(def);

  PDM_surf_deform_t  *_surf_deform  = def->surf_deform;
  PDM_cloud_deform_t *_cloud_deform = def->cloud_deform;

  free(_surf_deform->n_face       );
  free(_surf_deform->n_vtx        );
  free(_surf_deform->coords       );
  free(_surf_deform->vtx_ln_to_gn );
  free(_surf_deform->face_ln_to_gn);
  free(_surf_deform->face_vtx_idx );
  free(_surf_deform->face_vtx     );
  free(_surf_deform->dcoords      );
  free(_surf_deform->aux_geom     );

  free(_cloud_deform->n_points);
  free(_cloud_deform->coords  );
  free(_cloud_deform->gnum    );

  if (def->is_blk == 1) {
    free(_cloud_deform->blk_coords );
    free(_cloud_deform->blk_dcoords);
    free(_cloud_deform->blk_coords_from_surf[0] );
    free(_cloud_deform->blk_dcoords_from_surf[0]);
    free(_cloud_deform->blk_coords_from_surf    );
    free(_cloud_deform->blk_dcoords_from_surf   );
    if (def->n_aux_geom > 0) {
      free(_cloud_deform->blk_aux_geom_from_surf[0]);
      free(_cloud_deform->blk_aux_geom_from_surf   );
    }
    def->is_blk = 0;
  }

  if (def->is_res == 1) {
    for (int i_part = 0; i_part < _cloud_deform->n_part; i_part++) {
      free(_cloud_deform->dcoords[i_part]);
    }
    free(_cloud_deform->dcoords);
    def->is_res = 0;
  }

  free(def->surf_deform );
  free(def->cloud_deform);
  free(def);

}

static void
_PDM_mesh_deform_compute
(
  PDM_mesh_deform_t *def
)
{

  _PDM_mesh_deform_partial_free(def);

  PDM_surf_deform_t  *_surf_deform  = def->surf_deform;
  PDM_cloud_deform_t *_cloud_deform = def->cloud_deform;

  int n_rank;
  PDM_MPI_Comm_size(def->comm, &n_rank);
  int i_rank;
  PDM_MPI_Comm_rank(def->comm, &i_rank);

  PDM_g_num_t _max_g_num = 0;
  PDM_g_num_t  max_g_num = 0;

  for (int i_part = 0; i_part < _surf_deform->n_part; i_part++) {
    for (int i_vtx = 0; i_vtx < _surf_deform->n_vtx[i_part]; i_vtx++) {
      _max_g_num = PDM_MAX(_max_g_num, _surf_deform->vtx_ln_to_gn[i_part][i_vtx]);
    }
  }

  PDM_MPI_Allreduce(&_max_g_num, &max_g_num, 1, PDM_MPI_LONG, PDM_MPI_MAX, def->comm);

  _max_g_num = max_g_num + 1;

  _surf_deform->octree_layer = malloc(sizeof(PDM_para_octree_t  *) * def->n_layer);
  _surf_deform->ptp_layer    = malloc(sizeof(PDM_part_to_part_t *) * def->n_layer);

  PDM_g_num_t **distri_layer = malloc(sizeof(PDM_g_num_t *) * def->n_layer);
  PDM_g_num_t  *gnum_layer   = malloc(sizeof(PDM_g_num_t  ) * def->n_layer);

  int depth_max = 50;

  for (int i_layer = 0; i_layer < def->n_layer; i_layer++) {

    _surf_deform->octree_layer[i_layer] = PDM_para_octree_create(1,
                                                                 depth_max,
                                                                 def->n_leaf_per_layer[i_layer],
                                                                 0,
                                                                 def->comm);

    for (int i_part = 0; i_part < _surf_deform->n_part; i_part++) {

      PDM_para_octree_point_cloud_set(_surf_deform->octree_layer[i_layer],
                                      i_part,
                                      _surf_deform->n_vtx       [i_part],
                                      _surf_deform->coords      [i_part],
                                      _surf_deform->vtx_ln_to_gn[i_part]);

    }

    PDM_para_octree_build(_surf_deform->octree_layer[i_layer], NULL);

    int          n_pts_octree     = 0;
    double      *pts_coord_octree = NULL;
    PDM_g_num_t *pts_gnum_octree  = NULL;

    PDM_para_octree_points_get(_surf_deform->octree_layer[i_layer],
                              &n_pts_octree,
                              &pts_coord_octree,
                              &pts_gnum_octree);

    int *pts_gnum_octree_idx = PDM_array_new_idx_from_const_stride_int(1, n_pts_octree);

    _surf_deform->ptp_layer[i_layer] = PDM_part_to_part_create((const PDM_g_num_t **) &pts_gnum_octree,
                                                               (const int          *) &n_pts_octree,
                                                                                       1,
                                                               (const PDM_g_num_t **)  _surf_deform->vtx_ln_to_gn,
                                                               (const int          *)  _surf_deform->n_vtx,
                                                                                       _surf_deform->n_part,
                                                               (const int         **) &pts_gnum_octree_idx,
                                                               (const PDM_g_num_t **) &pts_gnum_octree,
                                                                                       def->comm);

    free(pts_gnum_octree_idx);

    int  n_explicit_nodes = 0;
    int *n_points         = NULL;
    int *range            = NULL;
    int *leaf_id          = NULL;
    int *leaf_nodes_id    = NULL;
    int *children_id      = NULL;
    int *ancestor_id      = NULL;
    int  n_child          = 0;
    int  n_leaf           = 0;
    int  stack_size       = 0;

    PDM_para_octree_explicit_node_get(_surf_deform->octree_layer[i_layer],
                                     &n_explicit_nodes,
                                     &n_points,
                                     &range,
                                     &leaf_id,
                                     &children_id,
                                     &ancestor_id,
                                     &n_child,
                                     &stack_size);

    _leaf_get(n_explicit_nodes,
              leaf_id,
             &n_leaf,
             &leaf_nodes_id);

    distri_layer[i_layer] = malloc (sizeof(PDM_g_num_t) * (n_rank+1));
    PDM_distrib_compute(n_leaf, distri_layer[i_layer], -1, def->comm);

    gnum_layer[i_layer] = _max_g_num;

    _surf_deform->n_vtx       [_surf_deform->n_part+i_layer] = n_leaf;
    _surf_deform->vtx_ln_to_gn[_surf_deform->n_part+i_layer] = malloc (sizeof(PDM_g_num_t)     * n_leaf);
    _surf_deform->coords      [_surf_deform->n_part+i_layer] = malloc (sizeof(double     ) * 3 * n_leaf);
    _surf_deform->dcoords     [_surf_deform->n_part+i_layer] = malloc (sizeof(double     ) * 3 * n_leaf);
    if (def->n_aux_geom > 0) {
      _surf_deform->aux_geom[_surf_deform->n_part+i_layer] = malloc (sizeof(double) * def->n_aux_geom * n_leaf);
    }

    for (int i_vtx = 0; i_vtx < n_leaf; i_vtx++) {
      _surf_deform->vtx_ln_to_gn[_surf_deform->n_part+i_layer][i_vtx] = _max_g_num + distri_layer[i_layer][i_rank] + i_vtx;
    }

    _max_g_num = _max_g_num + distri_layer[i_layer][n_rank];

    free(leaf_nodes_id);

  }

  PDM_dist_cloud_surf_t *dcs = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED,
                                                          1,
                                                          def->comm,
                                                          PDM_OWNERSHIP_KEEP);

  PDM_dist_cloud_surf_surf_mesh_global_data_set(dcs,
                                                _surf_deform->n_part);

  PDM_dist_cloud_surf_n_part_cloud_set(dcs,
                                       0,
                                       _cloud_deform->n_part);

  for (int i_part = 0; i_part < _surf_deform->n_part; i_part++) {

    PDM_dist_cloud_surf_surf_mesh_part_set(dcs,
                                           i_part,
                                           _surf_deform->n_face       [i_part],
                                           _surf_deform->face_vtx_idx [i_part],
                                           _surf_deform->face_vtx     [i_part],
                                           _surf_deform->face_ln_to_gn[i_part],
                                           _surf_deform->n_vtx        [i_part],
                                           _surf_deform->coords       [i_part],
                                           _surf_deform->vtx_ln_to_gn [i_part]);

  }

  for (int i_part = 0; i_part < _cloud_deform->n_part; i_part++) {

    PDM_dist_cloud_surf_cloud_set(dcs,
                                  0,
                                  i_part,
                                  _cloud_deform->n_points[i_part],
                                  _cloud_deform->coords  [i_part],
                                  _cloud_deform->gnum    [i_part]);

  }

  PDM_dist_cloud_surf_compute(dcs);

  double _max_dist = 0.0;
  double  max_dist = 0.0;

  for (int i_part = 0; i_part < _cloud_deform->n_part; i_part++) {

    PDM_g_num_t *surf_gnum = NULL;
    double      *surf_dist = NULL;
    double      *surf_proj = NULL;

    PDM_dist_cloud_surf_get(dcs,
                            0,
                            i_part,
                           &surf_dist,
                           &surf_proj,
                           &surf_gnum);

    for (int i_pts = 0; i_pts < _cloud_deform->n_points[i_part]; i_pts++) {
      _max_dist = PDM_MAX(_max_dist, surf_dist[i_pts]);
    }

  }

  PDM_MPI_Allreduce(&_max_dist, &max_dist, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, def->comm);

  int          *n_min_dist        = malloc (sizeof(int          ) * _cloud_deform->n_part);
  PDM_g_num_t **min_dist_gnum     = malloc (sizeof(PDM_g_num_t *) * _cloud_deform->n_part);
  double      **min_dist_coords   = malloc (sizeof(double      *) * _cloud_deform->n_part);
  int         **cloud_to_min_dist = malloc (sizeof(int         *) * _cloud_deform->n_part);

  int         **n_cloud_to_surf   = malloc (sizeof(int         *) * _cloud_deform->n_part);
  int         **cloud_to_surf_idx = malloc (sizeof(int         *) * _cloud_deform->n_part);
  PDM_g_num_t **cloud_to_surf     = malloc (sizeof(PDM_g_num_t *) * _cloud_deform->n_part);
  double      **cloud_weight      = malloc (sizeof(double      *) * _cloud_deform->n_part);

  for (int i_part = 0; i_part < _cloud_deform->n_part; i_part++) {

    n_min_dist       [i_part] = 0;
    min_dist_gnum    [i_part] = malloc (sizeof(PDM_g_num_t)     * _cloud_deform->n_points[i_part]);
    min_dist_coords  [i_part] = malloc (sizeof(double     ) * 3 * _cloud_deform->n_points[i_part]);
    cloud_to_min_dist[i_part] = malloc (sizeof(int        )     * _cloud_deform->n_points[i_part]);

    n_cloud_to_surf[i_part] = malloc (sizeof(int   ) * _cloud_deform->n_points[i_part]);
    cloud_weight   [i_part] = malloc (sizeof(double) * _cloud_deform->n_points[i_part]);

    PDM_g_num_t *sorted_gnum = malloc (sizeof(PDM_g_num_t) * _cloud_deform->n_points[i_part]);
    int         *order       = malloc (sizeof(int        ) * _cloud_deform->n_points[i_part]);

    PDM_g_num_t *surf_gnum = NULL;
    double      *surf_dist = NULL;
    double      *surf_proj = NULL;

    PDM_dist_cloud_surf_get(dcs,
                            0,
                            i_part,
                           &surf_dist,
                           &surf_proj,
                           &surf_gnum);

    for (int i_pts = 0; i_pts < _cloud_deform->n_points[i_part]; i_pts++) {

      n_cloud_to_surf[i_part][i_pts] = 0;

      if (surf_dist[i_pts] <= def->min_dist_d) {
        n_cloud_to_surf[i_part][i_pts] = def->min_dist_n_vtx;
        cloud_weight   [i_part][i_pts] = (double) def->min_dist_n_vtx;
        for (int j_vtx = 0; j_vtx < n_min_dist[i_part]; j_vtx++) {
          sorted_gnum[j_vtx] = min_dist_gnum[i_part][j_vtx];
          order      [j_vtx] = j_vtx;
        }
        PDM_sort_long(sorted_gnum, order, n_min_dist[i_part]);
        int i_gnum = -1;
        if (n_min_dist[i_part] > 0) {
          i_gnum = PDM_binary_search_long(surf_gnum[i_pts], sorted_gnum, n_min_dist[i_part]);
        }
        if (i_gnum < 0) {
          min_dist_gnum    [i_part][  n_min_dist[i_part]    ] = surf_gnum[i_pts];
          min_dist_coords  [i_part][3*n_min_dist[i_part]    ] = _cloud_deform->coords[i_part][3*i_pts    ];
          min_dist_coords  [i_part][3*n_min_dist[i_part] + 1] = _cloud_deform->coords[i_part][3*i_pts + 1];
          min_dist_coords  [i_part][3*n_min_dist[i_part] + 2] = _cloud_deform->coords[i_part][3*i_pts + 2];
          cloud_to_min_dist[i_part][i_pts] = n_min_dist[i_part];
          n_min_dist[i_part]++;
        } else {
          cloud_to_min_dist[i_part][i_pts] = order[i_gnum];
        }
      } else {
        for (int i_layer = 0; i_layer < def->n_layer; i_layer++) {
          double a_dist = def->min_dist_d + (max_dist - def->min_dist_d)/def->n_layer* i_layer;
          double b_dist = def->min_dist_d + (max_dist - def->min_dist_d)/def->n_layer*(i_layer+1);
          if (surf_dist[i_pts] > a_dist && surf_dist[i_pts] <= b_dist) {
            n_cloud_to_surf  [i_part][i_pts] = 1;
            cloud_weight     [i_part][i_pts] = (double) distri_layer[i_layer][n_rank];
            cloud_to_min_dist[i_part][i_pts] = i_layer;
          }
        }
      }

      assert(n_cloud_to_surf[i_part][i_pts] != 0);

    }

    min_dist_gnum  [i_part] = realloc (min_dist_gnum  [i_part], sizeof(PDM_g_num_t)     * n_min_dist[i_part]);
    min_dist_coords[i_part] = realloc (min_dist_coords[i_part], sizeof(double     ) * 3 * n_min_dist[i_part]);

    cloud_to_surf_idx[i_part] = PDM_array_new_idx_from_sizes_int(n_cloud_to_surf[i_part], _cloud_deform->n_points[i_part]);
    cloud_to_surf    [i_part] = malloc (sizeof(PDM_g_num_t) * cloud_to_surf_idx[i_part][_cloud_deform->n_points[i_part]]);

    free(sorted_gnum);
    free(order);

  }

  _compute_gnum_in_place(_cloud_deform->n_part,
                         n_min_dist,
                         min_dist_gnum,
                         def->comm);

  PDM_closest_point_t* cls = PDM_closest_points_create(def->comm,
                                                       def->min_dist_n_vtx,
                                                       PDM_OWNERSHIP_KEEP);

  PDM_closest_points_n_part_cloud_set(cls,
                                      _surf_deform->n_part,
                                      _cloud_deform->n_part);

  for (int i_part = 0; i_part < _surf_deform->n_part; i_part++) {

    PDM_closest_points_src_cloud_set(cls,
                                     i_part,
                                     _surf_deform->n_vtx        [i_part],
                                     _surf_deform->coords       [i_part],
                                     _surf_deform->vtx_ln_to_gn [i_part]);

  }

  for (int i_part = 0; i_part < _cloud_deform->n_part; i_part++) {

    PDM_closest_points_tgt_cloud_set(cls,
                                     i_part,
                                     n_min_dist     [i_part],
                                     min_dist_coords[i_part],
                                     min_dist_gnum  [i_part]);

  }

  PDM_closest_points_compute(cls);

  for (int i_part = 0; i_part < _cloud_deform->n_part; i_part++) {

    PDM_g_num_t *surf_gnum = NULL;
    double      *surf_dist = NULL;
    double      *surf_proj = NULL;

    PDM_dist_cloud_surf_get(dcs,
                            0,
                            i_part,
                           &surf_dist,
                           &surf_proj,
                           &surf_gnum);

    PDM_g_num_t *closest_src_gnum = NULL;
    double      *closest_src_dist = NULL;

    PDM_closest_points_get(cls,
                           i_part,
                          &closest_src_gnum,
                          &closest_src_dist);

    for (int i_pts = 0; i_pts < _cloud_deform->n_points[i_part]; i_pts++) {

      int idx_write = 0;

      if (surf_dist[i_pts] <= def->min_dist_d) {
        idx_write = cloud_to_surf_idx[i_part][i_pts];
        for (int j_vtx = 0; j_vtx < n_cloud_to_surf[i_part][i_pts]; j_vtx++) {
          cloud_to_surf[i_part][idx_write++] = closest_src_gnum[def->min_dist_n_vtx*cloud_to_min_dist[i_part][i_pts]+j_vtx];
        }
      } else {
        idx_write = cloud_to_surf_idx[i_part][i_pts];
        cloud_to_surf[i_part][idx_write] = -(cloud_to_min_dist[i_part][i_pts]+1);
      }

    }

  }

  PDM_dist_cloud_surf_free(dcs);
  PDM_closest_points_free(cls);

  for (int i_part = 0; i_part < _cloud_deform->n_part; i_part++) {
    free(min_dist_gnum    [i_part]);
    free(min_dist_coords  [i_part]);
    free(cloud_to_min_dist[i_part]);
  }

  free(n_min_dist       );
  free(min_dist_gnum    );
  free(min_dist_coords  );
  free(cloud_to_min_dist);

  _cloud_deform->ptb_blk = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                    PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                    1.,
                                                    _cloud_deform->gnum,
                                                    cloud_weight,
                                                    _cloud_deform->n_points,
                                                    _cloud_deform->n_part,
                                                    def->comm);

  _cloud_deform->blk_n_points = PDM_part_to_block_n_elt_block_get(_cloud_deform->ptb_blk);

  int         *blk_n_cloud_to_surf = NULL;
  PDM_g_num_t *blk_cloud_to_surf   = NULL;

  int blk_size = PDM_part_to_block_exch(_cloud_deform->ptb_blk,
                                        sizeof(PDM_g_num_t),
                                        PDM_STRIDE_VAR_INTERLACED,
                                        1,
                                        n_cloud_to_surf,
                             (void **)  cloud_to_surf,
                                       &blk_n_cloud_to_surf,
                             (void **) &blk_cloud_to_surf);

  PDM_UNUSED(blk_size);

  int *blk_cloud_to_surf_idx = PDM_array_new_idx_from_sizes_int(blk_n_cloud_to_surf, _cloud_deform->blk_n_points);

  _cloud_deform->blk_buffer_from_surf_idx = malloc (sizeof(int) * (_cloud_deform->blk_n_points+1));
  _cloud_deform->blk_buffer_from_surf_idx[0] = 0;

  for (int i_pts = 0; i_pts < _cloud_deform->blk_n_points; i_pts++) {

    _cloud_deform->blk_buffer_from_surf_idx[i_pts+1] = _cloud_deform->blk_buffer_from_surf_idx[i_pts];

    if (blk_cloud_to_surf[blk_cloud_to_surf_idx[i_pts]] < 0) {
      int i_layer = -(blk_cloud_to_surf[blk_cloud_to_surf_idx[i_pts]]+1);
      _cloud_deform->blk_buffer_from_surf_idx[i_pts+1] += distri_layer[i_layer][n_rank];
    } else {
      _cloud_deform->blk_buffer_from_surf_idx[i_pts+1] += blk_cloud_to_surf_idx[i_pts+1] - blk_cloud_to_surf_idx[i_pts];
    }

  }

  _cloud_deform->blk_buffer_from_surf = malloc (sizeof(int) * _cloud_deform->blk_buffer_from_surf_idx[_cloud_deform->blk_n_points]);

  int         *blk_cloud_to_surf_buffer_idx = malloc (sizeof(int        ) * (_cloud_deform->blk_n_points+1));
  PDM_g_num_t *blk_cloud_to_surf_buffer     = malloc (sizeof(PDM_g_num_t) * _cloud_deform->blk_buffer_from_surf_idx[_cloud_deform->blk_n_points]);
  blk_cloud_to_surf_buffer_idx[0] = 0;

  PDM_g_num_t *sorted_gnum = malloc (sizeof(PDM_g_num_t) * _cloud_deform->blk_buffer_from_surf_idx[_cloud_deform->blk_n_points]);
  int         *order       = malloc (sizeof(int        ) * _cloud_deform->blk_buffer_from_surf_idx[_cloud_deform->blk_n_points]);

  for (int i_pts = 0; i_pts < _cloud_deform->blk_n_points; i_pts++) {

    blk_cloud_to_surf_buffer_idx[i_pts+1] = blk_cloud_to_surf_buffer_idx[i_pts];

    int idx_write = _cloud_deform->blk_buffer_from_surf_idx[i_pts];

    if (blk_cloud_to_surf[blk_cloud_to_surf_idx[i_pts]] < 0) {

      int i_layer = -(blk_cloud_to_surf[blk_cloud_to_surf_idx[i_pts]]+1);

      for (int k_pts = 0; k_pts < blk_cloud_to_surf_buffer_idx[i_pts+1]; k_pts++) {
        sorted_gnum[k_pts] = blk_cloud_to_surf_buffer[k_pts];
        order      [k_pts] = k_pts;
      }
      PDM_sort_long(sorted_gnum, order, blk_cloud_to_surf_buffer_idx[i_pts+1]);

      int i_gnum = -1;
      if (blk_cloud_to_surf_buffer_idx[i_pts+1] > 0) {
        i_gnum = PDM_binary_search_long(gnum_layer[i_layer], sorted_gnum, blk_cloud_to_surf_buffer_idx[i_pts+1]);
      }

      int idx_write2 = 0;
      if (i_gnum < 0) {
        for (int j_pts = 0; j_pts < distri_layer[i_layer][n_rank]; j_pts++) {
          idx_write2 = blk_cloud_to_surf_buffer_idx[i_pts+1];
          blk_cloud_to_surf_buffer[idx_write2] = gnum_layer[i_layer]+j_pts;
          blk_cloud_to_surf_buffer_idx[i_pts+1]++;
          _cloud_deform->blk_buffer_from_surf[idx_write++] = idx_write2;
        }
      } else {
        for (int j_pts = 0; j_pts < distri_layer[i_layer][n_rank]; j_pts++) {
          _cloud_deform->blk_buffer_from_surf[idx_write++] = order[i_gnum]+j_pts;
        }
      }

    } else {

      for (int j_pts = blk_cloud_to_surf_idx[i_pts]; j_pts < blk_cloud_to_surf_idx[i_pts+1]; j_pts++) {

        for (int k_pts = 0; k_pts < blk_cloud_to_surf_buffer_idx[i_pts+1]; k_pts++) {
          sorted_gnum[k_pts] = blk_cloud_to_surf_buffer[k_pts];
          order      [k_pts] = k_pts;
        }
        PDM_sort_long(sorted_gnum, order, blk_cloud_to_surf_buffer_idx[i_pts+1]);

        int i_gnum = -1;
        if (blk_cloud_to_surf_buffer_idx[i_pts+1] > 0) {
          i_gnum = PDM_binary_search_long(blk_cloud_to_surf[j_pts], sorted_gnum, blk_cloud_to_surf_buffer_idx[i_pts+1]);
        }

        int idx_write2 = 0;
        if (i_gnum < 0) {
          idx_write2 = blk_cloud_to_surf_buffer_idx[i_pts+1];
          blk_cloud_to_surf_buffer[idx_write2] = blk_cloud_to_surf[j_pts];
          blk_cloud_to_surf_buffer_idx[i_pts+1]++;
          _cloud_deform->blk_buffer_from_surf[idx_write++] = idx_write2;
        } else {
          _cloud_deform->blk_buffer_from_surf[idx_write++] = order[i_gnum];
        }

      }

    }

  }

  free(sorted_gnum);
  free(order);

  blk_cloud_to_surf_buffer = realloc (blk_cloud_to_surf_buffer, blk_cloud_to_surf_buffer_idx[_cloud_deform->blk_n_points] * sizeof(PDM_g_num_t));

  printf("Reduction of data size : ini = %i ; end = %i\n", _cloud_deform->blk_buffer_from_surf_idx[_cloud_deform->blk_n_points],
                                                           blk_cloud_to_surf_buffer_idx[_cloud_deform->blk_n_points]);

  PDM_g_num_t *blk_gnum = PDM_part_to_block_block_gnum_get(_cloud_deform->ptb_blk);

  _cloud_deform->ptp_blk = PDM_part_to_part_create((const PDM_g_num_t **) &blk_gnum,
                                                   (const int          *) &_cloud_deform->blk_n_points,
                                                                           1,
                                                   (const PDM_g_num_t **)  _surf_deform->vtx_ln_to_gn,
                                                   (const int          *)  _surf_deform->n_vtx,
                                                                           _surf_deform->n_part+def->n_layer,
                                                   (const int         **) &blk_cloud_to_surf_buffer_idx,
                                                   (const PDM_g_num_t **) &blk_cloud_to_surf_buffer,
                                                                           def->comm);

  for (int i_part = 0; i_part < _cloud_deform->n_part; i_part++) {
    free(cloud_weight     [i_part]);
    free(n_cloud_to_surf  [i_part]);
    free(cloud_to_surf_idx[i_part]);
    free(cloud_to_surf    [i_part]);
  }

  free(cloud_weight     );
  free(n_cloud_to_surf  );
  free(cloud_to_surf_idx);
  free(cloud_to_surf    );
  free(blk_n_cloud_to_surf  );
  free(blk_cloud_to_surf_idx);
  free(blk_cloud_to_surf    );
  free(blk_cloud_to_surf_buffer_idx);
  free(blk_cloud_to_surf_buffer    );

  for (int i_layer = 0; i_layer < def->n_layer; i_layer++) {
    free(distri_layer[i_layer]);
  }

  free(distri_layer);
  free(gnum_layer  );

  def->is_cpt = 1;

}

static void
_PDM_mesh_deform_cloud_block_get
(
  PDM_mesh_deform_t  *def,
  int               **blk_buffer_from_surf_idx,
  int               **blk_buffer_from_surf,
  double            **blk_coords_from_surf,
  double            **blk_dcoords_from_surf,
  double            **blk_aux_geom_from_surf,
  int                *blk_cloud_n_points,
  double            **blk_cloud_coords,
  double            **blk_cloud_dcoords
)
{

  assert(def->is_cpt == 1);

  PDM_surf_deform_t  *_surf_deform  = def->surf_deform;
  PDM_cloud_deform_t *_cloud_deform = def->cloud_deform;

  if (def->is_blk == 1) {
    free(_cloud_deform->blk_coords );
    free(_cloud_deform->blk_dcoords);
    free(_cloud_deform->blk_coords_from_surf[0] );
    free(_cloud_deform->blk_dcoords_from_surf[0]);
    free(_cloud_deform->blk_coords_from_surf    );
    free(_cloud_deform->blk_dcoords_from_surf   );
    if (def->n_aux_geom > 0) {
      free(_cloud_deform->blk_aux_geom_from_surf[0]);
      free(_cloud_deform->blk_aux_geom_from_surf   );
    }
    def->is_blk = 0;
  }

  for (int i_layer = 0; i_layer < def->n_layer; i_layer++) {

    int  n_explicit_nodes = 0;
    int *n_points         = NULL;
    int *range            = NULL;
    int *leaf_id          = NULL;
    int *leaf_nodes_id    = NULL;
    int *children_id      = NULL;
    int *ancestor_id      = NULL;
    int  n_child          = 0;
    int  n_leaf           = 0;
    int  stack_size       = 0;

    PDM_para_octree_explicit_node_get(_surf_deform->octree_layer[i_layer],
                                     &n_explicit_nodes,
                                     &n_points,
                                     &range,
                                     &leaf_id,
                                     &children_id,
                                     &ancestor_id,
                                     &n_child,
                                     &stack_size);

    _leaf_get(n_explicit_nodes,
              leaf_id,
             &n_leaf,
             &leaf_nodes_id);

    int      request    = 0;
    double **tmp_octree = NULL;

    PDM_part_to_part_reverse_iexch(_surf_deform->ptp_layer[i_layer],
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   3,
                                   sizeof(double),
                                   NULL,
                 (const void **)   _surf_deform->coords,
                                   NULL,
                       (void ***) &tmp_octree,
                                  &request);

    PDM_part_to_part_reverse_iexch_wait(_surf_deform->ptp_layer[i_layer], request);

    _mean_array_double_per_leaf(n_leaf,
                                leaf_nodes_id,
                                n_points,
                                range,
                                3,
                                tmp_octree[0],
                                _surf_deform->coords[_surf_deform->n_part+i_layer]);

    free(tmp_octree[0]);
    free(tmp_octree   );

    PDM_part_to_part_reverse_iexch(_surf_deform->ptp_layer[i_layer],
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   3,
                                   sizeof(double),
                                   NULL,
                 (const void **)   _surf_deform->dcoords,
                                   NULL,
                       (void ***) &tmp_octree,
                                  &request);

    PDM_part_to_part_reverse_iexch_wait(_surf_deform->ptp_layer[i_layer], request);

    _mean_array_double_per_leaf(n_leaf,
                                leaf_nodes_id,
                                n_points,
                                range,
                                3,
                                tmp_octree[0],
                                _surf_deform->dcoords[_surf_deform->n_part+i_layer]);

    free(tmp_octree[0]);
    free(tmp_octree   );

    if (def->n_aux_geom > 0) {

      PDM_part_to_part_reverse_iexch(_surf_deform->ptp_layer[i_layer],
                                     PDM_MPI_COMM_KIND_P2P,
                                     PDM_STRIDE_CST_INTERLACED,
                                     PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                     def->n_aux_geom,
                                     sizeof(double),
                                     NULL,
                   (const void **)   _surf_deform->aux_geom,
                                     NULL,
                         (void ***) &tmp_octree,
                                    &request);

      PDM_part_to_part_reverse_iexch_wait(_surf_deform->ptp_layer[i_layer], request);

      _mean_array_double_per_leaf(n_leaf,
                                  leaf_nodes_id,
                                  n_points,
                                  range,
                                  def->n_aux_geom,
                                  tmp_octree[0],
                                  _surf_deform->aux_geom[_surf_deform->n_part+i_layer]);

      free(tmp_octree[0]);
      free(tmp_octree   );

    }

    free(leaf_nodes_id);

  }

  int request = 0;

  PDM_part_to_part_reverse_iexch(_cloud_deform->ptp_blk,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 3,
                                 sizeof(double),
                                 NULL,
               (const void **)   _surf_deform->coords,
                                 NULL,
                     (void ***) &_cloud_deform->blk_coords_from_surf,
                                &request);

  PDM_part_to_part_reverse_iexch_wait(_cloud_deform->ptp_blk, request);

  PDM_part_to_part_reverse_iexch(_cloud_deform->ptp_blk,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 3,
                                 sizeof(double),
                                 NULL,
               (const void **)   _surf_deform->dcoords,
                                 NULL,
                     (void ***) &_cloud_deform->blk_dcoords_from_surf,
                                &request);

  PDM_part_to_part_reverse_iexch_wait(_cloud_deform->ptp_blk, request);

  if (def->n_aux_geom > 0) {

    PDM_part_to_part_reverse_iexch(_cloud_deform->ptp_blk,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   def->n_aux_geom,
                                   sizeof(double),
                                   NULL,
                 (const void **)   _surf_deform->aux_geom,
                                   NULL,
                       (void ***) &_cloud_deform->blk_aux_geom_from_surf,
                                  &request);

    PDM_part_to_part_reverse_iexch_wait(_cloud_deform->ptp_blk, request);

  }

  int blk_size = PDM_part_to_block_exch(_cloud_deform->ptb_blk,
                                        sizeof(double),
                                        PDM_STRIDE_CST_INTERLACED,
                                        3,
                                        NULL,
                             (void **)  _cloud_deform->coords,
                                        NULL,
                             (void **) &_cloud_deform->blk_coords);

  PDM_UNUSED(blk_size);

  _cloud_deform->blk_dcoords = malloc (sizeof(double) * 3 * _cloud_deform->blk_n_points);

  for (int i_pts = 0; i_pts < _cloud_deform->blk_n_points; i_pts++) {
    _cloud_deform->blk_dcoords[3*i_pts    ] = 0.0;
    _cloud_deform->blk_dcoords[3*i_pts + 1] = 0.0;
    _cloud_deform->blk_dcoords[3*i_pts + 2] = 0.0;
  }

  *blk_buffer_from_surf_idx = _cloud_deform->blk_buffer_from_surf_idx;
  *blk_buffer_from_surf     = _cloud_deform->blk_buffer_from_surf;
  *blk_coords_from_surf     = _cloud_deform->blk_coords_from_surf[0];
  *blk_dcoords_from_surf    = _cloud_deform->blk_dcoords_from_surf[0];
  if (def->n_aux_geom > 0) {
    *blk_aux_geom_from_surf = _cloud_deform->blk_aux_geom_from_surf[0];
  } else {
    *blk_aux_geom_from_surf = NULL;
  }

  *blk_cloud_n_points = _cloud_deform->blk_n_points;
  *blk_cloud_coords   = _cloud_deform->blk_coords;
  *blk_cloud_dcoords  = _cloud_deform->blk_dcoords;

  def->is_blk = 1;

}

static void
_PDM_mesh_deform_cloud_dcoords_part_get
(
  PDM_mesh_deform_t   *def,
  double            ***cloud_dcoords
)
{

  assert(def->is_cpt == 1);
  assert(def->is_blk == 1);

  PDM_cloud_deform_t *_cloud_deform = def->cloud_deform;

  if (def->is_res == 1) {
    for (int i_part = 0; i_part < _cloud_deform->n_part; i_part++) {
      free(_cloud_deform->dcoords[i_part]);
    }
    free(_cloud_deform->dcoords);
    def->is_res = 0;
  }

  PDM_part_to_block_reverse_exch(_cloud_deform->ptb_blk,
                                 sizeof(double),
                                 PDM_STRIDE_CST_INTERLACED,
                                 3,
                                 NULL,
                     (void   *)  _cloud_deform->blk_dcoords,
                                 NULL,
                     (void ***) &_cloud_deform->dcoords);

  *cloud_dcoords = _cloud_deform->dcoords;

  def->is_res = 1;

}

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init(&argc, &argv);

  /*
   * Parameters
   */

  PDM_bool_t  show_mesh        = PDM_FALSE;
  PDM_g_num_t n_vtx_seg        = 30;
  double      length           = 1e-2;
  double      zero_x           = -length/2;
  double      zero_y           = -length/2;
  double      zero_z           = -length/2;
  int         n_part           = 2;
  int         n_part_bnd       = 1;
  double      min_dist         = pow(length/20, 2);
  int         n_layer          = 4;
  int        *n_leaf_per_layer = malloc (n_layer*sizeof(int));
  int         n_vtx_min_dist   = 10;
  int         n_var            = 1;

  n_leaf_per_layer[0] = 31;
  n_leaf_per_layer[1] = 62;
  n_leaf_per_layer[2] = 125;
  n_leaf_per_layer[3] = 250;

  /*
   * Generate a cube
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int          dn_face_group;
  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  int          dsface_vtx;
  int          dsface_group;
  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  double      *dvtx_coord      = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  if (i_rank == 0) {
    printf("-- Generate cube\n");
    fflush(stdout);
  }

  PDM_dcube_t *dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          zero_x,
                                          zero_y,
                                          zero_z,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                       &dn_face_group,
                       &dn_cell,
                       &dn_face,
                       &dn_vtx,
                       &dsface_vtx,
                       &dsface_group);

  PDM_dcube_gen_data_get(dcube,
                        &dface_cell,
                        &dface_vtx_idx,
                        &dface_vtx,
                        &dvtx_coord,
                        &dface_group_idx,
                        &dface_group);

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  if (i_rank == 0) {
    printf("-- Part cube\n");
    fflush(stdout);
  }

  PDM_part_t* ppart = PDM_part_create(comm,
                                      PDM_PART_SPLIT_HILBERT,
                                      "PDM_PART_RENUM_CELL_NONE",
                                      "PDM_PART_RENUM_FACE_NONE",
                                      n_property_cell,
                                      renum_properties_cell,
                                      n_property_face,
                                      renum_properties_face,
                                      n_part,
                                      dn_cell,
                                      dn_face,
                                      dn_vtx,
                                      dn_face_group,
                                      NULL,
                                      NULL,
                                      NULL,
                                      NULL,
                                      have_dcell_part,
                                      dcell_part,
                                      dface_cell,
                                      dface_vtx_idx,
                                      dface_vtx,
                                      NULL,
                                      dvtx_coord,
                                      NULL,
                                      dface_group_idx,
                                      dface_group);

  free(dcell_part);
  PDM_dcube_gen_free(dcube);

  int **selected_face = malloc (sizeof(int *) * n_part);

  PDM_g_num_t **bnd_face_ln_to_gn = malloc (sizeof(PDM_g_num_t *) * n_part_bnd);
  int         **bnd_face_vtx_idx  = malloc (sizeof(int         *) * n_part_bnd);
  int         **bnd_face_vtx      = malloc (sizeof(int         *) * n_part_bnd);
  int          *n_bnd_face        = malloc (sizeof(int          ) * n_part_bnd);

  PDM_g_num_t **bnd_vtx_ln_to_gn = malloc (sizeof(PDM_g_num_t *) * n_part_bnd);
  double      **bnd_vtx          = malloc (sizeof(double      *) * n_part_bnd);
  double      **bnd_dvtx         = malloc (sizeof(double      *) * n_part_bnd);
  double      **bnd_vtx_aux_geom = malloc (sizeof(double      *) * n_part_bnd);
  int          *n_bnd_vtx        = malloc (sizeof(int          ) * n_part_bnd);

  PDM_g_num_t **int_vtx_gnum = malloc (sizeof(PDM_g_num_t *) * n_part);
  double      **int_vtx      = malloc (sizeof(double      *) * n_part);
  int          *n_int_vtx    = malloc (sizeof(int          ) * n_part);

  PDM_extract_part_t* extrp = PDM_extract_part_create(2,
                                                      n_part,
                                                      n_part_bnd,
                                                      PDM_EXTRACT_PART_KIND_REEQUILIBRATE,
                                                      PDM_SPLIT_DUAL_WITH_HILBERT,
                                                      PDM_TRUE,
                                                      PDM_OWNERSHIP_KEEP,
                                                      comm);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nface_group;

    PDM_part_part_dim_get(ppart,
                          i_part,
                         &n_cell,
                         &n_face,
                         &n_face_part_bound,
                         &n_vtx,
                         &n_proc,
                         &n_total_part,
                         &scell_face,
                         &sface_vtx,
                         &sface_group,
                         &nface_group);

    int         *cell_tag                 = NULL;
    int         *cell_face_idx            = NULL;
    int         *cell_face                = NULL;
    PDM_g_num_t *cell_ln_to_gn            = NULL;
    int         *face_tag                 = NULL;
    int         *face_cell                = NULL;
    int         *face_vtx_idx             = NULL;
    int         *face_vtx                 = NULL;
    PDM_g_num_t *face_ln_to_gn            = NULL;
    int         *face_part_bound_proc_idx = NULL;
    int         *face_part_bound_part_idx = NULL;
    int         *face_part_bound          = NULL;
    int         *vtx_tag                  = NULL;
    double      *vtx                      = NULL;
    PDM_g_num_t *vtx_ln_to_gn             = NULL;
    int         *face_group_idx           = NULL;
    int         *face_group               = NULL;
    PDM_g_num_t *face_group_ln_to_gn      = NULL;

    PDM_part_part_val_get(ppart,
                          i_part,
                         &cell_tag,
                         &cell_face_idx,
                         &cell_face,
                         &cell_ln_to_gn,
                         &face_tag,
                         &face_cell,
                         &face_vtx_idx,
                         &face_vtx,
                         &face_ln_to_gn,
                         &face_part_bound_proc_idx,
                         &face_part_bound_part_idx,
                         &face_part_bound,
                         &vtx_tag,
                         &vtx,
                         &vtx_ln_to_gn,
                         &face_group_idx,
                         &face_group,
                         &face_group_ln_to_gn);

    if (show_mesh) {
      char filename[999];
      sprintf(filename, "ini_face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_polydata(filename,
                             n_vtx,
                             vtx,
                             vtx_ln_to_gn,
                             n_face,
                             face_vtx_idx,
                             face_vtx,
                             face_ln_to_gn,
                             NULL);
    }

    selected_face[i_part] = malloc (sizeof(int) * sface_group);

    for (int i_group = 0; i_group < nface_group; i_group++) {
      for (int i_face_group = face_group_idx[i_group]; i_face_group < face_group_idx[i_group+1]; i_face_group++) {
        selected_face[i_part][i_face_group] = face_group[i_face_group]-1;
      }
    }

    PDM_extract_part_part_set(extrp,
                              i_part,
                              n_cell,
                              n_face,
                              -1,
                              n_vtx,
                              cell_face_idx,
                              cell_face,
                              NULL,
                              NULL,
                              NULL,
                              face_vtx_idx,
                              face_vtx,
                              cell_ln_to_gn,
                              face_ln_to_gn,
                              NULL,
                              vtx_ln_to_gn,
                              vtx);

    PDM_extract_part_selected_lnum_set(extrp,
                                       i_part,
                                       sface_group,
                                       selected_face[i_part]);

  }

  PDM_extract_part_compute(extrp);

  for (int i_part = 0; i_part < n_part_bnd; i_part++) {

    /*
     * Get boundary faces and vertices
     */

    n_bnd_face[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                       i_part,
                                                       PDM_MESH_ENTITY_FACE);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                 &bnd_face_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_USER);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                     &bnd_face_vtx    [i_part],
                                     &bnd_face_vtx_idx[i_part],
                                      PDM_OWNERSHIP_USER);

    n_bnd_vtx[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                      i_part,
                                                      PDM_MESH_ENTITY_VERTEX);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_VERTEX,
                                 &bnd_vtx_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_USER);

    PDM_extract_part_vtx_coord_get(extrp,
                                   i_part,
                                  &bnd_vtx[i_part],
                                   PDM_OWNERSHIP_USER);

    bnd_dvtx        [i_part] = malloc (sizeof(double) * 3     * n_bnd_vtx[i_part]);
    bnd_vtx_aux_geom[i_part] = malloc (sizeof(double) * n_var * n_bnd_vtx[i_part]);

    for (int i_vtx = 0; i_vtx < n_bnd_vtx[i_part]; i_vtx++) {
      for (int i_dir = 0; i_dir < 3; i_dir++) {
        bnd_dvtx[i_part][3*i_vtx+i_dir] = 0.0;
      }
      for (int i_var = 0; i_var < n_var; i_var++) {
        bnd_vtx_aux_geom[i_part][n_var*i_vtx+i_var] = 0.0;
      }
      bnd_dvtx        [i_part][3*i_vtx] = 2e-3*(1.0-(bnd_vtx[i_part][3*i_vtx]+length/2)/length);
      bnd_vtx_aux_geom[i_part][  i_vtx] = 1.0;
    }

  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(selected_face[i_part]);
  }

  free(selected_face);

  PDM_part_to_part_t *ptp_bnd_vtx = NULL;

  PDM_extract_part_part_to_part_get(extrp,
                                    PDM_MESH_ENTITY_VERTEX,
                                   &ptp_bnd_vtx,
                                    PDM_OWNERSHIP_USER);

  PDM_extract_part_free(extrp);

  PDM_part_to_block_t *ptb_bnd_vtx = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                              1.,
                                                              bnd_vtx_ln_to_gn,
                                                              NULL,
                                                              n_bnd_vtx,
                                                              n_part_bnd,
                                                              comm);

  int  *n_ref_vtx = NULL;
  int **ref_vtx   = NULL;

  PDM_part_to_part_ref_lnum2_get(ptp_bnd_vtx,
                                &n_ref_vtx,
                                &ref_vtx);

  int  *n_unref_vtx = NULL;
  int **unref_vtx   = NULL;

  PDM_part_to_part_unref_lnum2_get(ptp_bnd_vtx,
                                  &n_unref_vtx,
                                  &unref_vtx);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nface_group;

    PDM_part_part_dim_get(ppart,
                          i_part,
                         &n_cell,
                         &n_face,
                         &n_face_part_bound,
                         &n_vtx,
                         &n_proc,
                         &n_total_part,
                         &scell_face,
                         &sface_vtx,
                         &sface_group,
                         &nface_group);

    int         *cell_tag                 = NULL;
    int         *cell_face_idx            = NULL;
    int         *cell_face                = NULL;
    PDM_g_num_t *cell_ln_to_gn            = NULL;
    int         *face_tag                 = NULL;
    int         *face_cell                = NULL;
    int         *face_vtx_idx             = NULL;
    int         *face_vtx                 = NULL;
    PDM_g_num_t *face_ln_to_gn            = NULL;
    int         *face_part_bound_proc_idx = NULL;
    int         *face_part_bound_part_idx = NULL;
    int         *face_part_bound          = NULL;
    int         *vtx_tag                  = NULL;
    double      *vtx                      = NULL;
    PDM_g_num_t *vtx_ln_to_gn             = NULL;
    int         *face_group_idx           = NULL;
    int         *face_group               = NULL;
    PDM_g_num_t *face_group_ln_to_gn      = NULL;

    PDM_part_part_val_get(ppart,
                          i_part,
                         &cell_tag,
                         &cell_face_idx,
                         &cell_face,
                         &cell_ln_to_gn,
                         &face_tag,
                         &face_cell,
                         &face_vtx_idx,
                         &face_vtx,
                         &face_ln_to_gn,
                         &face_part_bound_proc_idx,
                         &face_part_bound_part_idx,
                         &face_part_bound,
                         &vtx_tag,
                         &vtx,
                         &vtx_ln_to_gn,
                         &face_group_idx,
                         &face_group,
                         &face_group_ln_to_gn);

   /*
    * Get interior vertices
    */

    n_int_vtx[i_part] = n_unref_vtx[i_part];

    int_vtx_gnum[i_part] = malloc (sizeof(PDM_g_num_t)     * n_int_vtx[i_part]);
    int_vtx     [i_part] = malloc (sizeof(double     ) * 3 * n_int_vtx[i_part]);

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {
      int_vtx_gnum[i_part][  i_vtx    ] = vtx_ln_to_gn[unref_vtx[i_part][i_vtx]-1];
      int_vtx     [i_part][3*i_vtx    ] = vtx[3*(unref_vtx[i_part][i_vtx]-1)    ];
      int_vtx     [i_part][3*i_vtx + 1] = vtx[3*(unref_vtx[i_part][i_vtx]-1) + 1];
      int_vtx     [i_part][3*i_vtx + 2] = vtx[3*(unref_vtx[i_part][i_vtx]-1) + 2];
    }

  }

  _compute_gnum_in_place(n_part,
                         n_int_vtx,
                         int_vtx_gnum,
                         comm);

  /*
   * Create PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Create PDM_mesh_deform\n");
    fflush(stdout);
  }

  PDM_mesh_deform_t *def = _PDM_mesh_deform_create(n_part_bnd,
                                                   n_part,
                                                   n_var,
                                                   n_layer,
                                                   n_leaf_per_layer,
                                                   min_dist,
                                                   n_vtx_min_dist,
                                                   comm);

  /*
   * Set PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Set PDM_mesh_deform\n");
    fflush(stdout);
  }

  for (int i_part = 0; i_part < n_part_bnd; i_part++) {

    _PDM_mesh_deform_surf_part_set(def,
                                   i_part,
                                   n_bnd_face       [i_part],
                                   n_bnd_vtx        [i_part],
                                   bnd_face_vtx_idx [i_part],
                                   bnd_face_vtx     [i_part],
                                   bnd_face_ln_to_gn[i_part],
                                   bnd_vtx_ln_to_gn [i_part],
                                   bnd_vtx          [i_part],
                                   bnd_dvtx         [i_part],
                                   bnd_vtx_aux_geom [i_part]);

  }

  for (int i_part = 0; i_part < n_part; i_part++) {

    _PDM_mesh_deform_cloud_part_set(def,
                                    i_part,
                                    n_int_vtx   [i_part],
                                    int_vtx_gnum[i_part],
                                    int_vtx     [i_part]);

  }

  /*
   * Compute PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Compute PDM_mesh_deform\n");
    fflush(stdout);
  }

  _PDM_mesh_deform_compute(def);

  /*
   * Get block PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Get block PDM_mesh_deform 1/4\n");
    fflush(stdout);
  }

  int    *blk_buffer_from_bnd_idx = NULL;
  int    *blk_buffer_from_bnd     = NULL;
  double *blk_vtx_from_bnd        = NULL;
  double *blk_dvtx_from_bnd       = NULL;
  double *blk_aux_geom_from_bnd   = NULL;
  int     blk_n_int_vtx           = 0;
  double *blk_int_vtx             = NULL;
  double *blk_int_dvtx            = NULL;

  _PDM_mesh_deform_cloud_block_get(def,
                                  &blk_buffer_from_bnd_idx,
                                  &blk_buffer_from_bnd,
                                  &blk_vtx_from_bnd,
                                  &blk_dvtx_from_bnd,
                                  &blk_aux_geom_from_bnd,
                                  &blk_n_int_vtx,
                                  &blk_int_vtx,
                                  &blk_int_dvtx);

  /*
   * Get part PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Get part PDM_mesh_deform 1/4\n");
    fflush(stdout);
  }

  double **int_dvtx = NULL;

  _PDM_mesh_deform_cloud_dcoords_part_get(def,
                                         &int_dvtx);

  /*
   * Get block PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Get block PDM_mesh_deform 2/4\n");
    fflush(stdout);
  }

  blk_buffer_from_bnd_idx = NULL;
  blk_buffer_from_bnd     = NULL;
  blk_vtx_from_bnd        = NULL;
  blk_dvtx_from_bnd       = NULL;
  blk_aux_geom_from_bnd   = NULL;
  blk_n_int_vtx           = 0;
  blk_int_vtx             = NULL;
  blk_int_dvtx            = NULL;

  _PDM_mesh_deform_cloud_block_get(def,
                                  &blk_buffer_from_bnd_idx,
                                  &blk_buffer_from_bnd,
                                  &blk_vtx_from_bnd,
                                  &blk_dvtx_from_bnd,
                                  &blk_aux_geom_from_bnd,
                                  &blk_n_int_vtx,
                                  &blk_int_vtx,
                                  &blk_int_dvtx);

  /*
   * Get part PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Get part PDM_mesh_deform 2/4\n");
    fflush(stdout);
  }

  int_dvtx = NULL;

  _PDM_mesh_deform_cloud_dcoords_part_get(def,
                                         &int_dvtx);

  /*
   * Recompute PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Recompute PDM_mesh_deform\n");
    fflush(stdout);
  }

  _PDM_mesh_deform_compute(def);

  /*
   * Get block PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Get block PDM_mesh_deform 3/4\n");
    fflush(stdout);
  }

  blk_buffer_from_bnd_idx = NULL;
  blk_buffer_from_bnd     = NULL;
  blk_vtx_from_bnd        = NULL;
  blk_dvtx_from_bnd       = NULL;
  blk_aux_geom_from_bnd   = NULL;
  blk_n_int_vtx           = 0;
  blk_int_vtx             = NULL;
  blk_int_dvtx            = NULL;

  _PDM_mesh_deform_cloud_block_get(def,
                                  &blk_buffer_from_bnd_idx,
                                  &blk_buffer_from_bnd,
                                  &blk_vtx_from_bnd,
                                  &blk_dvtx_from_bnd,
                                  &blk_aux_geom_from_bnd,
                                  &blk_n_int_vtx,
                                  &blk_int_vtx,
                                  &blk_int_dvtx);

  /*
   * Get part PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Get part PDM_mesh_deform 3/4\n");
    fflush(stdout);
  }

  int_dvtx = NULL;

  _PDM_mesh_deform_cloud_dcoords_part_get(def,
                                         &int_dvtx);

  /*
   * Get block PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Get block PDM_mesh_deform 4/4\n");
    fflush(stdout);
  }

  blk_buffer_from_bnd_idx = NULL;
  blk_buffer_from_bnd     = NULL;
  blk_vtx_from_bnd        = NULL;
  blk_dvtx_from_bnd       = NULL;
  blk_aux_geom_from_bnd   = NULL;
  blk_n_int_vtx           = 0;
  blk_int_vtx             = NULL;
  blk_int_dvtx            = NULL;

  _PDM_mesh_deform_cloud_block_get(def,
                                  &blk_buffer_from_bnd_idx,
                                  &blk_buffer_from_bnd,
                                  &blk_vtx_from_bnd,
                                  &blk_dvtx_from_bnd,
                                  &blk_aux_geom_from_bnd,
                                  &blk_n_int_vtx,
                                  &blk_int_vtx,
                                  &blk_int_dvtx);

  /*
   * Move vertices
   */

  if (i_rank == 0) {
    printf("-- Move vertices\n");
    fflush(stdout);
  }

  int     blk_n_bnd_vtx        = PDM_part_to_block_n_elt_block_get(ptb_bnd_vtx);
  double *blk_bnd_vtx          = NULL;
  double *blk_bnd_dvtx         = NULL;
  double *blk_bnd_vtx_aux_geom = NULL;

  int blk_size = PDM_part_to_block_exch(ptb_bnd_vtx,
                                        sizeof(double),
                                        PDM_STRIDE_CST_INTERLACED,
                                        3,
                                        NULL,
                             (void **)  bnd_vtx,
                                        NULL,
                             (void **) &blk_bnd_vtx);

  PDM_UNUSED(blk_size);

  blk_size = PDM_part_to_block_exch(ptb_bnd_vtx,
                                    sizeof(double),
                                    PDM_STRIDE_CST_INTERLACED,
                                    3,
                                    NULL,
                         (void **)  bnd_dvtx,
                                    NULL,
                         (void **) &blk_bnd_dvtx);

  blk_size = PDM_part_to_block_exch(ptb_bnd_vtx,
                                    sizeof(double),
                                    PDM_STRIDE_CST_INTERLACED,
                                    n_var,
                                    NULL,
                         (void **)  bnd_vtx_aux_geom,
                                    NULL,
                         (void **) &blk_bnd_vtx_aux_geom);

  double *dx = malloc (sizeof(double) * 3);
  double *du = malloc (sizeof(double) * 3);
  double *dr = malloc (sizeof(double) * 3);
  double  sdist;
  double  dist;

  PDM_g_num_t  m_vtx = 0;
  PDM_g_num_t _m_vtx = 0;
  double       aire  = 0.0;
  double      _aire  = 0.0;
  double       l1    = 0.0;
  double      _l1    = 0.0;
  double       l2    = 0.0;
  double      _l2    = 0.0;

  dx[0] = 0.0;
  dx[1] = 0.0;
  dx[2] = 0.0;

  for (int i_vtx = 0; i_vtx < blk_n_bnd_vtx; i_vtx++) {
    m_vtx++;
    dx[0] = dx[0] + blk_bnd_vtx[3*i_vtx    ];
    dx[1] = dx[1] + blk_bnd_vtx[3*i_vtx + 1];
    dx[2] = dx[2] + blk_bnd_vtx[3*i_vtx + 2];
  }

  PDM_MPI_Allreduce (&m_vtx, &_m_vtx, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce (dx, dr, 3, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  dx[0] = dr[0]/_m_vtx;
  dx[1] = dr[1]/_m_vtx;
  dx[2] = dr[2]/_m_vtx;

  for (int i_vtx = 0; i_vtx < blk_n_bnd_vtx; i_vtx++) {
    dr[0] = dx[0] - blk_bnd_vtx[3*i_vtx    ];
    dr[1] = dx[1] - blk_bnd_vtx[3*i_vtx + 1];
    dr[2] = dx[2] - blk_bnd_vtx[3*i_vtx + 2];
    l1 = PDM_MAX(l1, pow(dr[0],2)+pow(dr[1],2)+pow(dr[2],2));
  }

  PDM_MPI_Allreduce (&l1, &_l1, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
  l1 = sqrt(_l1);

  dx[0] = 0.0;
  dx[1] = 0.0;
  dx[2] = 0.0;

  for (int i_vtx = 0; i_vtx < blk_n_bnd_vtx; i_vtx++) {
    aire  = aire  + blk_bnd_vtx_aux_geom[n_var*i_vtx];
    dx[0] = dx[0] + blk_bnd_vtx_aux_geom[n_var*i_vtx]*blk_bnd_dvtx[3*i_vtx    ];
    dx[1] = dx[1] + blk_bnd_vtx_aux_geom[n_var*i_vtx]*blk_bnd_dvtx[3*i_vtx + 1];
    dx[2] = dx[2] + blk_bnd_vtx_aux_geom[n_var*i_vtx]*blk_bnd_dvtx[3*i_vtx + 2];
  }

  PDM_MPI_Allreduce (dx, dr, 3, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce (&aire, &_aire, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  dx[0] = dr[0]/_aire;
  dx[1] = dr[1]/_aire;
  dx[2] = dr[2]/_aire;

  for (int i_vtx = 0; i_vtx < blk_n_bnd_vtx; i_vtx++) {
    dr[0] = dx[0] - blk_bnd_dvtx[3*i_vtx    ];
    dr[1] = dx[1] - blk_bnd_dvtx[3*i_vtx + 1];
    dr[2] = dx[2] - blk_bnd_dvtx[3*i_vtx + 2];
    l2 = PDM_MAX(l2, pow(dr[0],2)+pow(dr[1],2)+pow(dr[2],2));
  }

  PDM_MPI_Allreduce (&l2, &_l2, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);
  l2 = PDM_MAX(0.1, 5.0*sqrt(_l2)/l1)*l1;

  free(blk_bnd_vtx );
  free(blk_bnd_dvtx);
  free(blk_bnd_vtx_aux_geom);

  for (int i_vtx = 0; i_vtx < blk_n_int_vtx; i_vtx++) {
    dx[0] = 0.0;
    dx[1] = 0.0;
    dx[2] = 0.0;
    sdist = 0.0;
    for (int j_vtx = blk_buffer_from_bnd_idx[i_vtx]; j_vtx < blk_buffer_from_bnd_idx[i_vtx+1]; j_vtx++) {
      int k_vtx = blk_buffer_from_bnd[j_vtx];
      du[0] = blk_dvtx_from_bnd[3*k_vtx    ];
      du[1] = blk_dvtx_from_bnd[3*k_vtx + 1];
      du[2] = blk_dvtx_from_bnd[3*k_vtx + 2];
      dr[0] = blk_vtx_from_bnd[3*k_vtx    ] - blk_int_vtx[3*i_vtx    ];
      dr[1] = blk_vtx_from_bnd[3*k_vtx + 1] - blk_int_vtx[3*i_vtx + 1];
      dr[2] = blk_vtx_from_bnd[3*k_vtx + 2] - blk_int_vtx[3*i_vtx + 2];
      dist  = sqrt(pow(dr[0],2)+pow(dr[1],2)+pow(dr[2],2));
      if (sqrt(pow(du[0],2)+pow(du[1],2)+pow(du[2],2)) > 1e-6) {
        dist = blk_aux_geom_from_bnd[n_var*k_vtx]*(pow(l1/dist,3) + pow(l2/dist,5));
      } else {
        dist = blk_aux_geom_from_bnd[n_var*k_vtx]*pow(l1/dist,3);
      }
      sdist = sdist + dist;
      dx[0] = dx[0] + blk_dvtx_from_bnd[3*k_vtx    ]*dist;
      dx[1] = dx[1] + blk_dvtx_from_bnd[3*k_vtx + 1]*dist;
      dx[2] = dx[2] + blk_dvtx_from_bnd[3*k_vtx + 2]*dist;
    }
    blk_int_dvtx[3*i_vtx    ] = dx[0]/sdist;
    blk_int_dvtx[3*i_vtx + 1] = dx[1]/sdist;
    blk_int_dvtx[3*i_vtx + 2] = dx[2]/sdist;
  }

  free(dx);
  free(du);
  free(dr);

  /*
   * Get part PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Get part PDM_mesh_deform 4/4\n");
    fflush(stdout);
  }

  int_dvtx = NULL;

  _PDM_mesh_deform_cloud_dcoords_part_get(def,
                                         &int_dvtx);

  for (int i_part = 0; i_part < n_part_bnd; i_part++) {

    for (int i_vtx = 0; i_vtx < n_bnd_vtx[i_part]; i_vtx++) {
      bnd_vtx[i_part][3*i_vtx   ] = bnd_vtx[i_part][3*i_vtx   ] + bnd_dvtx[i_part][3*i_vtx    ];
      bnd_vtx[i_part][3*i_vtx+ 1] = bnd_vtx[i_part][3*i_vtx+ 1] + bnd_dvtx[i_part][3*i_vtx + 1];
      bnd_vtx[i_part][3*i_vtx+ 2] = bnd_vtx[i_part][3*i_vtx+ 2] + bnd_dvtx[i_part][3*i_vtx + 2];
    }

  }

  for (int i_part = 0; i_part < n_part; i_part++) {

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {
      int_vtx[i_part][3*i_vtx   ] = int_vtx[i_part][3*i_vtx   ] + int_dvtx[i_part][3*i_vtx    ];
      int_vtx[i_part][3*i_vtx+ 1] = int_vtx[i_part][3*i_vtx+ 1] + int_dvtx[i_part][3*i_vtx + 1];
      int_vtx[i_part][3*i_vtx+ 2] = int_vtx[i_part][3*i_vtx+ 2] + int_dvtx[i_part][3*i_vtx + 2];
    }

  }

  /*
   * Save vertices new coordinates
   */

  if (i_rank == 0) {
    printf("-- Save vertices new coordinates\n");
    fflush(stdout);
  }

  int      request     = 0;
  double **tmp_bnd_vtx = NULL;

  PDM_part_to_part_iexch(ptp_bnd_vtx,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         3,
                         sizeof(double),
                         NULL,
         (const void **) bnd_vtx,
                         NULL,
             (void ***) &tmp_bnd_vtx,
                        &request);

  PDM_part_to_part_iexch_wait(ptp_bnd_vtx, request);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_total_part;
    int scell_face;
    int sface_vtx;
    int sface_group;
    int nface_group;

    PDM_part_part_dim_get(ppart,
                          i_part,
                         &n_cell,
                         &n_face,
                         &n_face_part_bound,
                         &n_vtx,
                         &n_proc,
                         &n_total_part,
                         &scell_face,
                         &sface_vtx,
                         &sface_group,
                         &nface_group);

    int         *cell_tag                 = NULL;
    int         *cell_face_idx            = NULL;
    int         *cell_face                = NULL;
    PDM_g_num_t *cell_ln_to_gn            = NULL;
    int         *face_tag                 = NULL;
    int         *face_cell                = NULL;
    int         *face_vtx_idx             = NULL;
    int         *face_vtx                 = NULL;
    PDM_g_num_t *face_ln_to_gn            = NULL;
    int         *face_part_bound_proc_idx = NULL;
    int         *face_part_bound_part_idx = NULL;
    int         *face_part_bound          = NULL;
    int         *vtx_tag                  = NULL;
    double      *vtx                      = NULL;
    PDM_g_num_t *vtx_ln_to_gn             = NULL;
    int         *face_group_idx           = NULL;
    int         *face_group               = NULL;
    PDM_g_num_t *face_group_ln_to_gn      = NULL;

    PDM_part_part_val_get(ppart,
                          i_part,
                         &cell_tag,
                         &cell_face_idx,
                         &cell_face,
                         &cell_ln_to_gn,
                         &face_tag,
                         &face_cell,
                         &face_vtx_idx,
                         &face_vtx,
                         &face_ln_to_gn,
                         &face_part_bound_proc_idx,
                         &face_part_bound_part_idx,
                         &face_part_bound,
                         &vtx_tag,
                         &vtx,
                         &vtx_ln_to_gn,
                         &face_group_idx,
                         &face_group,
                         &face_group_ln_to_gn);

    for (int i_vtx = 0; i_vtx < n_ref_vtx[i_part]; i_vtx++) {
      vtx[3*(ref_vtx[i_part][i_vtx]-1)    ] = tmp_bnd_vtx[i_part][3*i_vtx    ];
      vtx[3*(ref_vtx[i_part][i_vtx]-1) + 1] = tmp_bnd_vtx[i_part][3*i_vtx + 1];
      vtx[3*(ref_vtx[i_part][i_vtx]-1) + 2] = tmp_bnd_vtx[i_part][3*i_vtx + 2];
    }

    for (int i_vtx = 0; i_vtx < n_int_vtx[i_part]; i_vtx++) {
      vtx[3*(unref_vtx[i_part][i_vtx]-1)    ] = int_vtx[i_part][3*i_vtx    ];
      vtx[3*(unref_vtx[i_part][i_vtx]-1) + 1] = int_vtx[i_part][3*i_vtx + 1];
      vtx[3*(unref_vtx[i_part][i_vtx]-1) + 2] = int_vtx[i_part][3*i_vtx + 2];
    }

    if (show_mesh) {
      char filename[999];
      sprintf(filename, "end_face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_polydata(filename,
                            n_vtx,
                            vtx,
                            vtx_ln_to_gn,
                            n_face,
                            face_vtx_idx,
                            face_vtx,
                            face_ln_to_gn,
                            NULL);
    }

    free(tmp_bnd_vtx[i_part]);

  }

  free(tmp_bnd_vtx);

  /*
   * Free PDM_mesh_deform
   */

  if (i_rank == 0) {
    printf("-- Free PDM_mesh_deform\n");
    fflush(stdout);
  }

  _PDM_mesh_deform_free(def);

  /*
   * Free memory
   */

  if (i_rank == 0) {
    printf("-- Free memory\n");
    fflush(stdout);
  }

  PDM_part_to_block_free(ptb_bnd_vtx);
  PDM_part_to_part_free(ptp_bnd_vtx);
  PDM_part_free(ppart);

  for (int i_part = 0; i_part < n_part_bnd; i_part++) {
    free(bnd_vtx_ln_to_gn [i_part]);
    free(bnd_vtx          [i_part]);
    free(bnd_dvtx         [i_part]);
    free(bnd_vtx_aux_geom [i_part]);
    free(bnd_face_ln_to_gn[i_part]);
    free(bnd_face_vtx_idx [i_part]);
    free(bnd_face_vtx     [i_part]);
  }
  for (int i_part = 0; i_part < n_part; i_part++) {
    free(int_vtx_gnum[i_part]);
    free(int_vtx     [i_part]);
  }
  free(bnd_face_ln_to_gn);
  free(bnd_face_vtx_idx );
  free(bnd_face_vtx     );
  free(n_bnd_face       );
  free(bnd_vtx_ln_to_gn );
  free(bnd_vtx          );
  free(bnd_dvtx         );
  free(bnd_vtx_aux_geom );
  free(n_bnd_vtx        );
  free(int_vtx_gnum     );
  free(int_vtx          );
  free(n_int_vtx        );
  free(n_leaf_per_layer );

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
