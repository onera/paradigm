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

typedef struct _PDM_surf_deform_t {

  int           n_part;
  int          *n_face;
  int          *n_vtx;
  double      **coords;
  PDM_g_num_t **vtx_ln_to_gn;
  PDM_g_num_t **face_ln_to_gn;

  int         **face_vtx_idx;
  int         **face_vtx;

  double      **dcoords;
  double      **aux_geom;

} PDM_surf_deform_t;

typedef struct _PDM_cloud_deform_t {

  int           n_part;
  int          *n_points;
  double      **coords;
  PDM_g_num_t **gnum;

} PDM_cloud_deform_t;

typedef struct _PDM_mesh_deform_t {

  PDM_MPI_Comm        comm;                    /*!< MPI communicator */
  int                 is_built;                /*!< Construction flag */
  int                 n_aux_geom;
  int                 n_layer;
  int                *n_leaf_per_layer;

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
  const PDM_MPI_Comm  comm
)
{

  PDM_mesh_deform_t *def = (PDM_mesh_deform_t *) malloc(sizeof(PDM_mesh_deform_t));

  def->comm             = comm;
  def->is_built         = 0;
  def->n_aux_geom       = n_aux_geom;
  def->n_layer          = n_layer;
  def->n_leaf_per_layer = (int*) n_leaf_per_layer;

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

  if (def->is_built == 0) {
    return;
  }

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
}

static void
_PDM_mesh_deform_cloud_block_get
(
  PDM_mesh_deform_t *def
)
{
}

static void
_PDM_mesh_deform_cloud_dcoord_part_get
(
  PDM_mesh_deform_t *def
)
{
}

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
_mean_array_double_per_leaf
(
  const int      n_explicit_nodes,
  const int     *n_points,
  const int     *range,
  const int     *leaf_id,
  const int      stride,
  const double  *pts_array_double,
        int     *n_leaf,
        double **leaf_array_double
)
{

  int _n_leaf = 0;

  int *extract_leaf_node_id = malloc(n_explicit_nodes * sizeof(int));
  for(int i_node = 0; i_node < n_explicit_nodes; i_node++) {
    if(leaf_id[i_node] != -1){
      extract_leaf_node_id[_n_leaf++] = i_node;
    }
  }
  extract_leaf_node_id = realloc(extract_leaf_node_id, _n_leaf * sizeof(int));

  double *_leaf_array_double = malloc(stride * _n_leaf * sizeof(double));

  for (int i_leaf = 0; i_leaf < _n_leaf; i_leaf++) {
    _leaf_array_double[3*i_leaf    ] = 0.0;
    _leaf_array_double[3*i_leaf + 1] = 0.0;
    _leaf_array_double[3*i_leaf + 2] = 0.0;
  }

  for (int i_leaf = 0; i_leaf < _n_leaf; i_leaf++) {
    for (int i_pts = range[extract_leaf_node_id[i_leaf]]; i_pts < range[extract_leaf_node_id[i_leaf]] + n_points[extract_leaf_node_id[i_leaf]]; i_pts++) {
      for (int i = 0; i < stride; i++) {
        _leaf_array_double[stride*i_leaf+i] = _leaf_array_double[stride*i_leaf+i] + pts_array_double[stride*i_pts+i];
      }
    }
    for (int i = 0; i < stride; i++) {
      _leaf_array_double[stride*i_leaf+i] /= n_points[extract_leaf_node_id[i_leaf]];
    }
  }

  free(extract_leaf_node_id);

  *n_leaf            = _n_leaf;
  *leaf_array_double = _leaf_array_double;

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
  int         depth_max        = 50;
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
   * Move vertices
   */

  if (i_rank == 0) {
    printf("-- Move vertices\n");
    fflush(stdout);
  }

  /*
   * Get vertices new coordinates
   */

  if (i_rank == 0) {
    printf("-- Get vertices new coordinates\n");
    fflush(stdout);
  }

  for (int i_part = 0; i_part < n_part_bnd; i_part++) {

    for (int i_vtx = 0; i_vtx < n_bnd_vtx[i_part]; i_vtx++) {
      bnd_vtx[i_part][3*i_vtx   ] = bnd_vtx[i_part][3*i_vtx   ] + bnd_dvtx[i_part][3*i_vtx    ];
      bnd_vtx[i_part][3*i_vtx+ 1] = bnd_vtx[i_part][3*i_vtx+ 1] + bnd_dvtx[i_part][3*i_vtx + 1];
      bnd_vtx[i_part][3*i_vtx+ 2] = bnd_vtx[i_part][3*i_vtx+ 2] + bnd_dvtx[i_part][3*i_vtx + 2];
    }

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
