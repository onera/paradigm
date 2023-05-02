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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_extract.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_to_part.h"
#include "pdm_points_merge.h"
#include "pdm_vtk.h"
#include "pdm_gnum.h"
#include "pdm_order.h"
#include "pdm_unique.h"
#include "pdm_array.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                     argc,
           char                  **argv,
           PDM_g_num_t           *n_vtx_seg,
           PDM_Mesh_nodal_elt_t  *t_elt,
           double                *length,
           int                   *post)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *t_elt = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static
void
_generate_mesh
(
  PDM_MPI_Comm           comm,
  PDM_g_num_t            n_vtx_seg,
  double                 length,
  double                 zero_x,
  double                 zero_y,
  double                 zero_z,
  PDM_Mesh_nodal_elt_t   t_elt,
  PDM_dcube_nodal_t    **dcube_nodal_out,
  PDM_dmesh_t          **dmesh_out
)
{
  int i_rank;
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);

  int order = 1;
  if(PDM_Mesh_nodal_elmt_is_ho(t_elt) == 1) {
    order = 2;
  }

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        zero_x,
                                                        zero_y,
                                                        zero_z,
                                                        t_elt,
                                                        order,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  // PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  // // double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  // // int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_USER);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  PDM_dmesh_t* dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

  *dcube_nodal_out = dcube;
  *dmesh_out       = dmesh;

  // PDM_dmesh_free(dmesh);

  PDM_dmesh_nodal_to_dmesh_free(dmntodm);

}


static
void
_dmesh_extract_from_group_id
(
        PDM_MPI_Comm           comm,
        PDM_dmesh_t           *dmesh,
        int                    i_group,
        int                   *dn_extract_face_out,
        int                   *dn_extract_edge_out,
        int                   *dn_extract_vtx_out,
        int                  **dextract_face_edge_idx_out,
        PDM_g_num_t          **dextract_face_edge_out,
        PDM_g_num_t          **dextract_face_vtx_out,
        PDM_g_num_t          **dextract_edge_vtx_out,
        PDM_g_num_t          **dextract_face_parent_gnum_out,
        PDM_g_num_t          **dextract_edge_parent_gnum_out,
        PDM_g_num_t          **dextract_vtx_parent_gnum_out,
        double               **dextract_vtx_coord_out,
        int                    post,
  const char*                  pattern
)
{
  int i_rank;
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);

  int         *dbound_face_idx = NULL;
  PDM_g_num_t *dbound_face     = NULL;
  int n_group = PDM_dmesh_bound_get(dmesh,
                                    PDM_BOUND_TYPE_FACE,
                                    &dbound_face,
                                    &dbound_face_idx,
                                    PDM_OWNERSHIP_KEEP);
  PDM_UNUSED(n_group);

  int beg = dbound_face_idx[i_group];
  int dn_group_face = dbound_face_idx[i_group+1] - beg;

  PDM_g_num_t *selected_face_gnum = &dbound_face[beg];


  PDM_dmesh_extract_t *dme = PDM_dmesh_extract_create(2, comm);

  PDM_dmesh_extract_selected_gnum_set(dme,
                                      PDM_MESH_ENTITY_FACE,
                                      dn_group_face,
                                      selected_face_gnum);

  PDM_dmesh_extract_dmesh_set(dme, dmesh);
  PDM_dmesh_extract_compute(dme);

  /*
   * Post-traitement
   */
  PDM_dmesh_t* dmesh_extract = NULL;
  PDM_dmesh_extract_dmesh_get(dme,
                              &dmesh_extract,
                              PDM_OWNERSHIP_KEEP);

  int         *dextract_face_edge_idx = NULL;
  PDM_g_num_t *dextract_face_edge     = NULL;
  int dn_extract_face = PDM_dmesh_connectivity_get(dmesh_extract, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                   &dextract_face_edge,
                                                   &dextract_face_edge_idx,
                                                   PDM_OWNERSHIP_USER);

  int         *dextract_edge_vtx_idx = NULL;
  PDM_g_num_t *dextract_edge_vtx     = NULL;
  int dn_extract_edge = PDM_dmesh_connectivity_get(dmesh_extract, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                   &dextract_edge_vtx,
                                                   &dextract_edge_vtx_idx,
                                                   PDM_OWNERSHIP_USER);
  assert(dextract_edge_vtx_idx == NULL);

  double *dextract_vtx_coord = NULL;
  PDM_dmesh_vtx_coord_get(dmesh_extract,
                          &dextract_vtx_coord,
                          PDM_OWNERSHIP_USER);

  /* Hook all information */
  *dextract_face_edge_out     = dextract_face_edge;
  *dextract_face_edge_idx_out = dextract_face_edge_idx;

  *dextract_vtx_coord_out = dextract_vtx_coord;
  *dextract_edge_vtx_out  = dextract_edge_vtx;

  PDM_dmesh_extract_parent_gnum_get(dme,
                                    PDM_MESH_ENTITY_FACE,
                                    dn_extract_face_out,
                                    dextract_face_parent_gnum_out,
                                    PDM_OWNERSHIP_USER);
  PDM_dmesh_extract_parent_gnum_get(dme,
                                    PDM_MESH_ENTITY_EDGE,
                                    dn_extract_edge_out,
                                    dextract_edge_parent_gnum_out,
                                    PDM_OWNERSHIP_USER);
  PDM_dmesh_extract_parent_gnum_get(dme,
                                    PDM_MESH_ENTITY_VERTEX,
                                    dn_extract_vtx_out,
                                    dextract_vtx_parent_gnum_out,
                                    PDM_OWNERSHIP_USER);

  /* Compute dface_vtx */
  int         *dextract_face_vtx_idx = dextract_face_edge_idx;
  PDM_g_num_t *dextract_face_vtx     = NULL;
  PDM_compute_dface_vtx_from_edges(comm,
                                   dn_extract_face,
                                   dn_extract_edge,
                                   dextract_face_edge_idx,
                                   dextract_face_edge,
                                   dextract_edge_vtx,
                                   &dextract_face_vtx);

  *dextract_face_vtx_out = dextract_face_vtx;


  PDM_g_num_t *distrib_face_extract = PDM_compute_entity_distribution(comm, dn_extract_face);

  PDM_g_num_t* extract_face_ln_to_gn = malloc(dn_extract_face * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_extract_face; ++i) {
    extract_face_ln_to_gn[i] = distrib_face_extract[i_rank] + i + 1;
  }

  int dn_extract_vtx  = PDM_dmesh_dn_entity_get(dmesh_extract, PDM_MESH_ENTITY_VERTEX);
  PDM_g_num_t *extract_vtx_distribution = PDM_compute_entity_distribution(comm, dn_extract_vtx);

  int pn_extract_vtx = -1;
  PDM_g_num_t *pextract_vtx_ln_to_gn = NULL;
  int         *pextract_face_vtx_idx = NULL;
  int         *pextract_face_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face_extract,
                                                           dextract_face_vtx_idx,
                                                           dextract_face_vtx,
                                                           dn_extract_face,
                                                           extract_face_ln_to_gn,
                                                           &pn_extract_vtx,
                                                           &pextract_vtx_ln_to_gn,
                                                           &pextract_face_vtx_idx,
                                                           &pextract_face_vtx);

  double** tmp_pextract_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        extract_vtx_distribution,
                                        dextract_vtx_coord,
                                        &pn_extract_vtx,
                 (const PDM_g_num_t **) &pextract_vtx_ln_to_gn,
                                        &tmp_pextract_vtx_coord);
  double* pextract_vtx_coord = tmp_pextract_vtx_coord[0];
  free(tmp_pextract_vtx_coord);
  free(extract_vtx_distribution);

  if (post) {
    char filename[999];
    sprintf(filename, "%s_%i.vtk", pattern, i_rank);
    PDM_vtk_write_polydata(filename,
                           pn_extract_vtx,
                           pextract_vtx_coord,
                           pextract_vtx_ln_to_gn,
                           dn_extract_face,
                           pextract_face_vtx_idx,
                           pextract_face_vtx,
                           extract_face_ln_to_gn,
                           NULL);
  }

  // free(dextract_face_vtx);
  free(pextract_vtx_ln_to_gn);
  free(pextract_face_vtx_idx);
  free(pextract_face_vtx    );
  free(pextract_vtx_coord    );

  free(distrib_face_extract);
  free(extract_face_ln_to_gn);

  PDM_dmesh_extract_free(dme);


}

static
void
_compute_characteristic_lenght_vtx
(
        PDM_MPI_Comm           comm,
        int                    dn_face,
        int                    dn_edge,
        int                    dn_vtx,
        int                   *dface_edge_idx,
        PDM_g_num_t           *dface_edge,
        PDM_g_num_t           *dface_vtx,
        PDM_g_num_t           *dedge_vtx,
        double                *dvtx_coord,
        double               **dchar_length_out
)
{
 PDM_UNUSED(dn_face);
 PDM_UNUSED(dface_edge_idx);
 PDM_UNUSED(dface_edge);
 PDM_UNUSED(dface_vtx);

  int i_rank;
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);

  PDM_g_num_t *distrib_vtx = PDM_compute_entity_distribution (comm, dn_vtx);
  PDM_g_num_t *distrib_edge = PDM_compute_entity_distribution(comm, dn_edge);

  // Compute graph of vtx
  int         *dvtx_vtx_idx = NULL;
  PDM_g_num_t *dvtx_vtx     = NULL;
  if(dedge_vtx == NULL) {
    abort();
  } else {

    int *dedge_vtx_idx = malloc((dn_edge+1) * sizeof(int));
    for(int i = 0; i < dn_edge+1; ++i) {
      dedge_vtx_idx[i] = 2 * i;
    }
    int         *dvtx_edge_idx = NULL;
    PDM_g_num_t *dvtx_edge     = NULL;
    PDM_dconnectivity_transpose(comm,
                                distrib_edge,
                                distrib_vtx,
                                dedge_vtx_idx,
                                dedge_vtx,
                                0,
                                &dvtx_edge_idx,
                                &dvtx_edge);

    PDM_deduce_combine_connectivity_dual(comm,
                                         distrib_vtx,
                                         distrib_edge,
                                         dvtx_edge_idx,
                                         dvtx_edge,
                                         dedge_vtx_idx,
                                         dedge_vtx,
                                         0,
                                         &dvtx_vtx_idx,
                                         &dvtx_vtx);

    free(dedge_vtx_idx);
    free(dvtx_edge_idx);
    free(dvtx_edge    );
  }

  // PDM_log_trace_connectivity_long(dvtx_vtx_idx, dvtx_vtx, dn_vtx, "dvtx_vtx ::");


  /*
   * Partitionnement du pauvre
   */
  PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_vtx,
                              (const PDM_g_num_t **)  &dvtx_vtx,
                                                      &dvtx_vtx_idx[dn_vtx],
                                                      1,
                                                      comm);

  int stride_one = 1;
  double **tmp_vtx_vtx_coord = NULL;
  PDM_block_to_part_exch(btp,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                         dvtx_coord,
                         NULL,
          (void ***)    &tmp_vtx_vtx_coord);
  double *pvtx_vtx_coord = tmp_vtx_vtx_coord[0];
  free(tmp_vtx_vtx_coord);
  PDM_block_to_part_free(btp);
  free(dvtx_vtx    );

  double *char_length = malloc(dn_vtx * sizeof(double));

  double tol = 1e-6;
  const double eps_base = 1e-12;
  for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {

    char_length[i_vtx] = HUGE_VAL;

    for(int idx_vtx = dvtx_vtx_idx[i_vtx]; idx_vtx < dvtx_vtx_idx[i_vtx+1]; ++idx_vtx) {

      double length2 = 0.;
      for(int k = 0; k < 3; ++k) {
        double delta = dvtx_coord[3*i_vtx + k] - pvtx_vtx_coord[3*idx_vtx + k];
        length2 += delta * delta;
      }

      char_length[i_vtx] = PDM_MIN(char_length[i_vtx], length2);


    }

    char_length[i_vtx] = PDM_MAX(eps_base, tol*sqrt(char_length[i_vtx]));
  }


  // PDM_log_trace_array_double(char_length, dn_vtx, "char_length ::");


  *dchar_length_out = char_length;


  free(distrib_vtx );
  free(distrib_edge);
  free(pvtx_vtx_coord);

  free(dvtx_vtx_idx);


}

/**
 *
 * \brief  Main
 *
 */
int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */
  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  PDM_Mesh_nodal_elt_t t_elt   = PDM_MESH_NODAL_HEXA8;
  int post                     = 0;
  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &t_elt,
             &length,
             &post);

  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  /*
   *  Create distributed cube
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_dcube_nodal_t *dcube1 = NULL;
  PDM_dmesh_t       *dm1    = NULL;
  _generate_mesh(comm,
                 n_vtx_seg,
                 length,
                 0.,
                 0.,
                 0.,
                 t_elt,
                 &dcube1,
                 &dm1);

  PDM_dcube_nodal_t *dcube2 = NULL;
  PDM_dmesh_t       *dm2    = NULL;
  _generate_mesh(comm,
                 n_vtx_seg,
                 length/2,
                 1.,
                 0.,
                 0.,
                 t_elt,
                 &dcube2,
                 &dm2);


  int          dn_extract_face_m1 = 0;
  int          dn_extract_edge_m1 = 0;
  int          dn_extract_vtx_m1  = 0;
  int         *dextract_m1_face_edge_idx    = NULL;
  PDM_g_num_t *dextract_m1_face_edge        = NULL;
  PDM_g_num_t *dextract_m1_face_vtx         = NULL;
  PDM_g_num_t *dextract_m1_edge_vtx         = NULL;
  PDM_g_num_t *dextract_m1_face_parent_gnum = NULL;
  PDM_g_num_t *dextract_m1_edge_parent_gnum = NULL;
  PDM_g_num_t *dextract_m1_vtx_parent_gnum  = NULL;
  double      *dextract_m1_vtx_coord        = NULL;

  _dmesh_extract_from_group_id(comm,
                               dm1,
                               3,
                               &dn_extract_face_m1,
                               &dn_extract_edge_m1,
                               &dn_extract_vtx_m1,
                               &dextract_m1_face_edge_idx,
                               &dextract_m1_face_edge,
                               &dextract_m1_face_vtx,
                               &dextract_m1_edge_vtx,
                               &dextract_m1_face_parent_gnum,
                               &dextract_m1_edge_parent_gnum,
                               &dextract_m1_vtx_parent_gnum,
                               &dextract_m1_vtx_coord,
                               post,
                               "extract_mesh1_");

  int          dn_extract_face_m2 = 0;
  int          dn_extract_edge_m2 = 0;
  int          dn_extract_vtx_m2  = 0;
  int         *dextract_m2_face_edge_idx    = NULL;
  PDM_g_num_t *dextract_m2_face_edge        = NULL;
  PDM_g_num_t *dextract_m2_face_vtx         = NULL;
  PDM_g_num_t *dextract_m2_edge_vtx         = NULL;
  PDM_g_num_t *dextract_m2_face_parent_gnum = NULL;
  PDM_g_num_t *dextract_m2_edge_parent_gnum = NULL;
  PDM_g_num_t *dextract_m2_vtx_parent_gnum  = NULL;
  double      *dextract_m2_vtx_coord        = NULL;
  _dmesh_extract_from_group_id(comm,
                               dm2,
                               2,
                               &dn_extract_face_m2,
                               &dn_extract_edge_m2,
                               &dn_extract_vtx_m2,
                               &dextract_m2_face_edge_idx,
                               &dextract_m2_face_edge,
                               &dextract_m2_face_vtx,
                               &dextract_m2_edge_vtx,
                               &dextract_m2_face_parent_gnum,
                               &dextract_m2_edge_parent_gnum,
                               &dextract_m2_vtx_parent_gnum,
                               &dextract_m2_vtx_coord,
                               post,
                               "extract_mesh2_");


  // Compute carateristic lenght
  double *dchar_lenght_m1 = NULL;
  _compute_characteristic_lenght_vtx(comm,
                                     dn_extract_face_m1,
                                     dn_extract_edge_m1,
                                     dn_extract_vtx_m1,
                                     dextract_m1_face_edge_idx,
                                     dextract_m1_face_edge,
                                     dextract_m1_face_vtx,
                                     dextract_m1_edge_vtx,
                                     dextract_m1_vtx_coord,
                                     &dchar_lenght_m1);

  // Compute carateristic lenght
  double *dchar_lenght_m2 = NULL;
  _compute_characteristic_lenght_vtx(comm,
                                     dn_extract_face_m2,
                                     dn_extract_edge_m2,
                                     dn_extract_vtx_m2,
                                     dextract_m2_face_edge_idx,
                                     dextract_m2_face_edge,
                                     dextract_m2_face_vtx,
                                     dextract_m2_edge_vtx,
                                     dextract_m2_vtx_coord,
                                     &dchar_lenght_m2);

  /* Point merge */
  double tolerance = 1.e-8;
  PDM_points_merge_t* pts_merge = PDM_points_merge_create(2, tolerance, comm, PDM_OWNERSHIP_KEEP);


  PDM_points_merge_cloud_set(pts_merge, 0, dn_extract_vtx_m1, dextract_m1_vtx_coord, dchar_lenght_m1);
  PDM_points_merge_cloud_set(pts_merge, 1, dn_extract_vtx_m2, dextract_m2_vtx_coord, dchar_lenght_m2);

  PDM_points_merge_process(pts_merge);

  int *m1_candidates_idx  = NULL;
  int *m1_candidates_desc = NULL;
  PDM_points_merge_candidates_get(pts_merge,
                                  0,
                                  &m1_candidates_idx,
                                  &m1_candidates_desc); // (i_proc, i_cloud, i_point)

  log_trace("m1_candidates = \n");
  for(int i_vtx = 0; i_vtx < dn_extract_vtx_m1; ++i_vtx) {
    log_trace("i_vtx = %i link to : ", i_vtx);
    for(int j = m1_candidates_idx[i_vtx]; j < m1_candidates_idx[i_vtx+1]; ++j) {
      log_trace(" (%i / %i / %i )", m1_candidates_desc[3*j], m1_candidates_desc[3*j+1], m1_candidates_desc[3*j+2]);
    }
    log_trace("\n");
  }


  int *m2_candidates_idx  = NULL;
  int *m2_candidates_desc = NULL;
  PDM_points_merge_candidates_get(pts_merge,
                                  1,
                                  &m2_candidates_idx,
                                  &m2_candidates_desc); // (i_proc, i_cloud, i_point)

  log_trace("m2_candidates = \n");
  for(int i_vtx = 0; i_vtx < dn_extract_vtx_m2; ++i_vtx) {
    log_trace("i_vtx = %i link to : ", i_vtx);
    for(int j = m2_candidates_idx[i_vtx]; j < m2_candidates_idx[i_vtx+1]; ++j) {
      log_trace(" (%i / %i / %i )", m2_candidates_desc[3*j], m2_candidates_desc[3*j+1], m2_candidates_desc[3*j+2]);
    }
    log_trace("\n");
  }

  // Create a new global numbering to prepare ptp
  int n_point_cloud = 2;
  int           n_elt1            [2] = {dn_extract_vtx_m1, dn_extract_vtx_m2};
  PDM_g_num_t  *part1_gnum        [2] = {dextract_m1_vtx_parent_gnum, dextract_m2_vtx_parent_gnum};

  PDM_gen_gnum_t *gnum = PDM_gnum_create(3, n_point_cloud, PDM_TRUE, 1.e-3, comm, PDM_OWNERSHIP_USER);

  PDM_gnum_set_parents_nuplet(gnum, 2); // (i_cloud / Gnum )

  PDM_g_num_t **part1_nuplet = malloc(n_point_cloud * sizeof(PDM_g_num_t *));
  for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {
    part1_nuplet[i_cloud] = malloc( 2 * n_elt1[i_cloud] * sizeof(PDM_g_num_t));

    for(int i = 0; i < n_elt1[i_cloud]; ++i) {
      part1_nuplet[i_cloud][2*i  ] = i_cloud;
      part1_nuplet[i_cloud][2*i+1] = part1_gnum[i_cloud][i];
    }

    PDM_gnum_set_from_parents(gnum, i_cloud, n_elt1[i_cloud], part1_nuplet[i_cloud]);

  }
  PDM_gnum_compute(gnum);

  PDM_g_num_t **part1_concat_gnum =  malloc(n_point_cloud * sizeof(PDM_g_num_t *));
  int         **part1_cloud       =  malloc(n_point_cloud * sizeof(int         *));
  for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {
    part1_concat_gnum[i_cloud] = PDM_gnum_get(gnum, i_cloud);
    free(part1_nuplet[i_cloud]);
    PDM_log_trace_array_long(part1_concat_gnum[i_cloud], n_elt1[i_cloud], "part1_concat_gnum ::");

    part1_cloud[i_cloud] = malloc(n_elt1[i_cloud] * sizeof(int));
    for(int i_vtx = 0; i_vtx < n_elt1[i_cloud]; ++i_vtx) {
      part1_cloud[i_cloud][i_vtx] = i_cloud;
    }

  }
  free(part1_nuplet);

  PDM_gnum_free(gnum);


  // Create ptp to exchange links
  int          n_elt2            [2] = {dn_extract_vtx_m1, dn_extract_vtx_m2};
  int         *part1_to_part2_idx[2] = {m1_candidates_idx, m2_candidates_idx};
  int         *part1_to_part2    [2] = {m1_candidates_desc, m2_candidates_desc};

  for(int i_vtx = 0; i_vtx < dn_extract_vtx_m1+1; ++i_vtx) {
    m1_candidates_idx[i_vtx] *= 3;
  }
  for(int i_vtx = 0; i_vtx < dn_extract_vtx_m2+1; ++i_vtx) {
    m2_candidates_idx[i_vtx] *= 3;
  }

  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) part1_concat_gnum,
                                                                      (const int          *) n_elt1,
                                                                                             n_point_cloud,
                                                                      (const int          *) n_elt2,
                                                                                             n_point_cloud,
                                                                      (const int         **) part1_to_part2_idx,
                                                                                             NULL,
                                                                      (const int         **) part1_to_part2,
                                                                                             comm);

  int request = -1;
  PDM_g_num_t **vtx_opp_gnum = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(PDM_g_num_t),
                         NULL,
        (const void  **) part1_gnum,
                         NULL,
        (void       ***) &vtx_opp_gnum,
                         &request);

  PDM_part_to_part_iexch_wait(ptp, request);

  int **vtx_opp_cloud = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(int),
                         NULL,
        (const void  **) part1_cloud,
                         NULL,
        (void       ***) &vtx_opp_cloud,
                         &request);

  PDM_part_to_part_iexch_wait(ptp, request);


  int  *n_ref = NULL;
  int **ref   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp,
                                 &n_ref,
                                 &ref);

  int  *n_unref = NULL;
  int **unref   = NULL;
  PDM_part_to_part_unref_lnum2_get(ptp,
                                   &n_unref,
                                   &unref);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp,
                                       &gnum1_come_from_idx,
                                       &gnum1_come_from);

  if(1 == 1) {
    for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {
      printf("n_unref[%i] = %i \n", i_cloud, n_unref[i_cloud]);

      // for (int i = 0; i < n_ref[i_cloud]; i++) {
      //   int elt_a_id = ref[i_cloud][i] - 1;
      // }
      // vtx_opp_gnum == gnum1_come_from
      PDM_log_trace_array_long(vtx_opp_gnum   [i_cloud], n_ref [i_cloud], "vtx_opp_gnum    :");
      PDM_log_trace_array_long(part1_gnum     [i_cloud], n_elt1[i_cloud], "part1_gnum      :");
      PDM_log_trace_array_long(gnum1_come_from[i_cloud], gnum1_come_from_idx[i_cloud][n_ref[i_cloud]], "gnum1_come_from :");
    }
  }


  for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {
    free(part1_concat_gnum[i_cloud]);
    free(part1_cloud[i_cloud]);
  }
  free(part1_concat_gnum);
  free(part1_cloud);


  // Generate interface_gnum
  PDM_gen_gnum_t *gnum_itrf = PDM_gnum_create(3, n_point_cloud, PDM_TRUE, 1.e-3, comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_parents_nuplet(gnum_itrf, 2); // (max(i_cloud, i_cloud_opp) / min(i_cloud, i_cloud_opp) )

  int          *n_itrf    = malloc(n_point_cloud * sizeof(int          ));
  PDM_g_num_t **itrf_pair = malloc(n_point_cloud * sizeof(PDM_g_num_t *));
  PDM_g_num_t **itrf_gnum = malloc(n_point_cloud * sizeof(PDM_g_num_t *));
  for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {

    int *_gnum1_come_from_idx = gnum1_come_from_idx[i_cloud];
    int n_come_from = _gnum1_come_from_idx[n_ref[i_cloud]];

    n_itrf[i_cloud] = n_come_from;

    itrf_pair[i_cloud] = malloc(2 * n_come_from * sizeof(PDM_g_num_t));
    for(int idx_vtx = 0; idx_vtx < n_ref[i_cloud]; ++idx_vtx) {
      // int i_vtx = ref[i_cloud][idx_vtx] - 1;
      for(int k = _gnum1_come_from_idx[idx_vtx]; k < _gnum1_come_from_idx[idx_vtx+1]; ++k) {
        int i_cloud_opp = vtx_opp_cloud[i_cloud][k];
        itrf_pair[i_cloud][2*k  ] = PDM_MIN(i_cloud, i_cloud_opp);
        itrf_pair[i_cloud][2*k+1] = PDM_MAX(i_cloud, i_cloud_opp);
      }
    }
    PDM_gnum_set_from_parents(gnum_itrf, i_cloud, n_come_from, itrf_pair[i_cloud]);
    // PDM_log_trace_array_long(itrf_pair[i_cloud], 2 * n_elt1[i_cloud], "itrf_pair (Before):");
  }

  PDM_gnum_compute(gnum_itrf);

  for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {
    itrf_gnum[i_cloud] = PDM_gnum_get(gnum_itrf, i_cloud);
    PDM_log_trace_array_long(itrf_gnum[i_cloud], n_itrf[i_cloud], "itrf_gnum :");
  }
  PDM_gnum_free(gnum_itrf);


  /*
   * Create unique id
   */
  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                                      itrf_gnum,
                                                      NULL,
                                                      n_itrf,
                                                      n_point_cloud,
                                                      comm);

  int *ditrf_pair = NULL;
  PDM_part_to_block_exch(ptb,
                         2 * sizeof(int),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
           (void **)     itrf_pair,
                         NULL,
           (void **)     &ditrf_pair);


  free(n_itrf);
  for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {
    free(itrf_pair[i_cloud]);
  }
  free(itrf_pair);

  int dn_interface = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *distrib_itrf_gnum   = PDM_compute_entity_distribution(comm, dn_interface);
  int         *distrib_itrf   = malloc((n_rank+1) * sizeof(int));
  int         *distrib_itrf_n = malloc( n_rank    * sizeof(int));

  for(int i = 0; i < n_rank+1; ++i) {
    distrib_itrf[i] = distrib_itrf_gnum[i];
  }
  free(distrib_itrf_gnum);

  for(int i = 0; i < n_rank; ++i) {
    distrib_itrf_n[i] = distrib_itrf[i+1] - distrib_itrf[i];
  }

  PDM_g_num_t *ditrf_gnum = PDM_part_to_block_block_gnum_get(ptb);

  PDM_log_trace_array_int(distrib_itrf_n, n_rank  , "distrib_itrf_n ::");
  PDM_log_trace_array_int(distrib_itrf  , n_rank+1, "distrib_itrf   ::");


  PDM_g_num_t *all_itrf_gnum = malloc(    distrib_itrf[n_rank] * sizeof(PDM_g_num_t));
  int         *all_itrf_pair = malloc(2 * distrib_itrf[n_rank] * sizeof(int        ));
  PDM_MPI_Allgatherv(ditrf_gnum,
                     dn_interface,
                     PDM__PDM_MPI_G_NUM,
                     all_itrf_gnum,
                     distrib_itrf_n,
                     distrib_itrf,
                     PDM__PDM_MPI_G_NUM,
                     comm);

  for(int i = 0; i < n_rank; ++i) {
    distrib_itrf  [i] = distrib_itrf  [i] * 2;
    distrib_itrf_n[i] = distrib_itrf_n[i] * 2;
  }
  PDM_MPI_Allgatherv(ditrf_pair,
                     2 * dn_interface,
                     PDM_MPI_INT,
                     all_itrf_pair,
                     distrib_itrf_n,
                     distrib_itrf,
                     PDM_MPI_INT,
                     comm);
  for(int i = 0; i < n_rank; ++i) {
    distrib_itrf  [i] = distrib_itrf  [i] / 2;
    distrib_itrf_n[i] = distrib_itrf_n[i] / 2;
  }


  PDM_g_num_t n_g_interface = distrib_itrf[n_rank];
  free(distrib_itrf);
  free(distrib_itrf_n);
  free(ditrf_pair);
  PDM_part_to_block_free(ptb);

  printf("n_g_interface = %i | dn_interface = %i \n", n_g_interface, dn_interface);
  PDM_MPI_Barrier(comm);

  // Bucket sort by interface
  int *entity_itrf_n = PDM_array_zeros_int(n_g_interface);

  for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {

    int *_gnum1_come_from_idx = gnum1_come_from_idx[i_cloud];
    for(int idx_vtx = 0; idx_vtx < n_ref[i_cloud]; ++idx_vtx) {
      // int i_vtx = ref[i_cloud][idx_vtx] - 1;
      for(int k = _gnum1_come_from_idx[idx_vtx]; k < _gnum1_come_from_idx[idx_vtx+1]; ++k) {
        int i_cloud_opp = vtx_opp_cloud[i_cloud][k];
        if(i_cloud < i_cloud_opp) {
          entity_itrf_n[itrf_gnum[i_cloud][k]-1]++;
        }
      }
    }
  }

  int *entity_itrf_idx = malloc((n_g_interface+1) * sizeof(int));
  entity_itrf_idx[0] = 0;
  for(int i = 0; i < n_g_interface; ++i) {
    entity_itrf_idx[i+1] = entity_itrf_idx[i] + entity_itrf_n[i];
    entity_itrf_n  [i  ] = 0;
  }

  PDM_g_num_t *concat_vtx_cur = malloc(entity_itrf_idx[n_g_interface] * sizeof(PDM_g_num_t));
  PDM_g_num_t *concat_vtx_opp = malloc(entity_itrf_idx[n_g_interface] * sizeof(PDM_g_num_t));

  for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {
    int *_gnum1_come_from_idx = gnum1_come_from_idx[i_cloud];
    for(int idx_vtx = 0; idx_vtx < n_ref[i_cloud]; ++idx_vtx) {
      int i_vtx = ref[i_cloud][idx_vtx] - 1;
      for(int k = _gnum1_come_from_idx[idx_vtx]; k < _gnum1_come_from_idx[idx_vtx+1]; ++k) {
        int i_cloud_opp = vtx_opp_cloud[i_cloud][k];
        PDM_g_num_t g_itrf = itrf_gnum[i_cloud][k]-1;
        if(i_cloud < i_cloud_opp) {
          int idx_write = entity_itrf_idx[g_itrf] + entity_itrf_n[g_itrf]++;
          concat_vtx_cur[idx_write] = vtx_opp_gnum[i_cloud][k];
          concat_vtx_opp[idx_write] = part1_gnum  [i_cloud][i_vtx];
        }
      }
    }
  }


  PDM_part_to_part_free(ptp);

  /* Let's go */
  int          *dn_vtx_itrf   = malloc(n_g_interface * sizeof(int          ));
  PDM_g_num_t **itrf_gnum_cur = malloc(n_g_interface * sizeof(PDM_g_num_t *));
  PDM_g_num_t **itrf_gnum_opp = malloc(n_g_interface * sizeof(PDM_g_num_t *));
  for(int i_itrf = 0; i_itrf < n_g_interface; ++i_itrf) {

    int beg = entity_itrf_idx[i_itrf];
    int pn_vtx = entity_itrf_idx[i_itrf+1] - beg;

    PDM_g_num_t *lconcat_vtx_cur = &concat_vtx_cur[beg];
    PDM_g_num_t *lconcat_vtx_opp = &concat_vtx_opp[beg];

    PDM_part_to_block_t* ptb_itrf = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                             PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                             1.,
                                                             &lconcat_vtx_cur,
                                                             NULL,
                                                             &pn_vtx,
                                                             1,
                                                             comm);

    dn_vtx_itrf[i_itrf] = PDM_part_to_block_n_elt_block_get(ptb_itrf);

    PDM_part_to_block_exch(ptb_itrf,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
             (void **)     &lconcat_vtx_opp,
                           NULL,
             (void **)     &itrf_gnum_opp[i_itrf]);

    PDM_part_to_block_exch(ptb_itrf,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
             (void **)     &lconcat_vtx_cur,
                           NULL,
             (void **)     &itrf_gnum_cur[i_itrf]);


    if(1 == 1) {
      PDM_log_trace_array_long(itrf_gnum_cur[i_itrf], dn_vtx_itrf[i_itrf], "itrf_gnum_cur :");
      PDM_log_trace_array_long(itrf_gnum_opp[i_itrf], dn_vtx_itrf[i_itrf], "itrf_gnum_opp :");
    }


    PDM_part_to_block_free(ptb_itrf);
  }

  free(concat_vtx_cur);
  free(concat_vtx_opp);

  for(int i_itrf = 0; i_itrf < n_g_interface; ++i_itrf) {
    free(itrf_gnum_cur[i_itrf]);
    free(itrf_gnum_opp[i_itrf]);
  }
  free(itrf_gnum_cur);
  free(itrf_gnum_opp);
  free(dn_vtx_itrf);


  PDM_log_trace_array_long(all_itrf_gnum,     n_g_interface, "all_itrf_gnum :");
  PDM_log_trace_array_int (all_itrf_pair, 2 * n_g_interface, "all_itrf_pair :");

  free(all_itrf_gnum);
  free(all_itrf_pair);
  for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {
    free(itrf_gnum[i_cloud]);
  }
  free(itrf_gnum);

  free(entity_itrf_idx);
  free(entity_itrf_n  );

  for(int i_cloud = 0; i_cloud < n_point_cloud; ++i_cloud) {
    free(vtx_opp_gnum [i_cloud]);
    free(vtx_opp_cloud[i_cloud]);
  }
  free(vtx_opp_gnum );
  free(vtx_opp_cloud);

  // Deduce face

  // Update in gnum frame

  PDM_points_merge_free(pts_merge);

  free(dextract_m1_face_edge_idx   );
  free(dextract_m1_face_edge       );
  free(dextract_m1_face_vtx        );
  free(dextract_m1_edge_vtx        );
  free(dextract_m1_face_parent_gnum);
  free(dextract_m1_edge_parent_gnum);
  free(dextract_m1_vtx_parent_gnum );
  free(dextract_m1_vtx_coord       );

  free(dextract_m2_face_edge_idx   );
  free(dextract_m2_face_edge       );
  free(dextract_m2_face_vtx        );
  free(dextract_m2_edge_vtx        );
  free(dextract_m2_face_parent_gnum);
  free(dextract_m2_edge_parent_gnum);
  free(dextract_m2_vtx_parent_gnum );
  free(dextract_m2_vtx_coord       );


  free(dchar_lenght_m1);
  free(dchar_lenght_m2);


  PDM_dmesh_free(dm1);
  PDM_dcube_nodal_gen_free(dcube1);
  PDM_dmesh_free(dm2);
  PDM_dcube_nodal_gen_free(dcube2);

  PDM_MPI_Finalize();

  return 0;
}
