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
#include "pdm_vtk.h"
#include "pdm_array.h"
#include "pdm_logging.h"

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
        *t_elt = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
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
  PDM_g_num_t          n_vtx_seg = 10;
  double               length    = 1.;
  PDM_Mesh_nodal_elt_t t_elt     = PDM_MESH_NODAL_HEXA8;
  int                  post      = 0;
  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &t_elt,
             &length,
             &post);

  int order = 1;
  if(PDM_Mesh_nodal_elmt_is_ho(t_elt) == 1) {
    order = 2;
  }

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

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        t_elt,
                                                        order,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);

  if(post) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  PDM_dmesh_t* dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

  int         *dface_vtx_idx;
  PDM_g_num_t *dface_vtx;
  PDM_dmesh_connectivity_get(dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                             &dface_vtx,
                             &dface_vtx_idx,
                             PDM_OWNERSHIP_KEEP);

  /* Extract */
  PDM_g_num_t *dgroup_face     = NULL;
  int         *dgroup_face_idx = NULL;
  int n_face_group = PDM_dmesh_bound_get(dmesh,
                                         PDM_BOUND_TYPE_FACE,
                                         &dgroup_face,
                                         &dgroup_face_idx,
                                         PDM_OWNERSHIP_BAD_VALUE);

  PDM_g_num_t *distrib_face = NULL;
  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE, &distrib_face);

  // if (1 == 1){
  //   int _n_face = distrib_face[i_rank+1]-distrib_face[i_rank];
  //   log_trace("n_face       = %d\n", _n_face);
  //   log_trace("n_face_group = %d\n", n_face_group);
  //   PDM_log_trace_array_long(distrib_face   ,  n_rank+1                     , "distrib_face    ::");
  //   PDM_log_trace_array_int (dface_vtx_idx  , _n_face+1                     , "dface_vtx_idx   ::");
  //   PDM_log_trace_array_long(dface_vtx      ,  dface_vtx_idx[_n_face]       , "dface_vtx       ::");
  //   PDM_log_trace_array_int (dgroup_face_idx,  n_face_group+1               , "dgroup_face_idx ::");
  //   PDM_log_trace_array_long(dgroup_face    ,  dgroup_face_idx[n_face_group], "dgroup_face     ::");
  // }

  int           n_group_ridge         = 0;
  PDM_g_num_t  *distrib_ridge         = NULL;
  PDM_g_num_t  *dridge_vtx            = NULL;
  int          *dgroup_edge_idx       = NULL;
  PDM_g_num_t  *dgroup_edge           = NULL;
  int          *dridge_face_group_idx = NULL;
  int          *dridge_face_group     = NULL;
  PDM_dmesh_find_topological_ridges(comm,
                                    distrib_face,
                                    dface_vtx_idx,
                                    dface_vtx,
                                    n_face_group,
                                    dgroup_face_idx,
                                    dgroup_face,
                                    &distrib_ridge,
                                    &dridge_vtx,
                                    &n_group_ridge,
                                    &dgroup_edge_idx,
                                    &dgroup_edge,
                                    &dridge_face_group_idx,
                                    &dridge_face_group);


  if (post) {

    int dn_extract_edge = distrib_ridge[i_rank+1] - distrib_ridge[i_rank];

    int *dridge_vtx_idx;
    PDM_malloc(dridge_vtx_idx,(dn_extract_edge+1) ,int);
    PDM_g_num_t *ridge_ln_to_gn;
    PDM_malloc(ridge_ln_to_gn,dn_extract_edge ,PDM_g_num_t);

    dridge_vtx_idx[0] = 0;
    for(int i = 0; i < dn_extract_edge; ++i) {
      dridge_vtx_idx[i+1] = dridge_vtx_idx[i] + 2;
      ridge_ln_to_gn[i] = distrib_ridge[i_rank] + i + 1;
    }

    int          pn_extract_vtx = -1;
    PDM_g_num_t *pextract_vtx_ln_to_gn = NULL;
    int         *pextract_edge_vtx_idx = NULL;
    int         *pextract_edge_vtx     = NULL;
    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                             distrib_ridge,
                                                             dridge_vtx_idx,
                                                             dridge_vtx,
                                                             dn_extract_edge,
                                                             ridge_ln_to_gn,
                                                             &pn_extract_vtx,
                                                             &pextract_vtx_ln_to_gn,
                                                             &pextract_edge_vtx_idx,
                                                             &pextract_edge_vtx);

    double** tmp_pextract_vtx_coord = NULL;
    PDM_part_dcoordinates_to_pcoordinates(comm,
                                          1,
                                          vtx_distrib,
                                          dvtx_coord,
                                          &pn_extract_vtx,
                   (const PDM_g_num_t **) &pextract_vtx_ln_to_gn,
                                          &tmp_pextract_vtx_coord);
    double* pextract_vtx_coord = tmp_pextract_vtx_coord[0];
   PDM_free(tmp_pextract_vtx_coord);

    /* Recuperation des couleurs de ridges */
    int         **tmp_pgroup_ridge_idx = NULL;
    int         **tmp_pgroup_ridge     = NULL;
    PDM_g_num_t **tmp_pgroup_ln_to_gn  = NULL;
    PDM_part_distgroup_to_partgroup(comm,
                                    distrib_ridge,
                                    n_group_ridge,
                                    dgroup_edge_idx,
                                    dgroup_edge,
                                    1,
                                    &dn_extract_edge,
            (const PDM_g_num_t **)  &ridge_ln_to_gn,
                                    &tmp_pgroup_ridge_idx,
                                    &tmp_pgroup_ridge,
                                    &tmp_pgroup_ln_to_gn);
    int         *pgroup_ridge_idx = tmp_pgroup_ridge_idx[0];
    int         *pgroup_ridge     = tmp_pgroup_ridge    [0];
    PDM_g_num_t *pgroup_ln_to_gn  = tmp_pgroup_ln_to_gn [0];

   PDM_free(tmp_pgroup_ridge_idx);
   PDM_free(tmp_pgroup_ridge    );
   PDM_free(tmp_pgroup_ln_to_gn );

    int *pedge_id = PDM_array_const_int(dn_extract_edge, -1);

    for(int i_group = 0; i_group < n_group_ridge; ++i_group) {
      for(int i = pgroup_ridge_idx[i_group]; i < pgroup_ridge_idx[i_group+1]; ++i) {
        pedge_id[pgroup_ridge[i]-1] = i_group;
      }
    }

    const char* field_name[] = {"group", 0 };
    const int*  field     [] = {pedge_id, 0 };
    char filename[999];
    sprintf(filename, "topological_ridges_%i.vtk", i_rank);
    PDM_vtk_write_std_elements(filename,
                               pn_extract_vtx,
                               pextract_vtx_coord,
                               pextract_vtx_ln_to_gn,
                               PDM_MESH_NODAL_BAR2,
                               dn_extract_edge,
                               pextract_edge_vtx,
                               ridge_ln_to_gn,
                               1,
                               field_name,
                               field);

   PDM_free(dridge_vtx_idx);
   PDM_free(ridge_ln_to_gn);
   PDM_free(pextract_vtx_coord);
   PDM_free(pextract_vtx_ln_to_gn);
   PDM_free(pextract_edge_vtx_idx);
   PDM_free(pextract_edge_vtx);

   PDM_free(pgroup_ridge_idx);
   PDM_free(pgroup_ridge    );
   PDM_free(pgroup_ln_to_gn );
   PDM_free(pedge_id        );

  }


 PDM_free(distrib_ridge  );
 PDM_free(dridge_vtx     );
 PDM_free(dgroup_edge_idx);
 PDM_free(dgroup_edge    );
 PDM_free(dridge_face_group_idx);
 PDM_free(dridge_face_group    );


  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_dcube_nodal_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
