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
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_para_graph_dual.h"
#include "pdm_part.h"
#include "pdm_multipart.h"
#include "pdm_dconnectivity_transform.h"

#include "pdm_writer.h"
#include "pdm_geom_elem.h"
#include "pdm_part_to_block.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_overlay.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi_node_first_rank.h"

#include "pdm_triangulate.h"

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
_usage
(
 int exit_code
 )
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -nA       <value>  Cube A - Number of vertices on the side.\n\n"
     "  -lA       <value>  Cube A - length.\n\n"
     "  -xminA    <value>  Cube A - Origin X.\n\n"
     "  -yminA    <value>  Cube A - Origin Y.\n\n"
     "  -yminA    <value>  Cube A - Origin Y.\n\n"
     "  -n_partA  <level>  Cube A - Number of partitions par process.\n\n"
     "  -nB       <value>  Cube A - Number of vertices on the side (default : nA+4).\n\n"
     "  -lB       <value>  Cube A - length.\n\n"
     "  -xminB    <value>  Cube A - Origin X.\n\n"
     "  -yminB    <value>  Cube A - Origin Y.\n\n"
     "  -yminB    <value>  Cube A - Origin Y.\n\n"
     "  -n_partB  <level>  Cube A - Number of partitions par process.\n\n"
     "  -post              Ensight output.\n\n"
     "  -no_random         No random to define points coordinates.\n\n"
     "  -random_time_init  Intialize random with the current time.\n\n"
     "  -parmetis          Call ParMETIS.\n\n"
     "  -pt-scocth         Call PT-Scotch.\n\n"
     "  -n_proc_data       Number of processses where there are data.\n\n"
     "  -h                 This message.\n\n");

  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   nVtxSeg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *nVtxSegA,
 double        *lengthA,
 double        *xminA,
 double        *yminA,
 int           *n_partA,
 PDM_g_num_t   *nVtxSegB,
 double        *lengthB,
 double        *xminB,
 double        *yminB,
 int           *n_partB,
 double        *depth,
 int           *post,
 int           *method,
 int           *haveRandom,
 int           *randomTimeInit,
 int           *randomMeshAInit,
 int           *randomMeshBInit,
 int           *nProcData
 )
{
  int i = 1;

  /* Parse and check command line */

  int passB = 0;

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp (argv[i], "-nA") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _nVtxSeg = atol (argv[i]);
        *nVtxSegA = (PDM_g_num_t) _nVtxSeg;
        if (!passB) {
          *nVtxSegB = *nVtxSegA +4;
        }
      }
    }
    else if (strcmp (argv[i], "-lA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *lengthA = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-xminA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *xminA = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-yminA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *yminA = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-n_partA") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *n_partA = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-randomMeshAInit") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *randomMeshAInit = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-nB") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _nVtxSeg = atol (argv[i]);
        *nVtxSegB = (PDM_g_num_t) _nVtxSeg;
        passB = 1;
      }
    }
    else if (strcmp (argv[i], "-lB") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *lengthB = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-xminB") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *xminB = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-yminB") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *yminB = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-n_partB") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *n_partB = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-depth") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else
        *depth = atof (argv[i]);
    }
    else if (strcmp (argv[i], "-randomMeshBInit") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *randomMeshBInit = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-no_random") == 0) {
      *haveRandom = 0;
    }
    else if (strcmp (argv[i], "-random_time_init") == 0) {
      *randomTimeInit = 1;
    }
    else if (strcmp (argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp (argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp (argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else if (strcmp (argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *nProcData = atoi (argv[i]);
      }
    }
    else
      _usage (EXIT_FAILURE);
    i++;
  }
}



static void _add_depth (const double  x_min,
                        const double  x_max,
                        const double  scale,
                        const int     n_pts,
                        double       *coord)
{
  PDM_UNUSED (x_min);
  PDM_UNUSED (x_max);

  for (int i = 0; i < n_pts; i++) {
    double x = coord[3*i];
    double y = coord[3*i+1];
    //coord[3*i+2] = scale * (x*x + y*y);
    coord[3*i+2] = 0.5 * scale * (cos(6*x + .2) + sin(5*y + .1));
  }
}




static void
_compute_face_vtx
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
 )
{
  *face_vtx = malloc (sizeof(int) * face_edge_idx[n_face]);

  int max_n_vtx = 0;
  for (int i = 0; i < n_face; i++) {
    int n_vtx = face_edge_idx[i+1] - face_edge_idx[i];
    max_n_vtx = PDM_MAX (max_n_vtx, n_vtx);
  }

  int *edge_used = malloc (sizeof(int) * max_n_vtx);

  for (int i = 0; i < n_face; i++) {
    int *fe = face_edge + face_edge_idx[i];
    int *fv = *face_vtx + face_edge_idx[i];

    int n_vtx = face_edge_idx[i+1] - face_edge_idx[i];
    for (int j = 0; j < n_vtx; j++) {
      edge_used[j] = 0;
    }

    int v0, v1;
    int iedge = fe[0];
    if (iedge < 0) {
      iedge = -iedge - 1;
      v0 = edge_vtx[2*iedge+1];
      v1 = edge_vtx[2*iedge];
    } else {
      iedge = iedge - 1;
      v0 = edge_vtx[2*iedge];
      v1 = edge_vtx[2*iedge+1];
    }
    fv[0] = v0;
    edge_used[0] = 1;

    for (int j = 1; j < n_vtx; j++) {
      fv[j] = v1;
      if (j == n_vtx-1) break;

      for (int k = 0; k < n_vtx; k++) {
        if (edge_used[k]) continue;

        iedge = fe[k];
        if (iedge < 0) {
          iedge = -iedge - 1;
          if (edge_vtx[2*iedge+1] == v1) {
            v1 = edge_vtx[2*iedge];
            edge_used[k] = 1;
            break;
          }
        } else {
          iedge = iedge - 1;
          if (edge_vtx[2*iedge] == v1) {
            v1 = edge_vtx[2*iedge+1];
            edge_used[k] = 1;
            break;
          }
        }
      }
    }
  }
}


static void
_split_multipart
(
 const PDM_MPI_Comm        comm,
 const PDM_split_dual_t    method,
 const int                 n_part,
 const int                 dn_face,
 const int                 dn_edge,
 const int                 dn_vtx,
 const int                 n_edge_group,
 PDM_g_num_t              *dedge_face,
 PDM_g_num_t              *dedge_vtx,
 double                   *dvtx_coord,
 int                      *dedge_group_idx,
 PDM_g_num_t              *dedge_group,
 int                     **n_face,
 int                    ***face_vtx_idx,
 int                    ***face_vtx,
 PDM_g_num_t            ***face_g_num,
 int                     **n_vtx,
 double                 ***vtx_coord,
 PDM_g_num_t            ***vtx_g_num
 )
{
  /* Initialize multipart */
  int n_zone = 1;
  int *n_part_zones  = (int *) malloc(n_zone * sizeof(int));
  for (int i_zone = 0; i_zone < n_zone; i_zone++){
    n_part_zones[i_zone] = n_part;
  }
  int mpart_id = PDM_multipart_create (n_zone,
                                       n_part_zones,
                                       PDM_FALSE,
                                       method,
                                       PDM_PART_SIZE_HOMOGENEOUS,
                                       NULL,
                                       comm,
                                       PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options (mpart_id,
                                        -1,
                                        "PDM_PART_RENUM_CELL_NONE",
                                        NULL,
                                        "PDM_PART_RENUM_FACE_NONE");




  int n_bound = n_edge_group;
  int n_join = 0;

  int *djoins_ids     = (int *) malloc(n_join * sizeof(int));
  int *dedge_bnd_idx  = (int *) malloc((n_bound + 1) * sizeof(int));
  int *dedge_join_idx = (int *) malloc((n_join  + 1) * sizeof(int));
  dedge_bnd_idx[0] = 0;
  dedge_join_idx[0] = 0;

  // First pass to count and allocate
  int i_bnd = 1;
  for (int igroup = 0; igroup < n_edge_group; igroup++) {
    int group_size = dedge_group_idx[igroup+1] - dedge_group_idx[igroup];
    dedge_bnd_idx[i_bnd++] = group_size;
  }
  for (int i = 0; i < n_bound; i++) {
    dedge_bnd_idx[i+1] = dedge_bnd_idx[i+1] + dedge_bnd_idx[i];
  }

  // Second pass to copy
  PDM_g_num_t *dedge_bnd  = (PDM_g_num_t *) malloc(dedge_bnd_idx[n_bound] * sizeof(PDM_g_num_t));
  PDM_g_num_t *dedge_join = (PDM_g_num_t *) malloc(dedge_join_idx[n_join] * sizeof(PDM_g_num_t));

  i_bnd = 0;
  for (int igroup = 0; igroup < n_edge_group; igroup++) {
    for (int i = dedge_group_idx[igroup]; i < dedge_group_idx[igroup+1]; i++) {
      dedge_bnd[i_bnd++] = dedge_group[i];
    }
  }

  PDM_dmesh_t *dmesh = PDM_dmesh_create (PDM_OWNERSHIP_KEEP,
                                         dn_face,
                                         dn_edge,
                                         -1,
                                         dn_vtx,
                                         n_bound,
                                         n_join,
                                         comm);

  int *dedge_vtx_idx = (int *) malloc((dn_edge + 1) * sizeof(int));
  dedge_vtx_idx[0] = 0;
  for (int i = 0; i < dn_edge; i++) {
    dedge_vtx_idx[i+1] = dedge_vtx_idx[i] + 2;
  }

  PDM_dmesh_set (dmesh,
                 dvtx_coord,
                 dedge_vtx_idx,
                 dedge_vtx,
                 dedge_face,
                 dedge_bnd_idx,
                 dedge_bnd,
                 djoins_ids,
                 dedge_join_idx,
                 dedge_join);

  PDM_multipart_register_block (mpart_id, 0, dmesh);

  /* Connection between zones */
  int n_total_joins = 0;
  int *join_to_opposite = (int *) malloc(n_total_joins*sizeof(int));
  PDM_multipart_register_joins (mpart_id, n_total_joins, join_to_opposite);

  /* Run */
  PDM_multipart_run_ppart (mpart_id);

  /* Get parts */
  *n_face = (int *) malloc (sizeof(int) * n_part);
  *n_vtx  = (int *) malloc (sizeof(int) * n_part);

  *face_vtx_idx = (int **)         malloc (sizeof(int *)         * n_part);
  *face_vtx     = (int **)         malloc (sizeof(int *)         * n_part);
  *face_g_num   = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part);
  *vtx_coord    = (double **)      malloc (sizeof(double *)      * n_part);
  *vtx_g_num    = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_proc, tn_part;
    int _n_face, _n_edge, _n_vtx, n_bounds, n_joins, n_part_joins;
    int sface_edge, sedge_vtx, sedge_bound, sedge_join;
    int  n_section;
    int* n_elt;

    PDM_multipart_part_dim_get(mpart_id, 0, i_part, &n_section, &n_elt,
                               &_n_face, &_n_edge, &n_part_joins, &_n_vtx, &n_proc, &tn_part,
                               &sface_edge, &sedge_vtx, &sedge_bound, &n_bounds, &sedge_join, &n_joins);

    double       *_vtx;
    int          *_face_edge_idx, *_face_edge, *_edge_face, *_edge_vtx_idx, *_edge_vtx;
    int          *edge_bound_idx, *edge_bound, *edge_join_idx, *edge_join;
    int          *edge_part_bound_proc_idx, *edge_part_bound_part_idx, *edge_part_bound;
    PDM_g_num_t  *_face_ln_to_gn, *_edge_ln_to_gn, *_vtx_ln_to_gn, *edge_bound_ln_to_gn, *edge_join_ln_to_gn;
    int          *face_tag, *edge_tag, *vtx_tag;
    int         **elt_vtx_idx;
    int         **elt_vtx;
    PDM_g_num_t **elt_section_ln_to_gn;

    PDM_multipart_part_val_get(mpart_id, 0, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                               &face_tag, &_face_edge_idx, &_face_edge, &_face_ln_to_gn,
                               &edge_tag, &_edge_face, &_edge_vtx_idx, &_edge_vtx, &_edge_ln_to_gn,
                               &edge_part_bound_proc_idx, &edge_part_bound_part_idx, &edge_part_bound,
                               &vtx_tag, &_vtx, &_vtx_ln_to_gn, &edge_bound_idx, &edge_bound,
                               &edge_bound_ln_to_gn, &edge_join_idx, &edge_join, &edge_join_ln_to_gn);


    (*n_face)[i_part] = _n_face;
    (*n_vtx)[i_part]  = _n_vtx;

    (*face_vtx_idx)[i_part] = (int *) malloc (sizeof(int) * (_n_face + 1));
    memcpy ((*face_vtx_idx)[i_part], _face_edge_idx, sizeof(int) * (_n_face + 1));

    (*face_g_num)[i_part] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * _n_face);
    memcpy ((*face_g_num)[i_part], _face_ln_to_gn, sizeof(PDM_g_num_t) * _n_face);

    (*vtx_coord)[i_part] = (double *) malloc (sizeof(double) * _n_vtx * 3);
    memcpy ((*vtx_coord)[i_part], _vtx, sizeof(double) * _n_vtx * 3);

    (*vtx_g_num)[i_part] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * _n_vtx);
    memcpy ((*vtx_g_num)[i_part], _vtx_ln_to_gn, sizeof(PDM_g_num_t) * _n_vtx);


    /* Get face-vtx connectivity */
    _compute_face_vtx (_n_face,
                       _face_edge_idx,
                       _face_edge,
                       _edge_vtx,
                       *face_vtx + i_part);
  }

  PDM_multipart_free (mpart_id);
  free (dmesh);
}



static void
_create_split_mesh_multipart
(
 int               active_rank,
 PDM_MPI_Comm      comm,
 double            xmin,
 double            ymin,
 PDM_g_num_t       n_vtx_seg,
 double            length,
 double            depth,
 int               n_part,
 PDM_part_split_t  method,
 int               have_random,
 int               init_random,
 PDM_g_num_t      *ng_face,
 PDM_g_num_t      *ng_vtx,
 int             **n_face,
 int            ***face_vtx_idx,
 int            ***face_vtx,
 PDM_g_num_t    ***face_g_num,
 int             **n_vtx,
 double         ***vtx_coord,
 PDM_g_num_t    ***vtx_g_num
 )
{
  struct timeval t_elaps_debut;

  int i_rank;
  int n_rank;

  if (active_rank) {
    PDM_MPI_Comm_rank (comm, &i_rank);
    PDM_MPI_Comm_size (comm, &n_rank);

    double xmax = xmin + length;
    double ymax = ymin + length;
    PDM_g_num_t nx = n_vtx_seg;
    PDM_g_num_t ny = n_vtx_seg;

    int dn_face;
    int dn_edge;
    int dn_vtx;
    int n_edge_group;

    double      *dvtx_coord      = NULL;
    int         *dface_vtx_idx   = NULL;
    PDM_g_num_t *dface_vtx       = NULL;
    PDM_g_num_t *dface_edge      = NULL;
    PDM_g_num_t *dedge_face      = NULL;
    PDM_g_num_t *dedge_vtx       = NULL;
    int         *dedge_group_idx = NULL;
    PDM_g_num_t *dedge_group     = NULL;

    PDM_g_num_t ng_edge;

    /*
     *  Create mesh i
     */

    gettimeofday (&t_elaps_debut, NULL);

    PDM_poly_surf_gen (comm,
                       xmin,
                       xmax,
                       ymin,
                       ymax,
                       have_random,
                       init_random,
                       nx,
                       ny,
                       ng_face,
                       ng_vtx,
                       &ng_edge,
                       &dn_vtx,
                       &dvtx_coord,
                       &dn_face,
                       &dface_vtx_idx,
                       &dface_vtx,
                       &dface_edge,
                       &dn_edge,
                       &dedge_vtx,
                       &dedge_face,
                       &n_edge_group,
                       &dedge_group_idx,
                       &dedge_group);

    struct timeval t_elaps_fin;

    gettimeofday (&t_elaps_fin, NULL);

    long tranche_elapsed = (t_elaps_fin.tv_usec + 1000000 * t_elaps_fin.tv_sec) -
      (t_elaps_debut.tv_usec + 1000000 *
       t_elaps_debut.tv_sec);
    long tranche_elapsed_max = tranche_elapsed;
    double t_elapsed = (double) tranche_elapsed_max/1000000.;
    if (i_rank == 0)
      PDM_printf("[%d] Temps dans creeMaillagePolygone2D : %12.5e\n",
                 i_rank, t_elapsed);

    _add_depth (xmin,
                xmax,
                depth,
                dn_vtx,
                dvtx_coord);

    _split_multipart (comm,
                      method,
                      n_part,
                      dn_face,
                      dn_edge,
                      dn_vtx,
                      n_edge_group,
                      dedge_face,
                      dedge_vtx,
                      dvtx_coord,
                      dedge_group_idx,
                      dedge_group,
                      n_face,
                      face_vtx_idx,
                      face_vtx,
                      face_g_num,
                      n_vtx,
                      vtx_coord,
                      vtx_g_num);

    free (dvtx_coord);
    free (dface_vtx_idx);
    free (dface_vtx);
    free (dface_edge);
    free (dedge_face);
    free (dedge_vtx);
    free (dedge_group_idx);
    free (dedge_group);
  }

  else {
    *n_face       = (int *)          malloc (sizeof(int)           * n_part);
    *face_vtx_idx = (int **)         malloc (sizeof(int *)         * n_part);
    *face_vtx     = (int **)         malloc (sizeof(int *)         * n_part);
    *face_g_num   = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part);

    *n_vtx     = (int *)          malloc (sizeof(int)           * n_part);
    *vtx_coord = (double **)      malloc (sizeof(double *)      * n_part);
    *vtx_g_num = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {
      int _n_face = 0;
      int _sface_edge = 0;
      (*n_face)[ipart] = _n_face;
      (*face_vtx_idx)[ipart] = (int *) malloc (sizeof(int) * (_n_face + 1));
      (*face_vtx_idx)[ipart][0] = 0;
      (*face_vtx)[ipart]   = (int *)         malloc (sizeof(int)         * _sface_edge);
      (*face_g_num)[ipart] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * _n_face);

      int _n_vtx = 0;
      (*n_vtx)[ipart] = _n_vtx;
      (*vtx_coord)[ipart] = (double *)      malloc (sizeof(double)      * _n_vtx * 3);
      (*vtx_g_num)[ipart] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * _n_vtx);
    }
  }

  PDM_MPI_Bcast (ng_face, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
  PDM_MPI_Bcast (ng_vtx,  1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);

}

/**
 *
 * \brief  Create and split Mesh
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 *
 */

static void
_export_ini_mesh
(
 const PDM_MPI_Comm pdm_mpi_comm,
 const int      n_part,
 int            *nFace[2],
 int            **faceVtxIdx[2],
 int            **faceVtx[2],
 PDM_g_num_t    **faceLNToGN[2],
 int            *nVtx[2],
 double         **vtxCoord[2],
 PDM_g_num_t    **vtxLNToGN[2],
 double         **sFieldA,
 double         **rFieldB
 )
{

  int i_rank;
  int numProcs;

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  /*
   *  Export Mesh to Ensight
   */

  int id_cs[2];

  id_cs[0] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CONSTANTE,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "mesh1",
                                pdm_mpi_comm,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1.,
                                NULL);

  id_cs[1] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CONSTANTE,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "mesh2",
                                pdm_mpi_comm,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1,
                                NULL);

  /*
   * Creation des variables
   */

  int id_var_num_part[2];
  int id_var_field[2];
  int id_var_coo_x[2];
  int id_var_coo_xyz[2];
  int id_geom[2];

  for (int imesh = 0; imesh < 2; imesh++) {

    id_var_num_part[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                    PDM_WRITER_OFF,
                                                    PDM_WRITER_VAR_SCALAIRE,
                                                    PDM_WRITER_VAR_ELEMENTS,
                                                    "num_part");
    if (imesh == 0) {

      id_var_field[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                   PDM_WRITER_OFF,
                                                   PDM_WRITER_VAR_SCALAIRE,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "sfieldA");
    }
    else {

      id_var_field[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                   PDM_WRITER_OFF,
                                                   PDM_WRITER_VAR_SCALAIRE,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "rfieldB");
    }


    id_var_coo_x[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                 PDM_WRITER_ON,
                                                 PDM_WRITER_VAR_SCALAIRE,
                                                 PDM_WRITER_VAR_SOMMETS,
                                                 "coo_x");

    id_var_coo_xyz[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                   PDM_WRITER_ON,
                                                   PDM_WRITER_VAR_VECTEUR,
                                                   PDM_WRITER_VAR_SOMMETS,
                                                   "coo_xyz");

    /*
     * Creation de la geometrie
     */

    char nom_geom[6];
    if (imesh == 0)
      strcpy (nom_geom,"mesh1");
    else
      strcpy (nom_geom,"mesh2");

    id_geom[imesh] = PDM_writer_geom_create (id_cs[imesh],
                                             nom_geom,
                                             PDM_WRITER_OFF,
                                             PDM_WRITER_OFF,
                                             n_part);
    /*
     * Debut des ecritures
     */

    int *nsom_part  = (int *) malloc(sizeof(int) * n_part);

    int *n_part_procs = (int *) malloc(sizeof(int) * numProcs);

    PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
                       (void *) n_part_procs, 1, PDM_MPI_INT,
                       PDM_MPI_COMM_WORLD);

    int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

    debPartProcs[0] = 0;
    for (int i = 0; i < numProcs; i++) {
      debPartProcs[i+1] = debPartProcs[i] + n_part_procs[i];
    }

    free(n_part_procs);

    PDM_writer_step_beg (id_cs[imesh], 0.);

    int **_face_nb =  malloc(sizeof(int *) * n_part);
    int **_face_idx =  malloc(sizeof(int *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      PDM_writer_geom_coord_set (id_cs[imesh],
                                 id_geom[imesh],
                                 ipart,
                                 nVtx[imesh][ipart],
                                 vtxCoord[imesh][ipart],
                                 vtxLNToGN[imesh][ipart]);

      _face_nb[ipart] = malloc(sizeof(int) * nFace[imesh][ipart]);
      _face_idx[ipart] = malloc(sizeof(int) * nFace[imesh][ipart]);

      for (int j = 0; j < nFace[imesh][ipart]; j++) {
        _face_nb[ipart][j]  = faceVtxIdx[imesh][ipart][j+1] - faceVtxIdx[imesh][ipart][j];
        _face_idx[ipart][j] = faceVtxIdx[imesh][ipart][j] + 1;
      }

      PDM_writer_geom_faces_facesom_add (id_cs[imesh],
                                         id_geom[imesh],
                                         ipart,
                                         nFace[imesh][ipart],
                                         _face_idx[ipart],
                                         _face_nb[ipart],
                                         faceVtx[imesh][ipart],
                                         faceLNToGN[imesh][ipart]);
    }

    for (int ipart = 0; ipart < n_part; ipart++) {
      free (_face_nb[ipart]);
      free (_face_idx[ipart]);
    }

    free(_face_nb);
    free(_face_idx);

    PDM_writer_geom_write(id_cs[imesh],
                          id_geom[imesh]);

    /* Creation des variables :
       - numero de partition
       - scalaire
       - vecteur
       - tenseur
    */

    PDM_real_t **val_num_part = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_coo_x    = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_coo_xyz  = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      val_num_part[ipart] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nFace[imesh][ipart]);
      val_coo_x[ipart]    = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nVtx[imesh][ipart]);
      val_coo_xyz[ipart]  = (PDM_real_t *) malloc (sizeof(PDM_real_t) * 3 * nVtx[imesh][ipart]);
      nsom_part[ipart]    = nVtx[imesh][ipart];

      for (int i = 0; i < nFace[imesh][ipart]; i++) {
        val_num_part[ipart][i] = ipart + 1 + debPartProcs[i_rank];
      }

      for (int i = 0; i < nVtx[imesh][ipart]; i++) {
        val_coo_x[ipart][i]       = vtxCoord[imesh][ipart][3*i];
        val_coo_xyz[ipart][3*i  ] = vtxCoord[imesh][ipart][3*i  ];
        val_coo_xyz[ipart][3*i+1] = vtxCoord[imesh][ipart][3*i+1];
        val_coo_xyz[ipart][3*i+2] = vtxCoord[imesh][ipart][3*i+2];
      }

      if (imesh == 0) {
        PDM_writer_var_set (id_cs[imesh],
                            id_var_field[imesh],
                            id_geom[imesh],
                            ipart,
                            sFieldA[ipart]);
      }
      else {
        PDM_writer_var_set (id_cs[imesh],
                            id_var_field[imesh],
                            id_geom[imesh],
                            ipart,
                            rFieldB[ipart]);
      }


      PDM_writer_var_set (id_cs[imesh],
                          id_var_num_part[imesh],
                          id_geom[imesh],
                          ipart,
                          val_num_part[ipart]);

      PDM_writer_var_set (id_cs[imesh],
                          id_var_coo_x[imesh],
                          id_geom[imesh],
                          ipart,
                          val_coo_x[ipart]);

      PDM_writer_var_set (id_cs[imesh],
                          id_var_coo_xyz[imesh],
                          id_geom[imesh],
                          ipart,
                          val_coo_xyz[ipart]);

    }

    PDM_writer_var_write (id_cs[imesh],
                          id_var_field[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_field[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                          id_var_num_part[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_num_part[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                          id_var_coo_x[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_coo_x[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                          id_var_coo_xyz[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_coo_xyz[imesh]);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free (val_num_part[ipart]);
      free (val_coo_x[ipart]);
      free (val_coo_xyz[ipart]);
    }

    free (val_num_part);
    free (val_coo_x);
    free (val_coo_xyz);
    free (nsom_part);

    PDM_writer_step_end (id_cs[imesh]);
    PDM_writer_geom_data_free (id_cs[imesh],
                               id_geom[imesh]);

    PDM_writer_geom_free (id_cs[imesh],
                          id_geom[imesh]);
    PDM_writer_free (id_cs[imesh]);

    free (debPartProcs);

  }

}

/**
 *
 * \brief  Export overlay mesh
 *
 * \param [in]    pdm_id    PDM identifier
 *
 */

static void
_export_ol_mesh
(
 const PDM_MPI_Comm pdm_mpi_comm,
 const int pdm_id,
 int**     nFace,
 double**  sFieldOlA,
 double**  rFieldOlB,
 const int n_part
 )
{

  int i_rank;
  int numProcs;

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  /*
   *  Export Mesh to Ensight
   */

  int id_cs[2];

  id_cs[0] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CONSTANTE,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "olmesh1",
                                pdm_mpi_comm,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1.,
                                NULL);

  id_cs[1] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CONSTANTE,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "olmesh2",
                                pdm_mpi_comm,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1,
                                NULL);

  /*
   * Creation des variables
   */

  int id_var_field[2];
  int id_var_num_part[2];
  int id_var_match[2];
  int id_var_cell_match[2];
  int id_var_origin[2];
  int id_geom[2];

  for (int imesh = 0; imesh < 2; imesh++) {
    if (imesh == 0) {

      id_var_field[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                   PDM_WRITER_OFF,
                                                   PDM_WRITER_VAR_SCALAIRE,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "sOlField");
    }
    else {

      id_var_field[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                   PDM_WRITER_OFF,
                                                   PDM_WRITER_VAR_SCALAIRE,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "rOlField");
    }

    id_var_num_part[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                    PDM_WRITER_OFF,
                                                    PDM_WRITER_VAR_SCALAIRE,
                                                    PDM_WRITER_VAR_ELEMENTS,
                                                    "num_part");

    id_var_match[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                 PDM_WRITER_OFF,
                                                 PDM_WRITER_VAR_SCALAIRE,
                                                 PDM_WRITER_VAR_ELEMENTS,
                                                 "matching");

    id_var_cell_match[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                      PDM_WRITER_OFF,
                                                      PDM_WRITER_VAR_SCALAIRE,
                                                      PDM_WRITER_VAR_ELEMENTS,
                                                      "cell_matching");

    id_var_origin[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                  PDM_WRITER_OFF,
                                                  PDM_WRITER_VAR_SCALAIRE,
                                                  PDM_WRITER_VAR_ELEMENTS,
                                                  "origin");

    /*
     * Creation de la geometrie
     */

    char nom_geom[8];
    PDM_ol_mesh_t mesht;

    if (imesh == 0) {
      strcpy (nom_geom,"olmesh1");
      mesht = PDM_OL_MESH_A;
    }
    else {
      strcpy (nom_geom,"olmesh2");
      mesht = PDM_OL_MESH_B;
    }

    id_geom[imesh] = PDM_writer_geom_create (id_cs[imesh],
                                             nom_geom,
                                             PDM_WRITER_OFF,
                                             PDM_WRITER_OFF,
                                             n_part);
    int *n_part_procs = (int *) malloc(sizeof(int) * numProcs);

    PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
                       (void *) n_part_procs, 1, PDM_MPI_INT,
                       PDM_MPI_COMM_WORLD);

    int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

    debPartProcs[0] = 0;
    for (int i = 0; i < numProcs; i++) {
      debPartProcs[i+1] = debPartProcs[i] + n_part_procs[i];
    }

    free(n_part_procs);

    /*
     * Debut des ecritures
     */

    PDM_g_num_t    nGOlFace;
    PDM_g_num_t    nGOlVtx;

    PDM_ol_mesh_dim_get (pdm_id,
                         mesht,
                         &nGOlFace,
                         &nGOlVtx);

    int **_olface_nb =  malloc(sizeof(int *) * n_part);
    int **_olface_idx =  malloc(sizeof(int *) * n_part);
    PDM_real_t **val_num_part = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_match = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_cell_match = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_origin = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_writer_step_beg (id_cs[imesh], 0.);

    for (int ipart = 0; ipart < n_part; ipart++) {

      int           nOlFace;
      int           nOlLinkedFace;
      int           nOlVtx;
      int           sOlFaceIniVtx;
      int           sOlface_vtx;
      int           sInitToOlFace;

      PDM_ol_part_mesh_dim_get (pdm_id,
                                mesht,
                                ipart,
                                &nOlFace,
                                &nOlLinkedFace,
                                &nOlVtx,
                                &sOlFaceIniVtx,
                                &sOlface_vtx,
                                &sInitToOlFace);

      int            *olFaceIniVtxIdx;
      int            *olFaceIniVtx;
      int            *olface_vtx_idx;
      int            *olface_vtx;
      int            *olLinkedface_procIdx;
      int            *olLinkedFace;
      PDM_g_num_t     *olface_ln_to_gn;
      double         *olCoords;
      PDM_g_num_t     *olvtx_ln_to_gn;
      int            *initToOlFaceIdx;
      int            *initToOlFace;

      PDM_ol_mesh_entities_get (pdm_id,
                                mesht,
                                ipart,
                                &olFaceIniVtxIdx,
                                &olFaceIniVtx,
                                &olface_vtx_idx,
                                &olface_vtx,
                                &olLinkedface_procIdx,
                                &olLinkedFace,
                                &olface_ln_to_gn,
                                &olCoords,
                                &olvtx_ln_to_gn,
                                &initToOlFaceIdx,
                                &initToOlFace);


      val_num_part[ipart] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nOlFace);
      val_match[ipart] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nOlFace);
      val_cell_match[ipart] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nOlFace);
      val_origin[ipart] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nOlFace);

      PDM_writer_geom_coord_set (id_cs[imesh],
                                 id_geom[imesh],
                                 ipart,
                                 nOlVtx,
                                 olCoords,
                                 olvtx_ln_to_gn);

      _olface_nb[ipart] = malloc(sizeof(int) * nOlFace);
      _olface_idx[ipart] = malloc(sizeof(int) * nOlFace);

      for (int j = 0; j < nOlFace; j++) {
        _olface_nb[ipart][j]  = olface_vtx_idx[j+1] - olface_vtx_idx[j];
        _olface_idx[ipart][j] = olface_vtx_idx[j] + 1;
      }

      PDM_writer_geom_faces_facesom_add (id_cs[imesh],
                                         id_geom[imesh],
                                         ipart,
                                         nOlFace,
                                         _olface_idx[ipart],
                                         _olface_nb[ipart],
                                         olface_vtx,
                                         olface_ln_to_gn);

      for (int i = 0; i < nOlFace; i++) {
        val_num_part[ipart][i] = ipart + 1 + debPartProcs[i_rank];
        val_origin[ipart][i] = -1;
        val_match[ipart][i] = -1;
        val_cell_match[ipart][i] = -1;
      }

      for (int i = 0; i < nOlLinkedFace; i++) {
        val_match[ipart][olLinkedFace[4*i]-1] = 100;
        val_cell_match[ipart][olLinkedFace[4*i]-1] = olLinkedFace[4*i + 2];
      }

      for (int i = 0; i < nFace[imesh][ipart]; i++) {
        for (int j = initToOlFaceIdx[i]; j < initToOlFaceIdx[i+1]; j++) {
          val_origin[ipart][initToOlFace[j]-1] = i;
        }
      }

      if (imesh == 0) {
        PDM_writer_var_set (id_cs[imesh],
                            id_var_field[imesh],
                            id_geom[imesh],
                            ipart,
                            sFieldOlA[ipart]);
      }
      else {
        PDM_writer_var_set (id_cs[imesh],
                            id_var_field[imesh],
                            id_geom[imesh],
                            ipart,
                            rFieldOlB[ipart]);
      }

      PDM_writer_var_set (id_cs[imesh],
                          id_var_num_part[imesh],
                          id_geom[imesh],
                          ipart,
                          val_num_part[ipart]);

      PDM_writer_var_set (id_cs[imesh],
                          id_var_match[imesh],
                          id_geom[imesh],
                          ipart,
                          val_match[ipart]);

      PDM_writer_var_set (id_cs[imesh],
                          id_var_cell_match[imesh],
                          id_geom[imesh],
                          ipart,
                          val_cell_match[ipart]);

      PDM_writer_var_set (id_cs[imesh],
                          id_var_origin[imesh],
                          id_geom[imesh],
                          ipart,
                          val_origin[ipart]);


    }

    PDM_writer_geom_write(id_cs[imesh],
                          id_geom[imesh]);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free (_olface_nb[ipart]);
      free (_olface_idx[ipart]);
    }

    PDM_writer_var_write (id_cs[imesh],
                          id_var_num_part[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                          id_var_match[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                          id_var_cell_match[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                          id_var_origin[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                          id_var_field[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_num_part[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_match[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_cell_match[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_origin[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_field[imesh]);

    for (int ipart = 0; ipart < n_part; ipart++) {
      free (val_num_part[ipart]);
      free (val_match[ipart]);
      free (val_cell_match[ipart]);
      free (val_origin[ipart]);
    }

    free (val_num_part);
    free (val_match);
    free (val_cell_match);
    free (val_origin);

    PDM_writer_step_end (id_cs[imesh]);
    PDM_writer_geom_data_free (id_cs[imesh],
                               id_geom[imesh]);

    PDM_writer_geom_free (id_cs[imesh],
                          id_geom[imesh]);
    PDM_writer_free (id_cs[imesh]);

    free (_olface_nb);
    free (_olface_idx);
    free (debPartProcs);
  }

}


static void _export_ensight
(
 const PDM_MPI_Comm   pdm_mpi_comm,
 const int            n_part,
 int                 *nFace[2],
 int                **faceVtxIdx[2],
 int                **faceVtx[2],
 PDM_g_num_t        **faceLNToGN[2],
 int                 *nVtx[2],
 double             **vtxCoord[2],
 PDM_g_num_t        **vtxLNToGN[2]
 )
{
  int i_rank;
  int numProcs;

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  /*
   *  Export Mesh to Ensight
   */

  int id_cs[2];

  id_cs[0] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CONSTANTE,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "meshA",
                                pdm_mpi_comm,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1.,
                                NULL);

  id_cs[1] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CONSTANTE,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "meshB",
                                pdm_mpi_comm,
                                PDM_IO_ACCES_MPI_SIMPLE,
                                1,
                                NULL);

  /*
   * Creation des variables
   */

  int id_var_num_part[2];
  int id_var_face_gnum[2];
  int id_var_vtx_gnum[2];
  int id_geom[2];

  for (int imesh = 0; imesh < 2; imesh++) {

    id_var_num_part[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                    PDM_WRITER_OFF,
                                                    PDM_WRITER_VAR_SCALAIRE,
                                                    PDM_WRITER_VAR_ELEMENTS,
                                                    "num_part");

    id_var_face_gnum[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                    PDM_WRITER_OFF,
                                                    PDM_WRITER_VAR_SCALAIRE,
                                                    PDM_WRITER_VAR_ELEMENTS,
                                                    "face_gnum");

    id_var_vtx_gnum[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                    PDM_WRITER_OFF,
                                                    PDM_WRITER_VAR_SCALAIRE,
                                                    PDM_WRITER_VAR_SOMMETS,
                                                    "vtx_gnum");

    /*
     * Creation de la geometrie
     */

    char nom_geom[6];
    if (imesh == 0)
      strcpy (nom_geom,"mesh1");
    else
      strcpy (nom_geom,"mesh2");

    id_geom[imesh] = PDM_writer_geom_create (id_cs[imesh],
                                             nom_geom,
                                             PDM_WRITER_OFF,
                                             PDM_WRITER_OFF,
                                             n_part);

    /*
     * Debut des ecritures
     */

    int *nsom_part  = (int *) malloc(sizeof(int) * n_part);

    int *n_part_procs = (int *) malloc(sizeof(int) * numProcs);

    PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
                       (void *) n_part_procs, 1, PDM_MPI_INT,
                       PDM_MPI_COMM_WORLD);

    int *debPartProcs = (int *) malloc(sizeof(int) * (numProcs + 1));

    debPartProcs[0] = 0;
    for (int i = 0; i < numProcs; i++) {
      debPartProcs[i+1] = debPartProcs[i] + n_part_procs[i];
    }

    free(n_part_procs);

    PDM_writer_step_beg (id_cs[imesh], 0.);

    int **_face_nb =  malloc(sizeof(int *) * n_part);
    int **_face_idx =  malloc(sizeof(int *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      PDM_writer_geom_coord_set (id_cs[imesh],
                                 id_geom[imesh],
                                 ipart,
                                 nVtx[imesh][ipart],
                                 vtxCoord[imesh][ipart],
                                 vtxLNToGN[imesh][ipart]);

      _face_nb[ipart] = malloc(sizeof(int) * nFace[imesh][ipart]);
      _face_idx[ipart] = malloc(sizeof(int) * nFace[imesh][ipart]);

      for (int j = 0; j < nFace[imesh][ipart]; j++) {
        _face_nb[ipart][j]  = faceVtxIdx[imesh][ipart][j+1] - faceVtxIdx[imesh][ipart][j];
        _face_idx[ipart][j] = faceVtxIdx[imesh][ipart][j] + 1;
      }

      PDM_writer_geom_faces_facesom_add (id_cs[imesh],
                                         id_geom[imesh],
                                         ipart,
                                         nFace[imesh][ipart],
                                         _face_idx[ipart],
                                         _face_nb[ipart],
                                         faceVtx[imesh][ipart],
                                         faceLNToGN[imesh][ipart]);
    } // End loop on parts

    for (int ipart = 0; ipart < n_part; ipart++) {
      free (_face_nb[ipart]);
      free (_face_idx[ipart]);
    }

    free(_face_nb);
    free(_face_idx);

    PDM_writer_geom_write(id_cs[imesh],
                          id_geom[imesh]);

    /* Creation des variables :
       - numero de partition
       - normale aux sommets
    */
    PDM_real_t **val_num_part  = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_face_gnum = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_vtx_gnum  = (PDM_real_t **) malloc (sizeof(PDM_real_t *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {

      val_num_part[ipart]  = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nFace[imesh][ipart]);
      val_face_gnum[ipart] = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nFace[imesh][ipart]);
      val_vtx_gnum[ipart]  = (PDM_real_t *) malloc (sizeof(PDM_real_t) * nVtx[imesh][ipart]);
      nsom_part[ipart] = nVtx[imesh][ipart];

      for (int i = 0; i < nFace[imesh][ipart]; i++) {
        val_num_part[ipart][i] = ipart + 1 + debPartProcs[i_rank];
        val_face_gnum[ipart][i] = faceLNToGN[imesh][ipart][i];
      }

      for (int i = 0; i < nVtx[imesh][ipart]; i++) {
        val_vtx_gnum[ipart][i] = vtxLNToGN[imesh][ipart][i];
      }

      PDM_writer_var_set (id_cs[imesh],
                          id_var_num_part[imesh],
                          id_geom[imesh],
                          ipart,
                          val_num_part[ipart]);

      PDM_writer_var_set (id_cs[imesh],
                          id_var_face_gnum[imesh],
                          id_geom[imesh],
                          ipart,
                          val_face_gnum[ipart]);

      PDM_writer_var_set (id_cs[imesh],
                          id_var_vtx_gnum[imesh],
                          id_geom[imesh],
                          ipart,
                          val_vtx_gnum[ipart]);
    } // End loop on parts

    PDM_writer_var_write (id_cs[imesh],
                          id_var_num_part[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_num_part[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                          id_var_face_gnum[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_face_gnum[imesh]);

    PDM_writer_var_write (id_cs[imesh],
                          id_var_vtx_gnum[imesh]);

    PDM_writer_var_free (id_cs[imesh],
                         id_var_vtx_gnum[imesh]);

   for (int ipart = 0; ipart < n_part; ipart++) {
      free (val_num_part[ipart]);
      free (val_face_gnum[ipart]);
      free (val_vtx_gnum[ipart]);
    }

    free (val_num_part);
    free (val_face_gnum);
    free (val_vtx_gnum);
    free (nsom_part);

    PDM_writer_step_end (id_cs[imesh]);
    PDM_writer_geom_data_free (id_cs[imesh],
                               id_geom[imesh]);

    PDM_writer_geom_free (id_cs[imesh],
                          id_geom[imesh]);
    PDM_writer_free (id_cs[imesh]);

    free (debPartProcs);
  } // End loop on meshes
}







static void _export_vtk
(
 const PDM_MPI_Comm  pdm_mpi_comm,
 const char         *name,
 const int           n_part,
 int                *nFace,
 int               **faceVtxIdx,
 int               **faceVtx,
 PDM_g_num_t       **faceLNToGN,
 double            **face_vector,
 int                *nVtx,
 double            **vtxCoord,
 PDM_g_num_t       **vtxLNToGN,
 double            **vtx_vector
 )
{
  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  for (int i_part = 0; i_part < n_part; i_part++) {
    char filename[999];
    sprintf(filename, "%s_%3.3d.vtk", name, n_part*i_rank + i_part);

    FILE *f = fopen(filename, "w");

    fprintf(f, "# vtk DataFile Version 2.0\nmesh\nASCII\nDATASET POLYDATA\n");

    fprintf(f, "POINTS %d double\n", nVtx[i_part]);
    for (int i = 0; i < nVtx[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        fprintf(f, "%lf ", vtxCoord[i_part][3*i+j]);
      }
      fprintf(f, "\n");
    }

    fprintf(f, "POLYGONS %d %d\n", nFace[i_part], nFace[i_part] + faceVtxIdx[i_part][nFace[i_part]]);
    for (int i = 0; i < nFace[i_part]; i++) {
      fprintf(f, "%d ", faceVtxIdx[i_part][i+1] - faceVtxIdx[i_part][i]);
      for (int j = faceVtxIdx[i_part][i]; j < faceVtxIdx[i_part][i+1]; j++) {
        fprintf(f, "%d ", faceVtx[i_part][j] - 1);
      }
      fprintf(f, "\n");
    }

    if (vtx_vector != NULL) {
      fprintf(f, "POINT_DATA %d\n", nVtx[i_part]);
      fprintf(f, "VECTORS pvec double\n");
      for (int i = 0; i < nVtx[i_part]; i++) {
        for (int j = 0; j < 3; j++) {
          fprintf(f, "%lf ", vtx_vector[i_part][3*i+j]);
        }
        fprintf(f, "\n");
      }
    }

    if (1) {//face_vector != NULL || faceLNToGN != NULL) {
      fprintf(f, "CELL_DATA %d\n", nFace[i_part]);
      if (0) {//faceLNToGN != NULL) {
        fprintf(f, "SCALARS gnum int\n");
        fprintf(f, "LOOKUP_TABLE default\n");
        for (int i = 0; i < nFace[i_part]; i++) {
          fprintf(f, PDM_FMT_G_NUM"\n", faceLNToGN[i_part][i]);
        }
      }
      else {
        fprintf(f, "SCALARS rank int\n");
        fprintf(f, "LOOKUP_TABLE default\n");
        for (int i = 0; i < nFace[i_part]; i++) {
          fprintf(f, "%d\n", i_rank);
        }
      }


      if (face_vector != NULL) {
        fprintf(f, "VECTORS fvec double\n");
        for (int i = 0; i < nFace[i_part]; i++) {
          for (int j = 0; j < 3; j++) {
            fprintf(f, "%lf ", face_vector[i_part][3*i+j]);
          }
          fprintf(f, "\n");
        }
      }

    }



    fclose(f);
  }
}




static void
_export_ol_vtk
(
 const PDM_MPI_Comm   pdm_mpi_comm,
 const int            pdm_id,
 int                **nFace,
 double             **sFieldOlA,
 double             **rFieldOlB,
 const int            n_part
 )
{
  int i_rank;
  int numProcs;

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  char nom_geom[8];
  for (int imesh = 0; imesh < 2; imesh++) {

    double **field;
    PDM_ol_mesh_t mesht;

    if (imesh == 0) {
      strcpy (nom_geom,"olmesh1");
      mesht = PDM_OL_MESH_A;
      field = sFieldOlA;
    }
    else {
      strcpy (nom_geom,"olmesh2");
      mesht = PDM_OL_MESH_B;
      field = rFieldOlB;
    }

    for (int ipart = 0; ipart < n_part; ipart++) {

      int           nOlFace;
      int           nOlLinkedFace;
      int           nOlVtx;
      int           sOlFaceIniVtx;
      int           sOlface_vtx;
      int           sInitToOlFace;

      PDM_ol_part_mesh_dim_get (pdm_id,
                                mesht,
                                ipart,
                                &nOlFace,
                                &nOlLinkedFace,
                                &nOlVtx,
                                &sOlFaceIniVtx,
                                &sOlface_vtx,
                                &sInitToOlFace);

      int            *olFaceIniVtxIdx;
      int            *olFaceIniVtx;
      int            *olface_vtx_idx;
      int            *olface_vtx;
      int            *olLinkedface_procIdx;
      int            *olLinkedFace;
      PDM_g_num_t    *olface_ln_to_gn;
      double         *olCoords;
      PDM_g_num_t    *olvtx_ln_to_gn;
      int            *initToOlFaceIdx;
      int            *initToOlFace;

      PDM_ol_mesh_entities_get (pdm_id,
                                mesht,
                                ipart,
                                &olFaceIniVtxIdx,
                                &olFaceIniVtx,
                                &olface_vtx_idx,
                                &olface_vtx,
                                &olLinkedface_procIdx,
                                &olLinkedFace,
                                &olface_ln_to_gn,
                                &olCoords,
                                &olvtx_ln_to_gn,
                                &initToOlFaceIdx,
                                &initToOlFace);

      int n_tri = 0;
      int max_n_vtx_poly = 0;
      for (int ipoly = 0; ipoly < nOlFace; ipoly++) {
        int n_vtx_poly = olface_vtx_idx[ipoly+1] - olface_vtx_idx[ipoly];
        n_tri += (n_vtx_poly - 2)*3;
        max_n_vtx_poly = PDM_MAX (max_n_vtx_poly, n_vtx_poly);
      }

      PDM_triangulate_state_t *state = PDM_triangulate_state_create (max_n_vtx_poly);

      int *tri_vtx_idx = malloc (sizeof(int) * (n_tri+1));
      int *tri_vtx = malloc (sizeof(int) * n_tri * 3);
      int *tri_poly = malloc (sizeof(int) * n_tri);
      tri_vtx_idx[0] = 0;
      n_tri = 0;
      for (int ipoly = 0; ipoly < nOlFace; ipoly++) {
        int n_vtx_poly = olface_vtx_idx[ipoly+1] - olface_vtx_idx[ipoly];

        int _n_tri = PDM_triangulate_polygon (3,
                                              n_vtx_poly,
                                              olCoords,
                                              NULL,
                                              olface_vtx + olface_vtx_idx[ipoly],
                                              PDM_TRIANGULATE_MESH_DEF,
                                              tri_vtx + tri_vtx_idx[n_tri],
                                              state);

        for (int i = 0; i < _n_tri; i++) {
          tri_vtx_idx[n_tri+1] = tri_vtx_idx[n_tri] + 3;
          tri_poly[n_tri] = ipoly + 1;
          n_tri++;
        }
      } // End loop on polygons


      char filename[999];
      sprintf(filename, "%s_%3.3d.vtk", nom_geom, n_part*i_rank + ipart);

      FILE *f = fopen(filename, "w");

      fprintf(f, "# vtk DataFile Version 2.0\nmesh\nASCII\nDATASET POLYDATA\n");

      fprintf(f, "POINTS %d double\n", nOlVtx);
      for (int i = 0; i < nOlVtx; i++) {
        for (int j = 0; j < 3; j++) {
          fprintf(f, "%lf ", olCoords[3*i+j]);
        }
        fprintf(f, "\n");
      }

      if (1) {
        fprintf(f, "POLYGONS %d %d\n", n_tri, 4*n_tri);
        for (int i = 0; i < n_tri; i++) {
          fprintf(f, "%d ", tri_vtx_idx[i+1] - tri_vtx_idx[i]);
          for (int j = tri_vtx_idx[i]; j < tri_vtx_idx[i+1]; j++) {
            fprintf(f, "%d ", tri_vtx[j] - 1);
          }
          fprintf(f, "\n");
        }

        fprintf(f, "CELL_DATA %d\n", n_tri);
        fprintf(f, "SCALARS field double 1\n");
        fprintf(f, "LOOKUP_TABLE default\n");
        for (int i = 0; i < n_tri; i++) {
          int iface = tri_poly[i] - 1;
          fprintf(f, "%f\n", field[ipart][iface]);
        }

      } else {
        fprintf(f, "POLYGONS %d %d\n", nOlFace, nOlFace + olface_vtx_idx[nOlFace]);
        for (int i = 0; i < nOlFace; i++) {
          fprintf(f, "%d ", olface_vtx_idx[i+1] - olface_vtx_idx[i]);
          for (int j = olface_vtx_idx[i]; j < olface_vtx_idx[i+1]; j++) {
            fprintf(f, "%d ", olface_vtx[j] - 1);
          }
          fprintf(f, "\n");
        }
      }

      fclose(f);


      state = PDM_triangulate_state_destroy (state);
      free (tri_vtx_idx);
      free (tri_vtx);
      free (tri_poly);
    } // End loop on parts

  } // End loop on meshes
}


/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{
  PDM_MPI_Init (&argc, &argv);

  /*
   *  Set default values
   */

  PDM_g_num_t      n_vtx_segA = 4;
  double           lengthA  = 1.;
  double           xminA = 0.;
  double           yminA = 0.;
  int              n_partA   = 1;

  PDM_g_num_t      n_vtx_segB = 8;
  double           lengthB  = 1.;
  double           xminB = 0.;
  double           yminB = 0.;
  int              n_partB   = 1;

  double depth = 0.5;

  int              post    = 0;
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
  int              haveRandom = 1;
  int              randomTimeInit = 0;

  int              randomMeshAInit = -1;
  int              randomMeshBInit = -1;

  int              nProcData = -1;

  int              i_rank;
  int              numProcs;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_segA,
              &lengthA,
              &xminA,
              &yminA,
              &n_partA,
              &n_vtx_segB,
              &lengthB,
              &xminB,
              &yminB,
              &n_partB,
              &depth,
              &post,
              (int *) &method,
              &haveRandom,
              &randomTimeInit,
              &randomMeshAInit,
              &randomMeshBInit,
              &nProcData);

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_vtx_segA : %d\n", n_vtx_segA);
    PDM_printf ("  - lengthA : %f\n", lengthA);
    PDM_printf ("  - xminA : %d\n", xminA);
    PDM_printf ("  - yminA : %d\n", yminA);
    PDM_printf ("  - n_partA : %d\n", n_partA);
    PDM_printf ("  - n_vtx_segB : %d\n", n_vtx_segB);
    PDM_printf ("  - lengthB : %f\n", lengthB);
    PDM_printf ("  - xminB : %d\n", xminB);
    PDM_printf ("  - yminB : %d\n", yminB);
    PDM_printf ("  - n_partB : %d\n", n_partB);
    PDM_printf ("  - post : %d\n", post);
    PDM_printf ("  - method : %d\n", method);
    PDM_printf ("  - haveRandom : %d\n", haveRandom);
    PDM_printf ("  - randomTimeInit : %d\n", randomTimeInit);
  }

  /*
   *  Create a partitioned mesh
   */

  PDM_g_num_t      nGFace[2];
  PDM_g_num_t      nGVtx[2];
  int            *nFace[2];
  int            **faceVtxIdx[2];
  int            **faceVtx[2];
  PDM_g_num_t    **faceLNToGN[2];
  int            *nVtx[2];
  double         **vtxCoord[2];
  PDM_g_num_t    **vtxLNToGN[2];

  int initRandom = 0;

  assert (n_partA == n_partB);
  int n_part = n_partA;

  PDM_MPI_Comm meshComm = PDM_MPI_COMM_WORLD;
  int activeRankMesh = 1;

  if (nProcData > 0 && nProcData < numProcs) {
    int rankInNode = PDM_io_mpi_node_rank (PDM_MPI_COMM_WORLD);

    int nNode = 0;
    int iNode = -1;
    int masterRank = 0;
    if (rankInNode == 0) {
      masterRank = 1;
    }

    int *rankInNodes = malloc(sizeof(int) * numProcs);

    PDM_MPI_Allreduce (&masterRank, &nNode, 1, PDM_MPI_INT, PDM_MPI_SUM, PDM_MPI_COMM_WORLD);
    PDM_MPI_Allgather (&rankInNode, 1, PDM_MPI_INT, rankInNodes, 1, PDM_MPI_INT, PDM_MPI_COMM_WORLD);

    activeRankMesh = 0;

    for (int i = 0; i < i_rank; i++) {
      if (rankInNodes[i] == 0) {
        iNode += 1;
      }
    }

    if (nProcData <= nNode) {
      if (iNode < nProcData && rankInNode == 0) {
        activeRankMesh = 1;
      }
    }

    else {

      if (rankInNode < (nProcData / nNode)) {
        activeRankMesh = 1;
      }
      if ((rankInNode == (nProcData / nNode)) && (iNode < (nProcData % nNode))) {
        activeRankMesh = 1;
      }

    }

    PDM_MPI_Comm_split(PDM_MPI_COMM_WORLD, activeRankMesh, i_rank, &meshComm);

    free (rankInNodes);
  }

  for (int imesh = 0; imesh < 2; imesh++) {

    double xmin;
    double ymin;
    double length;
    PDM_g_num_t n_vtx_seg;

    if (imesh == 0) {
      n_vtx_seg = n_vtx_segA;
      length = lengthA;
      xmin = xminA;
      ymin = yminA;
      n_part = n_partA;
    }
    else {
      n_vtx_seg = n_vtx_segB;
      length = lengthB;
      xmin = xminB;
      ymin = yminB;
      n_part = n_partB;
    }

    if (randomTimeInit) {
      initRandom =  time( NULL );
    }

    if (imesh == 0 && randomMeshAInit != -1) {
      initRandom = randomMeshAInit;
    }

    if (imesh == 1 && randomMeshBInit != -1) {
      initRandom = randomMeshBInit;
    }

    _create_split_mesh_multipart (activeRankMesh,
                                  meshComm,
                                  xmin,
                                  ymin,
                                  n_vtx_seg,
                                  length,
                                  depth,
                                  n_part,
                                  method,
                                  haveRandom,
                                  initRandom,
                                  nGFace + imesh,
                                  nGVtx + imesh,
                                  nFace + imesh,
                                  faceVtxIdx + imesh,
                                  faceVtx + imesh,
                                  faceLNToGN  + imesh,
                                  nVtx + imesh,
                                  vtxCoord + imesh,
                                  vtxLNToGN + imesh);

    ++initRandom;

  }

  if (1) {
    _export_ensight (PDM_MPI_COMM_WORLD,
                     n_part,
                     nFace,
                     faceVtxIdx,
                     faceVtx,
                     faceLNToGN,
                     nVtx,
                     vtxCoord,
                     vtxLNToGN);
  }

  if (nProcData > 0 && nProcData < numProcs) {
    PDM_MPI_Comm_free(&meshComm);
  }

  /*
   *  Appel des fonctions d'intersection
   */

  double projectCoeff = 0.;

  /*
   *  Creation de l'objet PDM
   */

  int pdm_id = PDM_ol_create (n_part,
                              nGFace[0],
                              nGVtx[0],
                              n_part,
                              nGFace[1],
                              nGVtx[1],
                              projectCoeff,
                              PDM_MPI_COMM_WORLD);
  if (i_rank == 0){
    PDM_printf ("- n_part MeshA : %d \n", n_part);
    PDM_printf ("- nGFaceMeshA : %d \n", nGFace[0]);
    PDM_printf ("- nGVtxMeshA : %d \n", nGVtx[0]);
    PDM_printf ("- n_part MeshB : %d \n", n_part);
    PDM_printf ("- nGFaceMeshB : %d \n", nGFace[1]);
    PDM_printf ("- nGVtxMeshB : %d \n", nGVtx[1]);
    PDM_printf ("- projectCoeff : %d \n", projectCoeff);
  }

  PDM_ol_parameter_set (pdm_id,
                        PDM_OL_CAR_LENGTH_TOL,
                        1e-4);

  PDM_ol_parameter_set (pdm_id,
                        PDM_OL_EXTENTS_TOL,
                        1e-4);

  /*
   *  Create field
   */

  double **sFieldA = malloc(sizeof(double *) * n_part);
  double **centerA = malloc(sizeof(double *) * n_part);
  double **rFieldB = malloc(sizeof(double *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    centerA[ipart] = malloc(sizeof(double) * 3 * nFace[0][ipart]);
    for (int i = 0; i < 3*nFace[0][ipart]; i++) {
      centerA[ipart][i] = 0.;
    }
    sFieldA[ipart] = malloc(sizeof(double) * nFace[0][ipart]);
    rFieldB[ipart] = malloc(sizeof(double) * nFace[1][ipart]);
  }

  /*
   *  Initial meshes definition
   */

  for (int imesh = 0; imesh < 2; imesh++) {

    PDM_ol_mesh_t mesh;
    if (imesh == 0) {
      mesh = PDM_OL_MESH_A;
    }
    else {
      mesh = PDM_OL_MESH_B;
    }

    for (int ipart = 0; ipart < n_part; ipart++) {

      PDM_ol_input_mesh_set (pdm_id,
                             mesh,
                             ipart,
                             nFace[imesh][ipart],
                             faceVtxIdx[imesh][ipart],
                             faceVtx[imesh][ipart],
                             faceLNToGN[imesh][ipart],
                             nVtx[imesh][ipart],
                             vtxCoord[imesh][ipart],
                             vtxLNToGN[imesh][ipart]);
      if (imesh == 0) {

        // Calcul du barycentre des sommets pour la face. Appeler les fonctions de calcul de centre pour mieux faire

        for (int i = 0; i < nFace[imesh][ipart]; i++) {
          for (int j = faceVtxIdx[imesh][ipart][i]; j < faceVtxIdx[imesh][ipart][i+1]; j++) {
            int    k = faceVtx[imesh][ipart][j] - 1;
            double x = vtxCoord[imesh][ipart][3*k  ];
            double y = vtxCoord[imesh][ipart][3*k+1];
            double z = vtxCoord[imesh][ipart][3*k+2];
            centerA[ipart][3*i  ] += x;
            centerA[ipart][3*i+1] += y;
            centerA[ipart][3*i+2] += z;
          }

          int nVtxElt = faceVtxIdx[imesh][ipart][i+1] - faceVtxIdx[imesh][ipart][i];
          centerA[ipart][3*i  ] /= nVtxElt;
          centerA[ipart][3*i+1] /= nVtxElt;
          centerA[ipart][3*i+2] /= nVtxElt;

          double _x = centerA[ipart][3*i]  - (lengthA -xminA)/2;
          double _y = centerA[ipart][3*i+1]- (lengthA -yminA)/2;

          double r = PDM_MIN ((lengthA -xminA)/2, (lengthA -yminA)/2)/4;
          double A = 10.;

          sFieldA[ipart][i] = A * exp((-(_x*_x)-(_y*_y))/(2*r*r));

        }
      }
    }
  }


#if 0
  if (post) {
    _export_ini_mesh (PDM_MPI_COMM_WORLD,
                      n_part,
                      nFace,
                      faceVtxIdx,
                      faceVtx,
                      faceLNToGN,
                      nVtx,
                      vtxCoord,
                      vtxLNToGN,
                      sFieldA,
                      rFieldB);
    _export_vtk (PDM_MPI_COMM_WORLD,
                 "meshA",
                 n_part,
                 nFace[0],
                 faceVtxIdx[0],
                 faceVtx[0],
                 faceLNToGN[0],
                 NULL,
                 nVtx[0],
                 vtxCoord[0],
                 vtxLNToGN[0],
                 NULL);
  }

#else
  /*
   *  Calcul
   */

  PDM_ol_compute (pdm_id);

  if (i_rank == 0){
    PDM_ol_dump_times (pdm_id);
  }

  /*
   *  Check graph
   */

  int  *_nOlFace[2];
  int  *_nOlLinkedFace[2];
  int  *_nOlVtx[2];
  int  *_sOlface_vtx[2];
  int  *_sInitToOlFace[2];

  int **_olFaceIniVtxIdx[2];
  int **_olFaceIniVtx[2];
  int **_olface_vtx_idx[2];
  int **_olface_vtx[2];
  int **_olLinkedface_procIdx[2];
  int **_olLinkedFace[2];
  PDM_g_num_t **_olface_ln_to_gn[2];
  double **_olCoords[2];
  PDM_g_num_t **_olvtx_ln_to_gn[2];
  int **_initToOlFaceIdx[2];
  int **_initToOlFace[2];

  double **sFieldOlA = malloc(sizeof(double *) * n_part);
  double **rFieldOlB = malloc(sizeof(double *) * n_part);
  double **surfB = malloc(sizeof(double *) * n_part);
  double **surfOlB = malloc(sizeof(double *) * n_part);

  for (int imesh = 0; imesh < 2; imesh++) {
    _nOlFace[imesh] = malloc (sizeof(int) * n_part);
    _nOlLinkedFace[imesh] = malloc (sizeof(int) * n_part);
    _nOlVtx[imesh] = malloc (sizeof(int) * n_part);
    _sOlface_vtx[imesh] = malloc (sizeof(int) * n_part);
    _sInitToOlFace[imesh] = malloc (sizeof(int) * n_part);
    _olFaceIniVtxIdx[imesh] = malloc (sizeof(int*) * n_part);
    _olFaceIniVtx[imesh] = malloc (sizeof(int*) * n_part);
    _olface_vtx_idx[imesh] = malloc (sizeof(int*) * n_part);
    _olface_vtx[imesh] = malloc (sizeof(int*) * n_part);
    _olLinkedface_procIdx[imesh] = malloc (sizeof(int*) * n_part);
    _olLinkedFace[imesh] = malloc (sizeof(int*) * n_part);
    _olface_ln_to_gn[imesh] = malloc (sizeof(PDM_g_num_t*) * n_part);
    _olCoords[imesh] = malloc (sizeof(double*) * n_part);;
    _olvtx_ln_to_gn[imesh] = malloc (sizeof(PDM_g_num_t*) * n_part);;
    _initToOlFaceIdx[imesh] = malloc (sizeof(int*) * n_part);
    _initToOlFace[imesh] = malloc (sizeof(int*) * n_part);
  }

  for (int imesh = 0; imesh < 2; imesh++) {
    PDM_ol_mesh_t mesht;

    if (imesh == 0) {
      mesht = PDM_OL_MESH_A;
    }
    else {
      mesht = PDM_OL_MESH_B;
    }

    PDM_g_num_t nGOlFace;
    PDM_g_num_t nGOlVtx;

    PDM_ol_mesh_dim_get (pdm_id,
                         mesht,
                         &nGOlFace,
                         &nGOlVtx);

    if (i_rank == 0) {
      printf("New mesh %d nGFace : "PDM_FMT_G_NUM"\n", imesh + 1, nGOlFace);
      printf("New mesh %d nGVtx : "PDM_FMT_G_NUM"\n", imesh + 1, nGOlVtx);
    }

    for (int ipart = 0; ipart < n_part; ipart++) {

      int           nOlFace;
      int           nOlLinkedFace;
      int           nOlVtx;
      int           sOlFaceIniVtx;
      int           sOlface_vtx;
      int           sInitToOlFace;

      PDM_ol_part_mesh_dim_get (pdm_id,
                                mesht,
                                ipart,
                                &nOlFace,
                                &nOlLinkedFace,
                                &nOlVtx,
                                &sOlFaceIniVtx,
                                &sOlface_vtx,
                                &sInitToOlFace);

      _nOlFace[imesh][ipart] = nOlFace;
      _nOlLinkedFace[imesh][ipart] = nOlLinkedFace;
      _nOlVtx[imesh][ipart] = nOlVtx;
      _sOlface_vtx[imesh][ipart] = sOlface_vtx;
      _sInitToOlFace[imesh][ipart] = sInitToOlFace;

      int            *olFaceIniVtxIdx;
      int            *olFaceIniVtx;
      int            *olface_vtx_idx;
      int            *olface_vtx;
      int            *olLinkedface_procIdx;
      int            *olLinkedFace;
      PDM_g_num_t     *olface_ln_to_gn;
      double         *olCoords;
      PDM_g_num_t     *olvtx_ln_to_gn;
      int            *initToOlFaceIdx;
      int            *initToOlFace;

      PDM_ol_mesh_entities_get (pdm_id,
                                mesht,
                                ipart,
                                &olFaceIniVtxIdx,
                                &olFaceIniVtx,
                                &olface_vtx_idx,
                                &olface_vtx,
                                &olLinkedface_procIdx,
                                &olLinkedFace,
                                &olface_ln_to_gn,
                                &olCoords,
                                &olvtx_ln_to_gn,
                                &initToOlFaceIdx,
                                &initToOlFace);

      _olFaceIniVtxIdx[imesh][ipart] = olFaceIniVtxIdx;
      _olFaceIniVtx[imesh][ipart] = olFaceIniVtx;
      _olface_vtx_idx[imesh][ipart] = olface_vtx_idx;
      _olface_vtx[imesh][ipart] = olface_vtx;
      _olLinkedface_procIdx[imesh][ipart] = olLinkedface_procIdx;
      _olLinkedFace[imesh][ipart] = olLinkedFace;
      _olface_ln_to_gn[imesh][ipart] = olface_ln_to_gn;
      _olCoords[imesh][ipart] = olCoords;
      _olvtx_ln_to_gn[imesh][ipart] = olvtx_ln_to_gn;
      _initToOlFaceIdx[imesh][ipart] = initToOlFaceIdx;
      _initToOlFace[imesh][ipart] = initToOlFace;

      if (imesh == 0) {
        sFieldOlA[ipart] = malloc(sizeof(double) * nOlFace);

        for (int i = 0; i < nFace[imesh][ipart]; i++) {
          for (int j = initToOlFaceIdx[i]; j < initToOlFaceIdx[i+1]; j++) {
            sFieldOlA[ipart][initToOlFace[j]-1] = sFieldA[ipart][i];
          }
        }
      }

      if (imesh == 1) {
        rFieldOlB[ipart] = malloc(sizeof(double) * nOlFace);

        surfB[ipart] = malloc(sizeof(double) * nFace[imesh][ipart]);
        surfOlB[ipart] = malloc(sizeof(double) * nOlFace);

        double   *ol_surface_vector = malloc(sizeof(double) * 3 * nOlFace);
        double   *ol_center = malloc(sizeof(double) * 3 * nOlFace);
        double   *ol_characteristicLength = malloc(sizeof(double) * nOlFace);
        int      *ol_isDegenerated = malloc(sizeof(int) * nOlFace);

        PDM_geom_elem_polygon_properties(nOlFace,
                                         _olface_vtx_idx[imesh][ipart],
                                         _olface_vtx[imesh][ipart],
                                         _olCoords[imesh][ipart],
                                         ol_surface_vector,
                                         ol_center,
                                         ol_characteristicLength,
                                         ol_isDegenerated);

        for (int i = 0; i < nOlFace; i++) {
          surfOlB[ipart][i] = PDM_MODULE(ol_surface_vector + 3*i);
          rFieldOlB[ipart][i] = -1.;//
        }

        for (int i = 0; i < nFace[imesh][ipart]; i++) {
          surfB[ipart][i] = 0.;
          for (int j = initToOlFaceIdx[i]; j < initToOlFaceIdx[i+1]; j++) {
            surfB[ipart][i] += surfOlB[ipart][initToOlFace[j]-1];
          }
        }

        free (ol_surface_vector);
        free (ol_center);
        free (ol_characteristicLength);
        free (ol_isDegenerated);
      }
    }
  }

  int *n_send = malloc (sizeof(int) * numProcs);
  int *n_recv = malloc (sizeof(int) * numProcs);
  for (int iproc = 0; iproc < numProcs; iproc++) {
    n_send[iproc] = 0;
    n_recv[iproc] = 0;
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int iproc = 0; iproc < numProcs; iproc++) {
      n_send[iproc] += _olLinkedface_procIdx[0][ipart][iproc+1] - _olLinkedface_procIdx[0][ipart][iproc];
      n_recv[iproc] += _olLinkedface_procIdx[1][ipart][iproc+1] - _olLinkedface_procIdx[1][ipart][iproc];
    }
  }

  int *i_send = malloc (sizeof(int) * (numProcs+1));
  int *i_recv = malloc (sizeof(int) * (numProcs+1));

  i_send[0] = 0;
  i_recv[0] = 0;
  for (int iproc = 0; iproc < numProcs; iproc++) {
    i_send[iproc+1] = i_send[iproc] + 4 * n_send[iproc];
    i_recv[iproc+1] = i_recv[iproc] + 4 * n_recv[iproc];
    n_send[iproc] = 0;
    n_recv[iproc] *= 4;
  }

  int *b_send = malloc (sizeof(int) * i_send[numProcs]);
  int *b_recv = malloc (sizeof(int) * i_recv[numProcs]);

  double *d_send = malloc (sizeof(double) * i_send[numProcs]);
  double *d_recv = malloc (sizeof(double) * i_recv[numProcs]);

  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int ilinked = 0; ilinked < _nOlLinkedFace[0][ipart]; ilinked++) {
      int iproc = _olLinkedFace[0][ipart][4*ilinked+1];
      d_send[(i_send[iproc]+n_send[iproc])/4] = sFieldOlA[ipart][_olLinkedFace[0][ipart][4*ilinked]-1];
      b_send[i_send[iproc]+n_send[iproc]++]   = ipart;
      b_send[i_send[iproc]+n_send[iproc]++]   = _olLinkedFace[0][ipart][4*ilinked];
      b_send[i_send[iproc]+n_send[iproc]++]   = _olLinkedFace[0][ipart][4*ilinked+2];
      b_send[i_send[iproc]+n_send[iproc]++]   = _olLinkedFace[0][ipart][4*ilinked+3];
    }
  }

  PDM_MPI_Alltoallv(b_send, n_send, i_send, PDM_MPI_INT,
                    b_recv, n_recv, i_recv, PDM_MPI_INT,
                    PDM_MPI_COMM_WORLD);

  for (int iproc = 0; iproc < numProcs + 1; iproc++) {
    i_send[iproc] /= 4;
    i_recv[iproc] /= 4;
  }

  for (int iproc = 0; iproc < numProcs; iproc++) {
    n_send[iproc] /= 4;
    n_recv[iproc] /= 4;
  }

  PDM_MPI_Alltoallv(d_send, n_send, i_send, PDM_MPI_DOUBLE,
                    d_recv, n_recv, i_recv, PDM_MPI_DOUBLE,
                    PDM_MPI_COMM_WORLD);

  free (n_send);
  free (b_send);
  free (i_send);
  free (n_recv);
  free (d_send);

  int **check_graph = malloc (sizeof(int*) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    check_graph[ipart] = malloc (sizeof(int) * 2 * _nOlFace[1][ipart]);
    for (int j = 0; j < 2 * _nOlFace[1][ipart]; j++) {
      check_graph[ipart][j] = -1;
    }
  }

  int idx = 0;
  int idx1 = 0;
  for (int ilinked = 0; ilinked < i_recv[numProcs]; ilinked++) {
    double val = d_recv[idx1++];
    int ipart1 = b_recv[idx++];
    int ielt1  = b_recv[idx++];
    int ipart2 = b_recv[idx++];
    int ielt2  = b_recv[idx++]-1;
    rFieldOlB[ipart2][ielt2] = val;
    check_graph[ipart2][2*ielt2] = ipart1;
    check_graph[ipart2][2*ielt2+1] = ielt1;
  }

  free (d_recv);

  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int i = 0; i < nFace[1][ipart]; i++) {
      rFieldB[ipart][i] = 0.;
      for (int j = _initToOlFaceIdx[1][ipart][i]; j < _initToOlFaceIdx[1][ipart][i+1]; j++) {
        rFieldB[ipart][i] += rFieldOlB[ipart][_initToOlFace[1][ipart][j]-1] *
                             surfOlB[ipart][_initToOlFace[1][ipart][j]-1];
      }
      rFieldB[ipart][i] /= surfB[ipart][i];
    }
  }


  for (int ipart2 = 0; ipart2 < n_part; ipart2++) {
    for (int ilinked = 0; ilinked < _nOlLinkedFace[1][ipart2]; ilinked++) {
      int ielt2  = _olLinkedFace[1][ipart2][4*ilinked]-1;
      int iproc1 = _olLinkedFace[1][ipart2][4*ilinked+2];
      int ipart1 = _olLinkedFace[1][ipart2][4*ilinked+2];
      int ielt1  = _olLinkedFace[1][ipart2][4*ilinked+3];
      if ((check_graph[ipart2][2*ielt2] != ipart1) ||
          (check_graph[ipart2][2*ielt2+1] != ielt1)) {
        printf("Erreur graph de comm : m1 %d %d %d / m2 %d %d\n", iproc1,
               ipart1, ielt1, check_graph[ipart2][2*ielt2], check_graph[ipart2][2*ielt2+1]);
        abort();
      }
    }
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    free (check_graph[ipart]);
  }
  free (check_graph);

  free (i_recv);
  free (b_recv);

  /*
   *  Export overlay mesh
   */

  if (post) {

    _export_ini_mesh (PDM_MPI_COMM_WORLD,
                      n_part,
                      nFace,
                      faceVtxIdx,
                      faceVtx,
                      faceLNToGN,
                      nVtx,
                      vtxCoord,
                      vtxLNToGN,
                      sFieldA,
                      rFieldB);

    _export_ol_mesh (PDM_MPI_COMM_WORLD,
                     pdm_id,
                     nFace,
                     sFieldOlA,
                     rFieldOlB,
                     n_part);

    _export_ol_vtk (PDM_MPI_COMM_WORLD,
                    pdm_id,
                    nFace,
                    sFieldOlA,
                    rFieldOlB,
                    n_part);
  }

  /*
   *  Free meshes
   */

  for (int ipart = 0; ipart < n_part; ipart++) {
    free (sFieldOlA[ipart]);
    free (sFieldA[ipart]);
    free (centerA[ipart]);
    free (rFieldOlB[ipart]);
    free (rFieldB[ipart]);
    free (surfB[ipart]);
    free (surfOlB[ipart]);
  }

  free (sFieldA);
  free (sFieldOlA);
  free (centerA);
  free (rFieldB);
  free (rFieldOlB);
  free (surfB);
  free (surfOlB);

  for (int imesh = 0; imesh < 2; imesh++) {

    for (int ipart = 0; ipart < n_part; ipart++) {
      free (faceVtxIdx[imesh][ipart]);
      free (faceVtx[imesh][ipart]);
      free (faceLNToGN[imesh][ipart]);
      free (vtxCoord[imesh][ipart]);
      free (vtxLNToGN[imesh][ipart]);
      free (_olFaceIniVtxIdx[imesh][ipart]);
      free (_olFaceIniVtx[imesh][ipart]);
      free (_olface_vtx_idx[imesh][ipart]);
      free (_olface_vtx[imesh][ipart]);
      free (_olLinkedface_procIdx[imesh][ipart]);
      free (_olLinkedFace[imesh][ipart]);
      free (_olface_ln_to_gn[imesh][ipart]);
      free (_olCoords[imesh][ipart]);
      free (_olvtx_ln_to_gn[imesh][ipart]);
      free (_initToOlFaceIdx[imesh][ipart]);
      free (_initToOlFace[imesh][ipart]);
    }

    free (faceVtxIdx[imesh]);
    free (faceVtx[imesh]);
    free (faceLNToGN[imesh]);
    free (vtxCoord[imesh]);
    free (vtxLNToGN[imesh]);
    free (nFace[imesh]);
    free (nVtx[imesh]);
    free (_olFaceIniVtxIdx[imesh]);
    free (_olFaceIniVtx[imesh]);
    free (_olface_vtx_idx[imesh]);
    free (_olface_vtx[imesh]);
    free (_olLinkedface_procIdx[imesh]);
    free (_olLinkedFace[imesh]);
    free (_olface_ln_to_gn[imesh]);
    free (_olCoords[imesh]);
    free (_olvtx_ln_to_gn[imesh]);
    free (_initToOlFaceIdx[imesh]);
    free (_initToOlFace[imesh]);

    free (_nOlFace[imesh]);
    free (_nOlLinkedFace[imesh]);
    free (_nOlVtx[imesh]);
    free (_sOlface_vtx[imesh]);
    free (_sInitToOlFace[imesh]);
  }
#endif

  /*
   *  Free Pdm
   */

  PDM_ol_del (pdm_id);

  PDM_MPI_Finalize ();

  return 0;

}
