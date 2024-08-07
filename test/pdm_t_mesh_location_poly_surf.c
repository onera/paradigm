#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
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
     "  -t      <level>  Bounding boxes tolerance.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -p      <level>  Number of points to locate.\n\n"
     "  -octree          Use octree-based method.\n\n"
     "  -dbbree          Use dbbtree-based method.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}



/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   tolerance  Bounding boxes tolerance
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   post       Ensight outputs status
 * \param [inout]   method     Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                         argc,
           char                      **argv,
           PDM_g_num_t                *n_vtx_seg,
           double                     *length,
           double                     *depth,
           int                        *rotation,
           double                     *tolerance,
           double                     *marge,
           int                        *n_part,
           PDM_g_num_t                *n_pts,
           int                        *n_proc_data,
           int                        *post,
           int                        *have_random,
           int                        *init_random,
           int                        *part_method,
           PDM_mesh_location_method_t *loc_method,
           int                        *use_mesh_vtx)
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
    else if (strcmp(argv[i], "-d") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *depth = atof(argv[i]);
    }
    else if (strcmp (argv[i], "-rot") == 0) {
      *rotation = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-m") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *marge = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_pts = atol(argv[i]);
        *n_pts = (PDM_g_num_t) _n_pts;
      }
    }

    else if (strcmp (argv[i], "-no_random") == 0) {
      *have_random = 0;
    }

    else if (strcmp (argv[i], "-use_vtx") == 0) {
      *use_mesh_vtx = 1;
    }

    else if (strcmp (argv[i], "-random_init") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *init_random = atoi(argv[i]);
      }
    }

    else if (strcmp (argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *n_proc_data = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_PART_SPLIT_HILBERT;
    }
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = PDM_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DBBTREE;
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}




static void
_point_cloud_from_mesh_vtx
(
 PDM_MPI_Comm   pdm_mpi_comm,
 double         xmin,
 double         ymin,
 PDM_g_num_t    n_vtxSeg,
 double         length,
 int            have_random,
 int            init_random,
 double       **coord,
 int           *n_pts_l
 )
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_mpi_comm, &n_rank);

  double       xmax = xmin + length;
  double       ymax = ymin + length;
  PDM_g_num_t  nx = n_vtxSeg;
  PDM_g_num_t  ny = n_vtxSeg;

  int          dn_face;
  int          dn_vtx;
  int          dn_edge;
  int         *dface_vtx_idx;
  PDM_g_num_t *dface_vtx;
  double      *dvtx_coord;
  PDM_g_num_t *dface_edge;
  PDM_g_num_t *dedge_vtx;
  PDM_g_num_t *dedge_face;
  int          n_edge_group;
  int         *dedge_group_idx;
  PDM_g_num_t *dedge_group;

  /*
   *  Create mesh
   */
  PDM_g_num_t n_g_face;
  PDM_g_num_t n_g_edge;
  PDM_g_num_t n_g_vtx;

  PDM_poly_surf_gen (pdm_mpi_comm,
                     xmin,
                     xmax,
                     ymin,
                     ymax,
                     have_random,
                     init_random,
                     nx,
                     ny,
                     &n_g_face,
                     &n_g_vtx,
                     &n_g_edge,
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

  *n_pts_l = dn_vtx;
  *coord   = dvtx_coord;

  PDM_free(dface_vtx_idx);
  PDM_free(dface_vtx);
  PDM_free(dedge_vtx);
  PDM_free(dedge_face);
  PDM_free(dedge_group_idx);
  PDM_free(dedge_group);
}


static void _add_depth (const int     n_pts,
                        const double  length,
                        const double  depth,
                        double       *coord)
{
  double inv_length = 1.;
  if (PDM_ABS (length) > 1e-15) inv_length /= length;

  for (int i = 0; i < n_pts; i++) {
    double x = 2.*coord[3*i]   * inv_length;
    double y = 2.*coord[3*i+1] * inv_length;
    coord[3*i+2] = 0.5*depth*(1. - (x*x + y*y));
  }
}


static void _rotate (const int  n_pts,
                     double    *coord)
{
  double R[3][3] = {{0.9362934, -0.2896295, 0.1986693},
                    {0.3129918,  0.9447025, -0.0978434},
                    {-0.1593451,  0.1537920,  0.9751703}};

  for (int i = 0; i < n_pts; i++) {
    double x = coord[3*i];
    double y = coord[3*i+1];
    double z = coord[3*i+2];

    for (int j = 0; j < 3; j++) {
      coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
    }
  }
}


static int _set_rank_has_mesh
(
 const PDM_MPI_Comm  comm,
 const int           n_proc_data,
 PDM_MPI_Comm       *mesh_comm
 )
{
  int current_rank_has_mesh = 1;

  int rank;
  int n_rank;

  PDM_MPI_Comm_rank (comm, &rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  if (n_proc_data > 0 && n_proc_data < n_rank) {

    int rank_in_node = PDM_io_mpi_node_rank (comm);

    int n_node = 0;
    int i_node = -1;
    int master_rank = (rank_in_node == 0);

    int *rank_in_nodes;
    PDM_malloc(rank_in_nodes, n_rank, int);

    PDM_MPI_Allreduce (&master_rank, &n_node, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    PDM_MPI_Allgather (&rank_in_node, 1, PDM_MPI_INT, rank_in_nodes, 1, PDM_MPI_INT, comm);

    current_rank_has_mesh = 0;

    for (int i = 0; i < rank; i++) {
      if (rank_in_nodes[i] == 0) {
        i_node += 1;
      }
    }

    if (n_proc_data <= n_node) {
      if (i_node < n_proc_data && master_rank) {
        current_rank_has_mesh = 1;
      }
    }

    else {

      if (rank_in_node < (n_proc_data / n_node)) {
        current_rank_has_mesh = 1;
      }
      if ((rank_in_node == (n_proc_data / n_node)) && (i_node < (n_proc_data % n_node))) {
        current_rank_has_mesh = 1;
      }

    }

    PDM_MPI_Comm_split (comm,
                        current_rank_has_mesh,
                        rank,
                        mesh_comm);
    PDM_free(rank_in_nodes);
  }

  return current_rank_has_mesh;
}


static void
_get_connectivity
(
 PDM_part_t    *ppart,
 int            n_part,
 int          **n_face,
 int         ***face_edge_idx,
 int         ***face_edge,
 int         ***face_vtx_idx,
 int         ***face_vtx,
 PDM_g_num_t ***face_ln_to_gn,
 int          **n_edge,
 int         ***edge_vtx_idx,
 int         ***edge_vtx,
 int          **n_vtx,
 double      ***vtx_coord,
 PDM_g_num_t ***vtx_ln_to_gn
 )
{
  PDM_malloc(*n_face       , n_part, int          );
  PDM_malloc(*face_edge_idx, n_part, int         *);
  PDM_malloc(*face_edge    , n_part, int         *);
  PDM_malloc(*face_vtx_idx , n_part, int         *);
  PDM_malloc(*face_vtx     , n_part, int         *);
  PDM_malloc(*face_ln_to_gn, n_part, PDM_g_num_t *);
  PDM_malloc(*n_edge       , n_part, int          );
  PDM_malloc(*edge_vtx_idx , n_part, int         *);
  PDM_malloc(*edge_vtx     , n_part, int         *);
  PDM_malloc(*n_vtx        , n_part, int          );
  PDM_malloc(*vtx_coord    , n_part, double      *);
  PDM_malloc(*vtx_ln_to_gn , n_part, PDM_g_num_t *);

   // int id_ppart = ppartId;

  for (int ipart = 0; ipart < n_part; ipart++) {

    int _n_face;
    int _n_edge;
    int _n_edgePartBound;
    int _n_vtx;
    int _nProc;
    int _nTPart;
    int _s_face_edge;
    int _s_edge_vtx;
    int _sEdgeGroup;
    int _n_edge_group2;

    PDM_part_part_dim_get (ppart,
                           ipart,
                           &_n_face,
                           &_n_edge,
                           &_n_edgePartBound,
                           &_n_vtx,
                           &_nProc,
                           &_nTPart,
                           &_s_face_edge,
                           &_s_edge_vtx,
                           &_sEdgeGroup,
                           &_n_edge_group2);

    int         *_faceTag;
    int         *_face_edge_idx;
    int         *_face_edge;
    PDM_g_num_t *_face_ln_to_gn;
    int         *_edge_tag;
    int         *_edge_face;
    int         *_edge_vtx_idx;
    int         *_edge_vtx;
    PDM_g_num_t *_edgeLNToGN;
    int         *_edgePartBoundProcIdx;
    int         *_edgePartBoundPartIdx;
    int         *_edgePartBound;
    int         *_vtxTag;
    double      *_vtx;
    PDM_g_num_t *_vtx_ln_to_gn;
    int         *_edgeGroupIdx;
    int         *_edgeGroup;
    PDM_g_num_t *_edgeGroupLNToGN;

    PDM_part_part_val_get (ppart,
                           ipart,
                           &_faceTag,
                           &_face_edge_idx,
                           &_face_edge,
                           &_face_ln_to_gn,
                           &_edge_tag,
                           &_edge_face,
                           &_edge_vtx_idx,
                           &_edge_vtx,
                           &_edgeLNToGN,
                           &_edgePartBoundProcIdx,
                           &_edgePartBoundPartIdx,
                           &_edgePartBound,
                           &_vtxTag,
                           &_vtx,
                           &_vtx_ln_to_gn,
                           &_edgeGroupIdx,
                           &_edgeGroup,
                           &_edgeGroupLNToGN);

    /* Faces */
    (*n_face)[ipart] = _n_face;
    PDM_malloc((*face_edge_idx)[ipart],  _n_face + 1, int        );
    PDM_malloc((*face_edge    )[ipart], _s_face_edge, int        );
    PDM_malloc((*face_vtx_idx )[ipart],  _n_face + 1, int        );
    PDM_malloc((*face_vtx     )[ipart], _s_face_edge, int        );
    PDM_malloc((*face_ln_to_gn)[ipart], _n_face     , PDM_g_num_t);

    memcpy ((*face_edge_idx)[ipart], _face_edge_idx, (_n_face + 1) * sizeof(int        ));
    memcpy ((*face_edge    )[ipart], _face_edge    ,  _s_face_edge * sizeof(int        ));
    memcpy ((*face_vtx_idx )[ipart], _face_edge_idx, (_n_face + 1) * sizeof(int        ));
    memcpy ((*face_ln_to_gn)[ipart], _face_ln_to_gn,  _n_face      * sizeof(PDM_g_num_t));

    /* Edges */
    (*n_edge)[ipart] = _n_edge;
    PDM_malloc((*edge_vtx_idx)[ipart], _n_edge + 1, int);
    PDM_malloc((*edge_vtx    )[ipart], _s_edge_vtx, int);

    memcpy ((*edge_vtx_idx)[ipart], _edge_vtx_idx, (_n_edge + 1) * sizeof(int));
    memcpy ((*edge_vtx    )[ipart], _edge_vtx    ,  _s_edge_vtx  * sizeof(int));

    /* Vertices */
    (*n_vtx)[ipart] = _n_vtx;
    PDM_malloc((*vtx_coord   )[ipart], 3 * _n_vtx, double     );
    PDM_malloc((*vtx_ln_to_gn)[ipart],     _n_vtx, PDM_g_num_t);

    memcpy ((*vtx_coord   )[ipart], _vtx         , 3 *_n_vtx * sizeof(double     ));
    memcpy ((*vtx_ln_to_gn)[ipart], _vtx_ln_to_gn,    _n_vtx * sizeof(PDM_g_num_t));

    /* Compute face-vtx connectivity */
    int *_face_vtx = (*face_vtx)[ipart];

    int *vtx_edge_idx = NULL;
    PDM_malloc(vtx_edge_idx, _n_vtx + 1, int);

    for (int i = 0; i < _n_vtx + 1; i++) {
      vtx_edge_idx[i] = 0;
    }

    for (int i = 0; i < _n_edge; i++) {
      int ivtx1 = _edge_vtx[2*i    ];
      int ivtx2 = _edge_vtx[2*i + 1];

      vtx_edge_idx[ivtx1] += 1;
      vtx_edge_idx[ivtx2] += 1;
    }

    for (int i = 1; i < _n_vtx + 1; i++) {
      vtx_edge_idx[i] = vtx_edge_idx[i] + vtx_edge_idx[i-1];
    }

    int *vtx_edge   = NULL;
    int *vtx_edge_n = NULL;
    PDM_malloc(vtx_edge  , vtx_edge_idx[_n_vtx], int);
    PDM_malloc(vtx_edge_n, _n_vtx              , int);
    for (int i = 0; i < _n_vtx; i++) {
      vtx_edge_n[i] = 0;
    }

    for (int i = 0; i < _n_edge; i++) {
      int ivtx1 = _edge_vtx[2*i    ] - 1;
      int ivtx2 = _edge_vtx[2*i + 1] - 1;
      int iedge = i + 1;

      vtx_edge[vtx_edge_idx[ivtx1] + vtx_edge_n[ivtx1]] = iedge;
      vtx_edge[vtx_edge_idx[ivtx2] + vtx_edge_n[ivtx2]] = iedge;
      vtx_edge_n[ivtx1] += 1;
      vtx_edge_n[ivtx2] += 1;
    }
    PDM_free(vtx_edge_n);

    for (int i = 0; i < _n_face; i++) {
      int idx = _face_edge_idx[i];
      int __n_edge = _face_edge_idx[i+1] - idx;
      int *_edges = _face_edge + idx;
      int *_vertices = _face_vtx + idx;

      int edge_cur = _edges[0];
      int vtx_deb =  _edge_vtx[2*(edge_cur - 1)];
      _vertices[0] = vtx_deb;
      int vtx_cur =  _edge_vtx[2*(edge_cur - 1) + 1];
      int idxVtx = 0;

      while (vtx_deb != vtx_cur) {
        _vertices[++idxVtx] = vtx_cur;
        int find_vtx = 0;

        for (int j = vtx_edge_idx[vtx_cur - 1]; j <  vtx_edge_idx[vtx_cur]; j++) {
          for (int k = 0; k < __n_edge; k++) {
            if ((_edges[k] == vtx_edge[j]) && (_edges[k] != edge_cur)) {
              edge_cur = _edges[k];
              if (_edge_vtx[2*(_edges[k]-1)] == vtx_cur) {
                vtx_cur = _edge_vtx[2*(_edges[k]-1) + 1];
              }
              else {
                vtx_cur = _edge_vtx[2*(_edges[k]-1)];
              }
              find_vtx = 1;
              break;
            }
          }
          if (find_vtx)
            break;
        }
        if (!find_vtx) {
          PDM_error(__FILE__, __LINE__, 0,"Error to compute vtx_edge !!!!\n");
          abort();
        }
      }
    }

    PDM_free(vtx_edge);
    PDM_free(vtx_edge_idx);

  }
}


static void
_create_split_mesh
(
 int                 activeRank,
 PDM_MPI_Comm        pdm_mpi_comm,
 double              xmin,
 double              ymin,
 PDM_g_num_t         n_vtxSeg,
 double              length,
 double              depth,
 int                 rotation,
 int                 n_part,
 PDM_part_split_t    method,
 int                 have_random,
 int                 init_random,
 PDM_g_num_t        *n_g_face,
 PDM_g_num_t        *n_g_vtx,
 int               **n_face,
 int              ***face_edge_idx,
 int              ***face_edge,
 int              ***face_vtx_idx,
 int              ***face_vtx,
 PDM_g_num_t      ***face_ln_to_gn,
 int               **n_edge,
 int              ***edge_vtx_idx,
 int              ***edge_vtx,
 int               **n_vtx,
 double           ***vtx_coord,
 PDM_g_num_t      ***vtx_ln_to_gn
 )
{
  int i_rank;
  int n_rank;

  if (activeRank) {

    PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
    PDM_MPI_Comm_size (pdm_mpi_comm, &n_rank);

    double       xmax = xmin + length;
    double       ymax = ymin + length;
    PDM_g_num_t  nx = n_vtxSeg;
    PDM_g_num_t  ny = n_vtxSeg;

    int          dn_face;
    int          dn_vtx;
    int          dn_edge;
    int         *dface_vtx_idx;
    PDM_g_num_t *dface_vtx;
    double      *dvtx_coord;
    PDM_g_num_t *dface_edge;
    PDM_g_num_t *dedge_vtx;
    PDM_g_num_t *dedge_face;
    int          n_edge_group;
    int         *dedge_group_idx;
    PDM_g_num_t *dedge_group;

    /*
     *  Create mesh
     */
    PDM_g_num_t n_g_edge;

    PDM_poly_surf_gen (pdm_mpi_comm,
                       xmin,
                       xmax,
                       ymin,
                       ymax,
                       have_random,
                       init_random,
                       nx,
                       ny,
                       n_g_face,
                       n_g_vtx,
                       &n_g_edge,
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

    _add_depth (dn_vtx,
                length,
                depth,
                dvtx_coord);

    if (rotation) {
      _rotate (dn_vtx,
               dvtx_coord);
    }

    /*
     *  Create mesh partitions
     */
    int have_dcell_part = 0;

    int *dcell_part    = NULL;
    int *dedge_vtx_idx = NULL;
    PDM_malloc(dcell_part   , dn_face    , int);
    PDM_malloc(dedge_vtx_idx, dn_edge + 1, int);

    dedge_vtx_idx[0] = 0;
    for (int i = 0; i < dn_edge; i++) {
      dedge_vtx_idx[i+1] = 2 + dedge_vtx_idx[i];
    }

    /*
     *  Split mesh
     */
    // int ppartId;

    int  n_property_cell       = 0;
    int *renum_properties_cell = NULL;
    int  n_property_face       = 0;
    int *renum_properties_face = NULL;

    PDM_part_t *ppart = PDM_part_create (pdm_mpi_comm,
                                         method,
                                         "PDM_PART_RENUM_CELL_NONE",
                                         "PDM_PART_RENUM_FACE_NONE",
                                         n_property_cell,
                                         renum_properties_cell,
                                         n_property_face,
                                         renum_properties_face,
                                         n_part,
                                         dn_face,
                                         dn_edge,
                                         dn_vtx,
                                         n_edge_group,
                                         NULL,
                                         NULL,
                                         NULL,
                                         NULL,
                                         have_dcell_part,
                                         dcell_part,
                                         dedge_face,
                                         dedge_vtx_idx,
                                         dedge_vtx,
                                         NULL,
                                         dvtx_coord,
                                         NULL,
                                         dedge_group_idx,
                                         dedge_group);

    PDM_free(dcell_part);

    PDM_free(dvtx_coord);
    PDM_free(dface_vtx_idx);
    PDM_free(dface_vtx);
    PDM_free(dface_edge);
    PDM_free(dedge_vtx_idx);
    PDM_free(dedge_vtx);
    PDM_free(dedge_face);
    PDM_free(dedge_group_idx);
    PDM_free(dedge_group);

    _get_connectivity (ppart,
                       n_part,
                       n_face,
                       face_edge_idx,
                       face_edge,
                       face_vtx_idx,
                       face_vtx,
                       face_ln_to_gn,
                       n_edge,
                       edge_vtx_idx,
                       edge_vtx,
                       n_vtx,
                       vtx_coord,
                       vtx_ln_to_gn);

    PDM_part_free (ppart);
  }

  else {
    PDM_malloc(*n_face       , n_part, int          );
    PDM_malloc(*face_edge_idx, n_part, int         *);
    PDM_malloc(*face_edge    , n_part, int         *);
    PDM_malloc(*face_vtx_idx , n_part, int         *);
    PDM_malloc(*face_vtx     , n_part, int         *);
    PDM_malloc(*face_ln_to_gn, n_part, PDM_g_num_t *);
    PDM_malloc(*n_edge       , n_part, int          );
    PDM_malloc(*edge_vtx_idx , n_part, int         *);
    PDM_malloc(*edge_vtx     , n_part, int         *);
    PDM_malloc(*n_vtx        , n_part, int          );
    PDM_malloc(*vtx_coord    , n_part, double      *);
    PDM_malloc(*vtx_ln_to_gn , n_part, PDM_g_num_t *);

    for (int ipart = 0; ipart < n_part; ipart++) {
      int _n_face = 0;
      int _s_face_edge = 0;
      (*n_face)[ipart] = _n_face;
      PDM_malloc((*face_edge_idx)[ipart], _n_face + 1 , int        );
      PDM_malloc((*face_vtx_idx )[ipart], _n_face + 1 , int        );
      PDM_malloc((*face_edge    )[ipart], _s_face_edge, int        );
      PDM_malloc((*face_vtx     )[ipart], _s_face_edge, int        );
      PDM_malloc((*face_ln_to_gn)[ipart], _n_face     , PDM_g_num_t);
      (*face_edge_idx)[ipart][0] = 0;
      (*face_vtx_idx )[ipart][0] = 0;

      int _n_edge = 0;
      int _s_edge_vtx = 0;
      (*n_edge)[ipart] = _n_edge;
      PDM_malloc((*edge_vtx_idx)[ipart], _n_edge + 1, int);
      PDM_malloc((*edge_vtx    )[ipart], _s_edge_vtx, int);
      (*edge_vtx_idx)[ipart][0] = 0;

      int _n_vtx = 0;
      (*n_vtx)[ipart] = _n_vtx;
      PDM_malloc((*vtx_coord   )[ipart], 3 * _n_vtx, double     );
      PDM_malloc((*vtx_ln_to_gn)[ipart],     _n_vtx, PDM_g_num_t);
    }
  }

  PDM_MPI_Bcast (n_g_face, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
  PDM_MPI_Bcast (n_g_vtx , 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
}


static inline double
_eval_field
(
 double *xyz
 )
{
  return 1 + 2*xyz[0] + 3*xyz[1] + 4*xyz[2];
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

  PDM_g_num_t n_vtx_seg   = 10;
  double      length      = 1.;
  double      depth       = 0.;
  int         rotation    = 0;
  double      tolerance   = 1e-6;
  double      marge       = 0.;
  int         n_part      = 1;
  int         post        = 0;
  int         have_random = 1;
  int         init_random = 0;
  int         n_proc_data = -1;
  PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;

  PDM_g_num_t n_pts = 10;
  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;

  int use_vtx = 0;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_seg,
              &length,
              &depth,
              &rotation,
              &tolerance,
              &marge,
              &n_part,
              &n_pts,
              &n_proc_data,
              &post,
              &have_random,
              &init_random,
              (int *) &part_method,
              &loc_method,
              &use_vtx);


  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_rank      : %d\n", n_rank);
    PDM_printf ("  - n_vtx_seg   : "PDM_FMT_G_NUM"\n", n_vtx_seg);
    PDM_printf ("  - n_pts       : "PDM_FMT_G_NUM"\n", n_pts);
    PDM_printf ("  - length      : %f\n", length);
    PDM_printf ("  - depth       : %f\n", depth);
    PDM_printf ("  - tolerance   : %f\n", tolerance);
    PDM_printf ("  - part_method : %d\n", (int) part_method);
    PDM_printf ("  - loc_method  : %d\n", (int) loc_method);
    PDM_printf ("  - n_proc_data : %d\n", n_proc_data);
  }

  /*
   *  Create partitionned surface mesh
   */

  if (i_rank == 0) {
    printf("-- Build surface mesh\n");
    fflush(stdout);
  }

  PDM_MPI_Comm mesh_comm = PDM_MPI_COMM_WORLD;
  int current_rank_has_mesh = _set_rank_has_mesh (PDM_MPI_COMM_WORLD,
                                                  n_proc_data,
                                                  &mesh_comm);

  const double xmin = -0.5*length;
  const double ymin = -0.5*length;

  PDM_g_num_t   n_g_face;
  PDM_g_num_t   n_g_vtx;
  int          *n_face        = NULL;
  PDM_g_num_t **face_ln_to_gn = NULL;
  int         **face_edge_idx = NULL;
  int         **face_edge     = NULL;
  int         **face_vtx_idx  = NULL;
  int         **face_vtx      = NULL;
  int          *n_edge        = NULL;
  int         **edge_vtx_idx  = NULL;
  int         **edge_vtx      = NULL;
  int          *n_vtx         = NULL;
  double      **vtx_coord     = NULL;
  PDM_g_num_t **vtx_ln_to_gn  = NULL;

  _create_split_mesh (current_rank_has_mesh,
                      mesh_comm,
                      xmin,
                      ymin,
                      n_vtx_seg,
                      length,
                      depth,
                      rotation,
                      n_part,
                      part_method,
                      have_random,
                      init_random,
                      &n_g_face,
                      &n_g_vtx,
                      &n_face,
                      &face_edge_idx,
                      &face_edge,
                      &face_vtx_idx,
                      &face_vtx,
                      &face_ln_to_gn,
                      &n_edge,
                      &edge_vtx_idx,
                      &edge_vtx,
                      &n_vtx,
                      &vtx_coord,
                      &vtx_ln_to_gn);


  /************************
   *
   * Point cloud definition
   *
   ************************/
  if (i_rank == 0) {
    printf("-- Point cloud\n");
    fflush(stdout);
  }

  int n_pts_l;
  double *pts_coords    = NULL;
  PDM_g_num_t *pts_gnum = NULL;

  marge *= length;
  double _min = -0.5*length - marge;
  double _max = -_min;

  if (use_vtx) {
    _point_cloud_from_mesh_vtx (PDM_MPI_COMM_WORLD,
                                _min,
                                _min,
                                n_pts,
                                length + 2.*marge,
                                have_random,
                                init_random + 1,
                                &pts_coords,
                                &n_pts_l);
  }
  else {
    PDM_point_cloud_gen_random (PDM_MPI_COMM_WORLD,
                                0, // seed
                                0, // geometric_g_num
                                n_pts,
                                _min,
                                _min,
                                0.,
                                _max,
                                _max,
                                0.,
                                &n_pts_l,
                                &pts_coords,
                                &pts_gnum);
  }

  _add_depth (n_pts_l,
              length,
              depth,
              pts_coords);


  /* Point cloud global numbering */
  if (use_vtx) {
    if (i_rank == 0) {
      printf("-- Point cloud g_num\n");
      fflush(stdout);
    }
  #if 1
    PDM_gen_gnum_t* gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_USER);

    double *char_length;
    PDM_malloc(char_length, n_pts_l, double);

    for (int i = 0; i < n_pts_l; i++) {
      char_length[i] = length * 1.e-9;
    }

    PDM_gnum_set_from_coords (gen_gnum, 0, n_pts_l, pts_coords, char_length);

    if (i_rank == 0) {
      printf(">> PDM_gnum_compute\n");
      fflush(stdout);
    }
    PDM_gnum_compute (gen_gnum);
    if (i_rank == 0) {
      printf("<< PDM_gnum_compute\n");
      fflush(stdout);
    }

    pts_gnum = PDM_gnum_get(gen_gnum, 0);

    PDM_gnum_free (gen_gnum);
    PDM_free(char_length);
  #else
  PDM_g_num_t *distrib = PDM_compute_entity_distribution (PDM_MPI_COMM_WORLD,
                                                          n_pts_l);
  PDM_g_num_t *pts_gnum;
  PDM_malloc(pts_gnum, n_pts_l, PDM_g_num_t);
  for (int i = 0; i < n_pts_l; i++) {
    pts_gnum[i] = distrib[i_rank] + i + 1;
  }

  PDM_free(distrib);
  #endif
  }


  if (post) {
    char filename[999];
    sprintf(filename, "point_cloud_%3.3d.vtk", i_rank);

    PDM_vtk_write_point_cloud (filename,
                               n_pts_l,
                               pts_coords,
                               pts_gnum,
                               NULL);
  }


  /************************
   *
   * Mesh location struct initializaiton
   *
   ************************/
  if (i_rank == 0) {
    printf("-- Create mesh loc\n");
    fflush(stdout);
  }

  PDM_mesh_location_t* mesh_loc = PDM_mesh_location_create (1,//const int n_point_cloud,
                                                            PDM_MPI_COMM_WORLD,
                                                            PDM_OWNERSHIP_KEEP);

  /* Set point cloud(s) */
  PDM_mesh_location_n_part_cloud_set (mesh_loc,
                                      0,//i_point_cloud,
                                      1);//n_part

  PDM_mesh_location_cloud_set (mesh_loc,
                               0,//i_point_cloud,
                               0,//i_part,
                               n_pts_l,
                               pts_coords,
                               pts_gnum);

  PDM_mesh_location_mesh_n_part_set (mesh_loc,
                                          n_part);

  /* Set mesh */
  if (i_rank == 0) {
    printf("-- Set mesh\n");
    fflush(stdout);
  }
  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_mesh_location_part_set_2d (mesh_loc,
                                   ipart,
                                   n_face       [ipart],
                                   face_edge_idx[ipart],
                                   face_edge    [ipart],
                                   face_ln_to_gn[ipart],
                                   n_edge       [ipart],
                                   edge_vtx     [ipart],
                                   n_vtx        [ipart],
                                   vtx_coord    [ipart],
                                   vtx_ln_to_gn [ipart]);
  }




  /* Set location parameters */
  PDM_mesh_location_tolerance_set (mesh_loc,
                                   tolerance);

  PDM_mesh_location_method_set (mesh_loc,
                                loc_method);


  /*
   * Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  // PDM_mesh_location_compute (mesh_loc);
  PDM_mesh_location_compute(mesh_loc);

  PDM_mesh_location_dump_times (mesh_loc);





  /*
   * Check results
   */
  if (i_rank == 0) {
    printf("-- Check\n");
    fflush(stdout);
  }

  int n_located = PDM_mesh_location_n_located_get (mesh_loc,
                                                   0,//i_point_cloud,
                                                   0);//i_part,

  int *located = PDM_mesh_location_located_get (mesh_loc,
                                                0,//i_point_cloud,
                                                0);

  int n_unlocated = PDM_mesh_location_n_unlocated_get (mesh_loc,
                                                       0,//i_point_cloud,
                                                       0);

  int *unlocated = PDM_mesh_location_unlocated_get (mesh_loc,
                                                    0,//i_point_cloud,
                                                    0);

  PDM_g_num_t *p_location    = NULL;
  double      *p_dist2  = NULL;
  double      *p_proj_coord  = NULL;
  PDM_mesh_location_point_location_get (mesh_loc,
                                        0,//i_point_cloud,
                                        0,//i_part,
                                        &p_location,
                                        &p_dist2,
                                        &p_proj_coord);

  if (0) {
    printf("Unlocated %d :\n", n_unlocated);
    for (int k1 = 0; k1 < n_unlocated; k1++) {
      printf("%d\n", unlocated[k1]);
    }
    printf("\n");

    printf("Located %d :\n", n_located);
    for (int k1 = 0; k1 < n_located; k1++) {
      printf("%d\n", located[k1]);
    }
    printf("\n");

    printf("Located %d :\n", n_located);
    for (int k1 = 0; k1 < n_located; k1++) {
      int ipt = located[k1] - 1;
      printf(PDM_FMT_G_NUM" : "PDM_FMT_G_NUM" / %12.5e %12.5e %12.5e / %12.5e / %12.5e %12.5e %12.5e",
        pts_gnum[ipt],  p_location[k1],
        pts_coords[3*ipt], pts_coords[3*ipt+1], pts_coords[3*ipt+2],
        p_dist2[k1],
        p_proj_coord[3*k1], p_proj_coord[3*k1+1], p_proj_coord[3*k1+2]);
      printf("\n");
    }
  }


  // /* int p_n_points; */
  // PDM_g_num_t *p_location    = NULL;
  // int         *p_weights_idx = NULL;
  // double      *p_weights     = NULL;
  // double      *p_proj_coord  = NULL;

  // PDM_mesh_location_get (mesh_loc,
  //                        0,//i_point_cloud,
  //                        0,//i_part,
  //                        &p_location,
  //                        &p_weights_idx,
  //                        &p_weights,
  //                        &p_proj_coord);

  if (0) {
    printf("Unlocated %d :\n", n_unlocated);
    for (int k1 = 0; k1 < n_unlocated; k1++) {
      printf("%d\n", unlocated[k1]);
    }
    printf("\n");

    printf("Located %d :\n", n_located);
    for (int k1 = 0; k1 < n_located; k1++) {
      printf("%d\n", located[k1]);
    }
    printf("\n");

    printf("Located %d :\n", n_located);
    for (int k1 = 0; k1 < n_located; k1++) {
      int ipt = located[k1] - 1;
      printf(PDM_FMT_G_NUM" : "PDM_FMT_G_NUM" / %12.5e %12.5e %12.5e / %12.5e / %12.5e %12.5e %12.5e",
        pts_gnum[ipt],  p_location[k1],
        pts_coords[3*ipt], pts_coords[3*ipt+1], pts_coords[3*ipt+2],
        p_dist2[k1],
        p_proj_coord[3*k1], p_proj_coord[3*k1+1], p_proj_coord[3*k1+2]);
      printf("\n");
    }
  }

  //  if (0) {
  //   for (int ipt = 0; ipt < n_pts_l; ipt++) {
  //     printf("Point ("PDM_FMT_G_NUM") (%f %f %f), location = ("PDM_FMT_G_NUM"), proj = (%f %f %f), weights =",
  //            pts_gnum[ipt],
  //            pts_coords[3*ipt], pts_coords[3*ipt+1], pts_coords[3*ipt+2],
  //            p_location[ipt],
  //            p_proj_coord[3*ipt], p_proj_coord[3*ipt+1], p_proj_coord[3*ipt+2]);
  //     for (int i = p_weights_idx[ipt]; i < p_weights_idx[ipt+1]; i++) {
  //       printf(" %f", p_weights[i]);
  //     }
  //     printf("\n");
  //   }
  // }




  /*
   *  Check location (interpolation of an affine field)
   */
  double **src_field = NULL;
  PDM_malloc(src_field, n_part, double *);
  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_malloc(src_field[ipart], n_vtx[ipart], double);
    for (int i = 0; i < n_vtx[ipart]; i++) {
      src_field[ipart][i] = _eval_field(&vtx_coord[ipart][3*i]);
    }
  }


  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_USER);
  if (ptp == NULL) {
    int         **pelt_pts_idx  = NULL;
    PDM_g_num_t **pelt_pts_gnum = NULL;
    PDM_malloc(pelt_pts_idx , n_part, int         *);
    PDM_malloc(pelt_pts_gnum, n_part, PDM_g_num_t *);
    for (int ipart = 0; ipart < n_part; ipart++) {
      double *elt_pts_coord      = NULL;
      double *elt_pts_uvw        = NULL;
      int    *elt_pts_weight_idx = NULL;
      double *elt_pts_weight     = NULL;
      double *elt_pts_dist2      = NULL;
      double *elt_pts_proj_coord = NULL;
      PDM_mesh_location_points_in_elt_get(mesh_loc,
                                          ipart,
                                          0, // i_point_cloud,
                                          &pelt_pts_idx [ipart],
                                          &pelt_pts_gnum[ipart],
                                          &elt_pts_coord,
                                          &elt_pts_uvw,
                                          &elt_pts_weight_idx,
                                          &elt_pts_weight,
                                          &elt_pts_dist2,
                                          &elt_pts_proj_coord);
    }

    ptp = PDM_part_to_part_create((const PDM_g_num_t **) face_ln_to_gn,
                                  (const int          *) n_face,
                                                         n_part,
                                  (const PDM_g_num_t **) &pts_gnum,
                                  (const int          *) &n_pts_l,
                                                         1,
                                  (const int         **) pelt_pts_idx,
                                  (const PDM_g_num_t **) pelt_pts_gnum,
                                  PDM_MPI_COMM_WORLD);

    PDM_free(pelt_pts_idx );
    PDM_free(pelt_pts_gnum);
  }


  double **send_field = NULL;
  PDM_malloc(send_field, n_part, double *);
  for (int ipart = 0; ipart < n_part; ipart++) {
    int         *elt_pts_idx        = NULL;
    PDM_g_num_t *elt_pts_gnum       = NULL;
    double      *elt_pts_coord      = NULL;
    double      *elt_pts_uvw        = NULL;
    int         *elt_pts_weight_idx = NULL;
    double      *elt_pts_weight     = NULL;
    double      *elt_pts_dist2      = NULL;
    double      *elt_pts_proj_coord = NULL;
    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        ipart,
                                        0, // i_point_cloud,
                                        &elt_pts_idx,
                                        &elt_pts_gnum,
                                        &elt_pts_coord,
                                        &elt_pts_uvw,
                                        &elt_pts_weight_idx,
                                        &elt_pts_weight,
                                        &elt_pts_dist2,
                                        &elt_pts_proj_coord);

    PDM_malloc(send_field[ipart], elt_pts_idx[n_face[ipart]], double);
    for (int ielt = 0; ielt < n_face[ipart]; ielt++) {

      int *fv = face_vtx[ipart] + face_vtx_idx[ipart][ielt];

      for (int idx_pt = elt_pts_idx[ielt]; idx_pt < elt_pts_idx[ielt+1]; idx_pt++) {
        send_field[ipart][idx_pt] = 0.;
        int idx_vtx = 0;
        // double e[3] = {
        //   elt_pts_proj_coord[3*idx_pt  ],
        //   elt_pts_proj_coord[3*idx_pt+1],
        //   elt_pts_proj_coord[3*idx_pt+2]
        // };
        double e[3] = {
          elt_pts_coord[3*idx_pt  ],
          elt_pts_coord[3*idx_pt+1],
          elt_pts_coord[3*idx_pt+2]
        };

        for (int idx_w = elt_pts_weight_idx[idx_pt]; idx_w < elt_pts_weight_idx[idx_pt+1]; idx_w++) {
          int vtx_id = fv[idx_vtx++] - 1;
          send_field[ipart][idx_pt] += elt_pts_weight[idx_w] * src_field[ipart][vtx_id];
          for (int j = 0; j < 3; j++) {
            e[j] -= elt_pts_weight[idx_w] * vtx_coord[ipart][3*vtx_id+j];
          }
        }

        // log_trace("pt "PDM_FMT_G_NUM" (%f %f %f), in elt "PDM_FMT_G_NUM" : dist = %e\n",
        //           elt_pts_gnum[idx_pt],
        //           elt_pts_coord[3*idx_pt], elt_pts_coord[3*idx_pt+1], elt_pts_coord[3*idx_pt+2],
        //           face_ln_to_gn[ipart][ielt],
        //           PDM_MODULE(e));
      }

    }
  }


  double **recv_field = NULL;
  int request = -1;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) send_field,
                         NULL,
        (      void ***) &recv_field,
                         &request);

  PDM_part_to_part_iexch_wait(ptp, request);

  double lmax_err = 0.;
  for (int i = 0; i < n_located; i++) {
    int pt_id = located[i] - 1;

    double f = _eval_field(&pts_coords[3*pt_id]);

    double err = PDM_ABS(recv_field[0][i] - f);
    lmax_err = PDM_MAX(lmax_err, err);

    if (err > 1.e-12) {
      log_trace("point "PDM_FMT_G_NUM" (%f %f %f) located in elt "PDM_FMT_G_NUM" : error = %e (%20.16f / %20.16f)\n",
                pts_gnum[pt_id],
                pts_coords[3*pt_id], pts_coords[3*pt_id+1], pts_coords[3*pt_id+2],
                p_location[i], err, recv_field[0][i], f);
    }
  }

  PDM_free(recv_field[0]);
  PDM_free(recv_field);


  double gmax_err;
  PDM_MPI_Allreduce(&lmax_err, &gmax_err, 1, PDM_MPI_DOUBLE,
                    PDM_MPI_MAX, PDM_MPI_COMM_WORLD);


  if (i_rank == 0) {
    printf("global max interpolation error = %e\n", gmax_err);
  }


  /*
   * Finalize
   */
  PDM_mesh_location_free(mesh_loc);
  PDM_part_to_part_free (ptp);
                          

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_free(face_edge_idx[ipart]);
    PDM_free(face_edge    [ipart]);
    PDM_free(face_vtx_idx [ipart]);
    PDM_free(face_vtx     [ipart]);
    PDM_free(face_ln_to_gn[ipart]);
    PDM_free(edge_vtx_idx [ipart]);
    PDM_free(edge_vtx     [ipart]);
    PDM_free(vtx_coord    [ipart]);
    PDM_free(vtx_ln_to_gn [ipart]);

    PDM_free(src_field    [ipart]);
    PDM_free(send_field   [ipart]);
  }
  PDM_free(face_vtx_idx);
  PDM_free(face_vtx);
  PDM_free(n_face);
  PDM_free(face_edge_idx);
  PDM_free(face_edge);
  PDM_free(face_ln_to_gn);
  PDM_free(n_edge);
  PDM_free(edge_vtx_idx);
  PDM_free(edge_vtx);
  PDM_free(n_vtx);
  PDM_free(vtx_coord);
  PDM_free(vtx_ln_to_gn);

  PDM_free(src_field);
  PDM_free(send_field);
  /*PDM_part_free (ppart_id);


  PDM_free(dvtx_coord);
  PDM_free(dface_vtx_idx);
  PDM_free(dface_vtx);
  PDM_free(dface_edge);
  PDM_free(dedge_vtx);
  PDM_free(dedge_face);
  PDM_free(dedge_group_idx);
  PDM_free(dedge_group);
  PDM_free(dedge_vtx_idx);*/

  PDM_free(pts_coords);
  PDM_free(pts_gnum);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
