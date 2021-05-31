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
#include "pdm_part.h"

#include "pdm_writer.h"
#include "pdm_geom_elem.h"
#include "pdm_part_to_block.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_overlay.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_priv.h"

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
 PDM_g_num_t  *nVtxSegA,
 double        *lengthA,
 double        *xminA,
 double        *yminA,
 int           *n_partA,
 PDM_g_num_t  *nVtxSegB,
 double        *lengthB,
 double        *xminB,
 double        *yminB,
 int           *n_partB,
 int           *post,
 int           *method,
 int           *haveRandom,
 int           *randomTimeInit,
 int           *randomMeshAInit,
 int           *randomMeshBInit,
 int           *quadA,
 int           *quadB,
 int           *nProcData,
 int           *rotate
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
    else if (strcmp (argv[i], "-quadA") == 0) {
      *quadA = 1;
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
    else if (strcmp (argv[i], "-randomMeshBInit") == 0) {
      i++;
      if (i >= argc)
        _usage (EXIT_FAILURE);
      else {
        *randomMeshBInit = atoi (argv[i]);
      }
    }
    else if (strcmp (argv[i], "-quadB") == 0) {
      *quadB = 1;
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
    else if (strcmp (argv[i], "-rotate") == 0) {
      *rotate = 1;
    }
    else
      _usage (EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Compute faceVtx connectivity
 *
 * \param [in]      ppartId  ppart identifier
 * \param [in]      n_part    Number of partitions
 *
 * \return          faceVtx connectivity for each partition of each mesh
 */

static void
_compute_faceVtx
(
 int           ppartId,
 int            n_part,
 int          **nFace,
 int         ***faceVtxIdx,
 int         ***faceVtx,
 PDM_g_num_t ***faceLNToGN,
 int          **nVtx,
 double      ***vtxCoord,
 PDM_g_num_t ***vtxLNToGN
)
{

  *nFace = (int *) malloc(sizeof(int) * n_part);
  *faceVtxIdx = (int **) malloc(sizeof(int *) * n_part);
  *faceVtx = (int **) malloc(sizeof(int *) * n_part);
  *faceLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  *nVtx = (int *) malloc(sizeof(int) * n_part);
  *vtxCoord = (double **) malloc(sizeof(double *) * n_part);
  *vtxLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  int id_ppart = ppartId;

  for (int ipart = 0; ipart < n_part; ipart++) {

    int _nFace;
    int _nEdge;
    int _nEdgePartBound;
    int _nVtx;
    int _nProc;
    int _nTPart;
    int _sFaceEdge;
    int _sEdgeVtx;
    int _sEdgeGroup;
    int _nEdgeGroup2;

    PDM_part_part_dim_get (id_ppart,
                           ipart,
                           &_nFace,
                           &_nEdge,
                           &_nEdgePartBound,
                           &_nVtx,
                           &_nProc,
                           &_nTPart,
                           &_sFaceEdge,
                           &_sEdgeVtx,
                           &_sEdgeGroup,
                           &_nEdgeGroup2);

    int          *_faceTag;
    int          *_faceEdgeIdx;
    int          *_faceEdge;
    PDM_g_num_t *_faceLNToGN;
    int          *_edgeTag;
    int          *_edgeFace;
    int          *_edgeVtxIdx;
    int          *_edgeVtx;
    PDM_g_num_t *_edgeLNToGN;
    int          *_edgePartBoundProcIdx;
    int          *_edgePartBoundPartIdx;
    int          *_edgePartBound;
    int          *_vtxTag;
    double       *_vtx;
    PDM_g_num_t *_vtxLNToGN;
    int          *_edgeGroupIdx;
    int          *_edgeGroup;
    PDM_g_num_t  *_edgeGroupLNToGN;

    PDM_part_part_val_get (id_ppart,
                           ipart,
                           &_faceTag,
                           &_faceEdgeIdx,
                           &_faceEdge,
                           &_faceLNToGN,
                           &_edgeTag,
                           &_edgeFace,
                           &_edgeVtxIdx,
                           &_edgeVtx,
                           &_edgeLNToGN,
                           &_edgePartBoundProcIdx,
                           &_edgePartBoundPartIdx,
                           &_edgePartBound,
                           &_vtxTag,
                           &_vtx,
                           &_vtxLNToGN,
                           &_edgeGroupIdx,
                           &_edgeGroup,
                           &_edgeGroupLNToGN);

    (*nFace)[ipart] = _nFace;
    (*faceVtxIdx) [ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
    (*faceVtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
    (*faceLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);

    memcpy ((*faceVtxIdx)[ipart], _faceEdgeIdx, (_nFace + 1) * sizeof(int));
    memcpy ((*faceLNToGN)[ipart], _faceLNToGN, _nFace * sizeof(PDM_g_num_t));

    (*nVtx)[ipart] = _nVtx;
    (*vtxCoord)[ipart] = (double *) malloc(sizeof(double) * (3 * _nVtx));
    (*vtxLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);

    memcpy ((*vtxCoord)[ipart], _vtx, 3 *_nVtx * sizeof(double));
    memcpy ((*vtxLNToGN)[ipart], _vtxLNToGN, _nVtx * sizeof(PDM_g_num_t));

    int *_faceVtx = (*faceVtx)[ipart];

    int *vtxEdgeIdx = (int *) malloc(sizeof(int) * (_nVtx + 1));

    for (int i = 0; i < _nVtx + 1; i++) {
      vtxEdgeIdx[i] = 0;
    }

    for (int i = 0; i < _nEdge; i++) {
      int ivtx1 = _edgeVtx[2*i];
      int ivtx2 = _edgeVtx[2*i + 1];

      vtxEdgeIdx[ivtx1] += 1;
      vtxEdgeIdx[ivtx2] += 1;
    }

    for (int i = 1; i < _nVtx + 1; i++) {
      vtxEdgeIdx[i] = vtxEdgeIdx[i] + vtxEdgeIdx[i-1];
    }

    int *vtxEdge = (int *) malloc(sizeof(int) * vtxEdgeIdx[_nVtx]);
    int *vtxEdgeN = (int *) malloc(sizeof(int) * _nVtx);
    for (int i = 0; i < _nVtx; i++) {
      vtxEdgeN[i] = 0;
    }

    for (int i = 0; i < _nEdge; i++) {
      int ivtx1 = _edgeVtx[2*i] - 1;
      int ivtx2 = _edgeVtx[2*i + 1] - 1;
      int iedge = i + 1;

      vtxEdge[vtxEdgeIdx[ivtx1] + vtxEdgeN[ivtx1]] = iedge;
      vtxEdge[vtxEdgeIdx[ivtx2] + vtxEdgeN[ivtx2]] = iedge;
      vtxEdgeN[ivtx1] += 1;
      vtxEdgeN[ivtx2] += 1;
    }
    free(vtxEdgeN);

    for (int i = 0; i < _nFace; i++) {
      int idx = _faceEdgeIdx[i];
      int __nEdge = _faceEdgeIdx[i+1] - idx;
      int *_edges = _faceEdge + idx;
      int *_vertices = _faceVtx + idx;

      int edge_cur = _edges[0];
      int vtx_deb =  _edgeVtx[2*(edge_cur - 1)];
      _vertices[0] = vtx_deb;
      int vtx_cur =  _edgeVtx[2*(edge_cur - 1) + 1];
      int idxVtx = 0;

      while (vtx_deb != vtx_cur) {
        _vertices[++idxVtx] = vtx_cur;
        int find_vtx = 0;

        for (int j = vtxEdgeIdx[vtx_cur - 1]; j <  vtxEdgeIdx[vtx_cur]; j++) {
          for (int k = 0; k < __nEdge; k++) {
            if ((_edges[k] == vtxEdge[j]) && (_edges[k] != edge_cur)) {
              edge_cur = _edges[k];
              if (_edgeVtx[2*(_edges[k]-1)] == vtx_cur) {
                vtx_cur = _edgeVtx[2*(_edges[k]-1) + 1];
              }
              else {
                vtx_cur = _edgeVtx[2*(_edges[k]-1)];
              }
              find_vtx = 1;
              break;
            }
          }
          if (find_vtx)
            break;
        }
        if (!find_vtx) {
          PDM_error(__FILE__, __LINE__, 0,"Error to compute vtxedge !!!!\n");
          abort();
        }
      }
    }

    free (vtxEdge);
    free (vtxEdgeIdx);

  }

}



static void
_rotation_matrix
(
 double rot[3][3]
 )
{
  double q[4] = {1., -2., 3., -4};

  double l = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  q[0] /= l;
  q[1] /= l;
  q[2] /= l;
  q[3] /= l;

  rot[0][0] = 1. - 2.*(q[1]*q[1] + q[2]*q[2]);
  rot[0][1] =      2.*(q[0]*q[1] - q[2]*q[3]);
  rot[0][2] =      2.*(q[0]*q[2] + q[1]*q[3]);
  rot[1][0] =      2.*(q[0]*q[1] + q[2]*q[3]);
  rot[1][1] = 1. - 2.*(q[0]*q[0] + q[2]*q[2]);
  rot[1][2] =      2.*(q[1]*q[2] - q[0]*q[3]);
  rot[2][0] =      2.*(q[0]*q[2] - q[1]*q[3]);
  rot[2][1] =      2.*(q[1]*q[2] + q[0]*q[3]);
  rot[2][2] = 1. - 2.*(q[0]*q[0] + q[1]*q[1]);
}





static double random01(void)
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / PDM_ABS(rsigna - rsignb);
  double resultat = sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}

static void
_quad_mesh_gen
(
 PDM_MPI_Comm        pdm_comm,
 double              xmin,
 double              xmax,
 double              ymin,
 double              ymax,
 int                 have_random,
 int                 init_random,
 PDM_g_num_t         nx,
 PDM_g_num_t         ny,
 PDM_g_num_t        *ng_face,
 PDM_g_num_t        *ng_vtx,
 PDM_g_num_t        *ng_edge,
 int                *dn_vtx,
 double            **dvtx_coord,
 int                *dn_face,
 int               **dface_vtx_idx,
 PDM_g_num_t       **dface_vtx,
 PDM_g_num_t       **dface_edge,
 int                *dn_edge,
 PDM_g_num_t       **dedge_vtx,
 PDM_g_num_t       **dedge_face,
 int                *n_edge_group,
 int               **dedge_group_idx,
 PDM_g_num_t       **dedge_group
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  const double coef_rand = 0.3 * have_random;
  srand(init_random);

  const PDM_g_num_t nx1 = nx - 1;
  const PDM_g_num_t ny1 = ny - 1;

  *ng_face = nx1 * ny1;
  *ng_edge = nx1 * ny + nx * ny1;
  *ng_vtx  = nx * ny;
  *n_edge_group = 4;

  PDM_g_num_t ng_edge_lim = 2 * (nx1 + ny1);

  /* Define distributions */
  PDM_g_num_t *distrib_vtx      = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_edge     = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_face     = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_edge_lim = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  distrib_vtx[0]      = 0;
  distrib_edge[0]     = 0;
  distrib_face[0]     = 0;
  distrib_edge_lim[0] = 0;

  PDM_g_num_t step_vtx       = *ng_vtx / n_rank;
  PDM_g_num_t remainder_vtx  = *ng_vtx % n_rank;

  PDM_g_num_t step_edge      = *ng_edge / n_rank;
  PDM_g_num_t remainder_edge = *ng_edge % n_rank;

  PDM_g_num_t step_face      = *ng_face / n_rank;
  PDM_g_num_t remainder_face = *ng_face % n_rank;

  PDM_g_num_t step_edge_lim      = ng_edge_lim / n_rank;
  PDM_g_num_t remainder_edge_lim = ng_edge_lim % n_rank;

  for (int i = 0; i < n_rank; i++) {
    distrib_vtx[i+1] = distrib_vtx[i] + step_vtx;
    if (i < remainder_vtx) {
      distrib_vtx[i+1]++;
    }

    distrib_edge[i+1] = distrib_edge[i] + step_edge;
    if (i < remainder_edge) {
      distrib_edge[i+1]++;
    }

    distrib_face[i+1] = distrib_face[i] + step_face;
    if (i < remainder_face) {
      distrib_face[i+1]++;
    }

    distrib_edge_lim[i+1] = distrib_edge_lim[i] + step_edge_lim;
    if (i < remainder_edge_lim) {
      distrib_edge_lim[i+1]++;
    }
  }
  *dn_vtx  = (int) distrib_vtx[i_rank+1]  - distrib_vtx[i_rank];
  *dn_edge = (int) distrib_edge[i_rank+1] - distrib_edge[i_rank];
  *dn_face = (int) distrib_face[i_rank+1] - distrib_face[i_rank];
  int dn_edge_lim = (int) (distrib_edge_lim[i_rank+1] - distrib_edge_lim[i_rank]);


  /*
   *  Vertices
   */
  *dvtx_coord = malloc (sizeof(double) * (*dn_vtx) * 3);
  double *_dvtx_coord = *dvtx_coord;

  double step_x = (xmax - xmin) / (double) nx1;
  double step_y = (ymax - ymin) / (double) ny1;

  PDM_g_num_t b_vtx_y = distrib_vtx[i_rank] / nx;
  PDM_g_num_t r_vtx_y = distrib_vtx[i_rank] % nx;

  int ivtx = 0;
  for (int j = b_vtx_y; j < ny; j++) {

    PDM_g_num_t _b_vtx_x = 0;
    if (j == b_vtx_y) {
      _b_vtx_x = r_vtx_y;
    }

    for (int i = _b_vtx_x; i < nx; i++) {
      _dvtx_coord[3*ivtx]     = xmin + step_x * (i + (i > 0 && i < nx1) * coef_rand*random01());
      _dvtx_coord[3*ivtx + 1] = ymin + step_y * (j + (j > 0 && j < ny1) * coef_rand*random01());
      _dvtx_coord[3*ivtx + 2] = 0.;
      ivtx++;
      if (ivtx == *dn_vtx) break;
    }
    if (ivtx == *dn_vtx) break;
  }

  /*
   *  Edges
   */
  *dedge_vtx  = malloc (sizeof(PDM_g_num_t ) * (*dn_edge) * 2);
  *dedge_face = malloc (sizeof(PDM_g_num_t ) * (*dn_edge) * 2);
  PDM_g_num_t  *_dedge_vtx = *dedge_vtx;
  PDM_g_num_t  *_dedge_face = *dedge_face;

  PDM_g_num_t ng_edge_h = nx1 * ny;

  int iedg = 0;
  // Horizontal edges
  if (distrib_edge[i_rank] < ng_edge_h) {
    const PDM_g_num_t b_edge_y = distrib_edge[i_rank] / nx1;
    const PDM_g_num_t r_edge_y = distrib_edge[i_rank] % nx1;

    for (PDM_g_num_t j = b_edge_y; j < ny; j++) {

      PDM_g_num_t _b_edge_x = 0;
      if (j == b_edge_y) {
        _b_edge_x = r_edge_y;
      }

      for (PDM_g_num_t i = _b_edge_x; i < nx1; i++) {
        if (j < ny1) {
          _dedge_vtx[2*iedg    ] = 1 + i   + nx*j;
          _dedge_vtx[2*iedg + 1] = 1 + i+1 + nx*j;

          _dedge_face[2*iedg] = 1 + i + nx1*j;
          if (j == 0) {
            _dedge_face[2*iedg + 1] = 0;
          } else {
            _dedge_face[2*iedg + 1] = 1 + i + nx1*(j-1);
          }
        } else {
          _dedge_vtx[2*iedg    ] = 1 + i+1 + nx*j;
          _dedge_vtx[2*iedg + 1] = 1 + i   + nx*j;

          _dedge_face[2*iedg    ] = 1 + i + nx1*(j-1);
          _dedge_face[2*iedg + 1] = 0;
        }

        iedg++;
        if (iedg == *dn_edge) break;
      }
      if (iedg == *dn_edge) break;
    }
  }

  // Vertical edges
  if (iedg < *dn_edge) {
    const PDM_g_num_t b_edge_y = (distrib_edge[i_rank] + iedg - ng_edge_h) / nx;
    const PDM_g_num_t r_edge_y = (distrib_edge[i_rank] + iedg - ng_edge_h) % nx;

    for (PDM_g_num_t j = b_edge_y; j < ny1; j++) {

      PDM_g_num_t _b_edge_x = 0;
      if (j == b_edge_y) {
        _b_edge_x = r_edge_y;
      }

      for (PDM_g_num_t i = _b_edge_x; i < nx; i++) {
        if (i > 0) {
          _dedge_vtx[2*iedg    ] = 1 + i + nx*j;
          _dedge_vtx[2*iedg + 1] = 1 + i + nx*(j+1);

          _dedge_face[2*iedg] = 1 + i-1 + nx1*j;
          if (i < nx1) {
            _dedge_face[2*iedg + 1] = 1 + i + nx1*j;
          } else {
            _dedge_face[2*iedg + 1] = 0;
          }
        } else {
          _dedge_vtx[2*iedg    ] = 1 + i + nx*(j+1);
          _dedge_vtx[2*iedg + 1] = 1 + i + nx*j;

          _dedge_face[2*iedg    ] = 1 + i + nx1*j;
          _dedge_face[2*iedg + 1] = 0;
        }

        iedg++;
        if (iedg == *dn_edge) break;
      }
      if (iedg == *dn_edge) break;
    }
  }



  /*
   *  Edge lim
   */
  *dedge_group_idx = malloc (sizeof(int) * (*n_edge_group + 1));
  *dedge_group     = malloc (sizeof(PDM_g_num_t) * dn_edge_lim);
  int *_dedge_group_idx = *dedge_group_idx;
  PDM_g_num_t *_dedge_group = *dedge_group;

  _dedge_group_idx[0] = 0;
  iedg = 0;
  if (distrib_edge_lim[i_rank] < nx1) {
    for (PDM_g_num_t i = distrib_edge_lim[i_rank]; i < nx1; i++) {
      _dedge_group[iedg++] = 1 + i;
      if (iedg == dn_edge_lim) break;
    }
  }
  _dedge_group_idx[1] = iedg;


  if (iedg < dn_edge_lim && distrib_edge_lim[i_rank] < nx1 + ny1) {
    PDM_g_num_t b = PDM_MAX (distrib_edge_lim[i_rank] - nx1, 0);

    for (PDM_g_num_t j = b; j < ny1; j++) {
      _dedge_group[iedg++] = ng_edge_h + 1 + nx1 + nx*j;
      if (iedg == dn_edge_lim) break;
    }
  }
  _dedge_group_idx[2] = iedg;


  if (iedg < dn_edge_lim && distrib_edge_lim[i_rank] < 2*nx1 + ny1) {
    PDM_g_num_t b = PDM_MAX (distrib_edge_lim[i_rank] - nx1 - ny1, 0);

    for (PDM_g_num_t i = b; i < nx1; i++) {
      PDM_g_num_t _i = nx1 - i - 1;
      _dedge_group[iedg++] = 1 + _i + nx1*ny1;
      if (iedg == dn_edge_lim) break;
    }
  }
  _dedge_group_idx[3] = iedg;


  if (iedg < dn_edge_lim) {
    PDM_g_num_t b = PDM_MAX (distrib_edge_lim[i_rank] - 2*(nx1 + ny1), 0);

    for (PDM_g_num_t j = b; j < ny1; j++) {
      PDM_g_num_t _j = ny1 - j - 1;
      _dedge_group[iedg++] = ng_edge_h + 1 + nx*_j;
      if (iedg == dn_edge_lim) break;
    }
  }
  _dedge_group_idx[4] = iedg;


  /*
   *  Faces
   */
  *dface_vtx_idx = malloc (sizeof(int) * (*dn_face + 1));
  int *_dface_vtx_idx = *dface_vtx_idx;
  _dface_vtx_idx[0] = 0;

  for (int i = 0; i < *dn_face; i++) {
    _dface_vtx_idx[i+1] = _dface_vtx_idx[i] + 4;
  }

  *dface_vtx  = malloc (sizeof(PDM_g_num_t) * _dface_vtx_idx[*dn_face]);
  *dface_edge = malloc (sizeof(PDM_g_num_t) * _dface_vtx_idx[*dn_face]);
  PDM_g_num_t *_dface_vtx  = *dface_vtx;
  PDM_g_num_t *_dface_edge = *dface_edge;

  PDM_g_num_t b_face_y = distrib_face[i_rank] / nx1;
  PDM_g_num_t r_face_y = distrib_face[i_rank] % nx1;

  int ifac = 0;

  for (int j = b_face_y; j < ny1; j++) {

    PDM_g_num_t _b_face_x = 0;
    if (j == b_face_y) {
      _b_face_x = r_face_y;
    }

    for (int i = _b_face_x; i < nx1; i++) {
      _dface_vtx[4*ifac]     = 1 + i   + nx*j;
      _dface_vtx[4*ifac + 1] = 1 + i+1 + nx*j;
      _dface_vtx[4*ifac + 2] = 1 + i+1 + nx*(j+1);
      _dface_vtx[4*ifac + 3] = 1 + i   + nx*(j+1);

      _dface_edge[4*ifac]     = 1 + i + nx1*j;
      _dface_edge[4*ifac + 1] = ng_edge_h + 1 + i+1 + nx*j;
      _dface_edge[4*ifac + 2] = 1 + i + nx1*(j+1);
      _dface_edge[4*ifac + 3] = ng_edge_h + 1 + i + nx*j;

      if (j < ny1) _dface_edge[4*ifac + 2] = -_dface_edge[4*ifac + 2];
      if (i > 0)   _dface_edge[4*ifac + 3] = -_dface_edge[4*ifac + 3];

      ifac++;
      if (ifac == *dn_face) break;
    }
    if (ifac == *dn_face) break;
  }

}

/**
 *
 * \brief  Create and split Mesh
 *
 * \param [in]      nVtxSeg  Number of arguments
 * \param [in]      length   Lenght of square
 * \param [in]      n_part    Number to obtain on this processus
 * \param [in]      post     mesh export status
 * \param [in]      method   Split method
 *
 * \return ppartId  ppart identifier
 *
 */

static void
_create_split_mesh
(
 int               activeRank,
 PDM_MPI_Comm      pdm_mpi_comm,
 double            xmin,
 double            ymin,
 PDM_g_num_t       nVtxSeg,
 double            length,
 int               rotate,
 int               quad,
 int               n_part,
 PDM_part_split_t  method,
 int               haveRandom,
 int               initRandom,
 PDM_g_num_t      *nGFace,
 PDM_g_num_t      *nGVtx,
 int            **nFace,
 int            ***faceVtxIdx,
 int            ***faceVtx,
 PDM_g_num_t    ***faceLNToGN,
 int            **nVtx,
 double         ***vtxCoord,
 PDM_g_num_t    ***vtxLNToGN
)
{
  struct timeval t_elaps_debut;

  int i_rank;
  int numProcs;

  if (activeRank) {

    PDM_MPI_Comm_rank (pdm_mpi_comm, &i_rank);
    PDM_MPI_Comm_size (pdm_mpi_comm, &numProcs);

    double        xmax = xmin + length;
    double        ymax = ymin + length;
    PDM_g_num_t  nx = nVtxSeg;
    PDM_g_num_t  ny = nVtxSeg;

    int           dNFace;
    int           dNVtx;
    int           dNEdge;
    int          *dFaceVtxIdx;
    PDM_g_num_t *dFaceVtx;
    double       *dVtxCoord;
    PDM_g_num_t *dFaceEdge;
    PDM_g_num_t *dEdgeVtx;
    PDM_g_num_t *dEdgeFace;
    int           nEdgeGroup;
    int          *dEdgeGroupIdx;
    PDM_g_num_t   *dEdgeGroup;


    /*
     *  Create mesh i
     */

    gettimeofday(&t_elaps_debut, NULL);

    PDM_g_num_t nGEdge;

    if (quad) {
      _quad_mesh_gen (pdm_mpi_comm,
                      xmin,
                      xmax,
                      ymin,
                      ymax,
                      haveRandom,
                      initRandom,
                      nx,
                      ny,
                      nGFace,
                      nGVtx,
                      &nGEdge,
                      &dNVtx,
                      &dVtxCoord,
                      &dNFace,
                      &dFaceVtxIdx,
                      &dFaceVtx,
                      &dFaceEdge,
                      &dNEdge,
                      &dEdgeVtx,
                      &dEdgeFace,
                      &nEdgeGroup,
                      &dEdgeGroupIdx,
                      &dEdgeGroup);
    } else {
      PDM_poly_surf_gen (pdm_mpi_comm,
                         xmin,
                         xmax,
                         ymin,
                         ymax,
                         haveRandom,
                         initRandom,
                         nx,
                         ny,
                         nGFace,
                         nGVtx,
                         &nGEdge,
                         &dNVtx,
                         &dVtxCoord,
                         &dNFace,
                         &dFaceVtxIdx,
                         &dFaceVtx,
                         &dFaceEdge,
                         &dNEdge,
                         &dEdgeVtx,
                         &dEdgeFace,
                         &nEdgeGroup,
                         &dEdgeGroupIdx,
                         &dEdgeGroup);
    }

    if (1) {
      double angle = PDM_PI / 6.;
      double rmin = 0.3;
      double delta_r = 0.7;
      for (int i = 0; i < dNVtx; i++) {
        /*double x = dVtxCoord[3*i]   - 0.5;
        double y = dVtxCoord[3*i+1] - 0.5;
        double f = PDM_MAX (fabs(x), fabs(y)) / sqrt(x*x + y*y);
        x = 0.5 + f*x;
        y = 0.5 + f*y;
        dVtxCoord[3*i]   = x;
        dVtxCoord[3*i+1] = y;*/
        double x = dVtxCoord[3*i] - 0.5;
        double y = dVtxCoord[3*i+1];
        double t = angle * (x + 0.2*y);
        double r = rmin + delta_r * y;
        dVtxCoord[3*i]   = r * cos(t);
        dVtxCoord[3*i+1] = r * sin(t);
      }
    }

    if (rotate) {
      double rot[3][3];
      _rotation_matrix (rot);

      for (int i = 0; i < dNVtx; i++) {
        double x = dVtxCoord[3*i];
        double y = dVtxCoord[3*i+1];
        double z = dVtxCoord[3*i+2];

        for (int j = 0; j < 3; j++) {
          dVtxCoord[3*i+j] = x * rot[j][0] + y * rot[j][1] + z * rot[j][2];
        }
      }
    }

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

    /*
     *  Create mesh partitions
     */

    int have_dCellPart = 0;

    int *dCellPart = (int *) malloc (dNFace*sizeof(int));
    int *dEdgeVtxIdx = (int *) malloc ((dNEdge+1)*sizeof(int));

    dEdgeVtxIdx[0] = 0;
    for (int i = 0; i < dNEdge; i++) {
      dEdgeVtxIdx[i+1] = 2 + dEdgeVtxIdx[i];
    }

    /*
     *  Split mesh i
     */

    int ppartId;

    int nPropertyCell = 0;
    int *renum_properties_cell = NULL;
    int nPropertyFace = 0;
    int *renum_properties_face = NULL;

    // PDM_part_create (&ppartId,
    //                  pdm_mpi_comm,
    //                  method,
    //                  "PDM_PART_RENUM_CELL_NONE",
    //                  "PDM_PART_RENUM_FACE_NONE",
    //                  nPropertyCell,
    //                  renum_properties_cell,
    //                  nPropertyFace,
    //                  renum_properties_face,
    //                  n_part,
    //                  dNFace,
    //                  dNEdge,
    //                  dNVtx,
    //                  nEdgeGroup,
    //                  NULL,
    //                  NULL,
    //                  NULL,
    //                  NULL,
    //                  have_dCellPart,
    //                  dCellPart,
    //                  dEdgeFace,
    //                  dEdgeVtxIdx,
    //                  dEdgeVtx,
    //                  NULL,
    //                  dVtxCoord,
    //                  NULL,
    //                  dEdgeGroupIdx,
    //                  dEdgeGroup);

    printf("dNFace = %i | dNEdge = %i | dNVtx = %i \n", dNFace, dNEdge, dNVtx);
    PDM_part_create (&ppartId,
                     pdm_mpi_comm,
                     method,
                     "PDM_PART_RENUM_CELL_NONE",
                     "PDM_PART_RENUM_FACE_NONE",
                     nPropertyCell,
                     renum_properties_cell,
                     nPropertyFace,
                     renum_properties_face,
                     n_part,
                     dNFace,
                     dNEdge,
                     dNVtx,
                     nEdgeGroup,
                     dFaceVtxIdx,
                     dFaceEdge,
                     NULL,
                     NULL,
                     have_dCellPart,
                     dCellPart,
                     NULL,
                     dEdgeVtxIdx,
                     dEdgeVtx,
                     NULL,
                     dVtxCoord,
                     NULL,
                     dEdgeGroupIdx,
                     dEdgeGroup);
    free (dCellPart);

    double  *elapsed = NULL;
    double  *cpu = NULL;
    double  *cpu_user = NULL;
    double  *cpu_sys = NULL;

    PDM_part_time_get (ppartId,
                       &elapsed,
                       &cpu,
                       &cpu_user,
                       &cpu_sys);

    if (i_rank == 0)
      PDM_printf("[%d] Temps dans ppart : %12.5e\n",
                 i_rank,  elapsed[0]);

    /* Statistiques */

    int    cells_average;
    int    cells_median;
    double cells_std_deviation;
    int    cells_min;
    int    cells_max;
    int    bound_part_faces_average;
    int    bound_part_faces_median;
    double bound_part_faces_std_deviation;
    int    bound_part_faces_min;
    int    bound_part_faces_max;
    int    bound_part_faces_sum;

    PDM_part_stat_get (ppartId,
                       &cells_average,
                       &cells_median,
                       &cells_std_deviation,
                       &cells_min,
                       &cells_max,
                       &bound_part_faces_average,
                       &bound_part_faces_median,
                       &bound_part_faces_std_deviation,
                       &bound_part_faces_min,
                       &bound_part_faces_max,
                       &bound_part_faces_sum);

    if (i_rank == 0) {
      PDM_printf ("Statistics :\n");
      PDM_printf ("  - Number of cells :\n");
      PDM_printf ("       * average            : %i\n", cells_average);
      PDM_printf ("       * median             : %i\n", cells_median);
      PDM_printf ("       * standard deviation : %12.5e\n", cells_std_deviation);
      PDM_printf ("       * min                : %i\n", cells_min);
      PDM_printf ("       * max                : %i\n", cells_max);
      PDM_printf ("  - Number of faces exchanging with another partition :\n");
      PDM_printf ("       * average            : %i\n", bound_part_faces_average);
      PDM_printf ("       * median             : %i\n", bound_part_faces_median);
      PDM_printf ("       * standard deviation : %12.5e\n", bound_part_faces_std_deviation);
      PDM_printf ("       * min                : %i\n", bound_part_faces_min);
      PDM_printf ("       * max                : %i\n", bound_part_faces_max);
      PDM_printf ("       * total              : %i\n", bound_part_faces_sum);
    }

    free (dVtxCoord);
    free (dFaceVtxIdx);
    free (dFaceVtx);
    free (dFaceEdge);
    free (dEdgeVtxIdx);
    free (dEdgeVtx);
    free (dEdgeFace);
    free (dEdgeGroupIdx);
    free (dEdgeGroup);

    _compute_faceVtx (ppartId,
                      n_part,
                      nFace,
                      faceVtxIdx,
                      faceVtx,
                      faceLNToGN,
                      nVtx,
                      vtxCoord,
                      vtxLNToGN);

    PDM_part_free (ppartId);

  }
  else {
    *nFace = (int *) malloc(sizeof(int) * n_part);
    *faceVtxIdx = (int **) malloc(sizeof(int *) * n_part);
    *faceVtx = (int **) malloc(sizeof(int *) * n_part);
    *faceLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    *nVtx = (int *) malloc(sizeof(int) * n_part);
    *vtxCoord = (double **) malloc(sizeof(double *) * n_part);
    *vtxLNToGN = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

    for (int ipart = 0; ipart < n_part; ipart++) {
      int _nFace = 0;
      int _sFaceEdge = 0;
      (*nFace)[ipart] = _nFace;
      (*faceVtxIdx) [ipart] = (int *) malloc(sizeof(int) * (_nFace + 1));
      (*faceVtxIdx) [ipart][0] = 0;
      (*faceVtx)[ipart] = (int *) malloc(sizeof(int) * _sFaceEdge);
      (*faceLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nFace);

      int _nVtx = 0;
      (*nVtx)[ipart] = _nVtx;
      (*vtxCoord)[ipart] = (double *) malloc(sizeof(double) * (3 * _nVtx));
      (*vtxLNToGN)[ipart] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);
    }
  }

  PDM_MPI_Bcast (nGFace, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);
  PDM_MPI_Bcast (nGVtx, 1, PDM__PDM_MPI_G_NUM, 0, PDM_MPI_COMM_WORLD);

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
        val_match[ipart][olLinkedFace[4*i]-1] = i;//100;
        val_cell_match[ipart][olLinkedFace[4*i]-1] = olLinkedFace[4*i+3];//xolLinkedFace[4*i + 2];
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

  int              post    = 0;
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
  int              haveRandom = 1;
  int              randomTimeInit = 0;

  int              randomMeshAInit = -1;
  int              randomMeshBInit = -1;

  int              nProcData = -1;

  int              i_rank;
  int              numProcs;

  int rotate = 0;
  int quadA = 0;
  int quadB = 0;
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
              &post,
              (int *) &method,
              &haveRandom,
              &randomTimeInit,
              &randomMeshAInit,
              &randomMeshBInit,
              &quadA,
              &quadB,
              &nProcData,
              &rotate);

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
    int quad;

    if (imesh == 0) {
      n_vtx_seg = n_vtx_segA;
      length = lengthA;
      xmin = xminA;
      ymin = yminA;
      n_part = n_partA;
      quad = quadA;
    }
    else {
      n_vtx_seg = n_vtx_segB;
      length = lengthB;
      xmin = xminB;
      ymin = yminB;
      n_part = n_partB;
      quad = quadB;
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

    _create_split_mesh (activeRankMesh,
                        meshComm,
                        xmin,
                        ymin,
                        n_vtx_seg,
                        length,
                        rotate,
                        quad,
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

  if (nProcData > 0 && nProcData < numProcs) {
    PDM_MPI_Comm_free(&meshComm);
  }

  if (1 == 0) {
    for (int imesh = 0; imesh < 2; imesh++) {
      for (int ipart = 0; ipart < n_part; ipart++) {
        printf("vtxCoord :\n");
        for (int j = 0; j < nVtx[imesh][ipart]; j++) {
          printf(""PDM_FMT_G_NUM" : %12.5e %12.5e %12.5e\n", vtxLNToGN[imesh][ipart][j],
                 vtxCoord[imesh][ipart][3*j],
                 vtxCoord[imesh][ipart][3*j+1],
                 vtxCoord[imesh][ipart][3*j+2]);
        }
        printf("faceVtx :\n");
        for (int j = 0; j < nFace[imesh][ipart]; j++) {
          printf(""PDM_FMT_G_NUM" :", faceLNToGN[imesh][ipart][j]);
          for (int k = faceVtxIdx[imesh][ipart][j]; k < faceVtxIdx[imesh][ipart][j+1]; k++) {
            printf(" %d", faceVtx[imesh][ipart][k]);
          }
          printf("\n");
        }
      }
    }
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

  /*
   *  Free Pdm
   */

  PDM_ol_del (pdm_id);

  PDM_MPI_Finalize ();

  return 0;

}
