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

#include "pdm_predicate.h"

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
    else if (strcmp (argv[i], "-hilbert") == 0) {
      *method = 3;
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
 PDM_part_t    *ppart,
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

  PDM_malloc(*nFace,n_part,int);
  PDM_malloc(*faceVtxIdx,n_part,int *);
  PDM_malloc(*faceVtx,n_part,int *);
  PDM_malloc(*faceLNToGN,n_part,PDM_g_num_t *);

  PDM_malloc(*nVtx,n_part,int);
  PDM_malloc(*vtxCoord,n_part,double *);
  PDM_malloc(*vtxLNToGN,n_part,PDM_g_num_t *);

  // int id_ppart = ppartId;

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

    PDM_part_part_dim_get (ppart,
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

    PDM_part_part_val_get (ppart,
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
    (*faceVtxIdx) PDM_malloc([ipart],(_nFace + 1),int);
    ( *faceVtx)[ipart];
    PDM_malloc(faceVtx)[ipart],_sFaceEdge,int);
    ( *faceLNToGN)[ipart];
    PDM_malloc(faceLNToGN)[ipart],_nFace,PDM_g_num_t);

    memcpy ((*faceVtxIdx)[ipart], _faceEdgeIdx, (_nFace + 1) * sizeof(int));
    memcpy ((*faceLNToGN)[ipart], _faceLNToGN, _nFace * sizeof(PDM_g_num_t));

    (*nVtx)[ipart] = _nVtx;
    ( *vtxCoord)[ipart];
    PDM_malloc(vtxCoord)[ipart],(3 * _nVtx),double);
    ( *vtxLNToGN)[ipart];
    PDM_malloc(vtxLNToGN)[ipart],_nVtx,PDM_g_num_t);

    memcpy ((*vtxCoord)[ipart], _vtx, 3 *_nVtx * sizeof(double));
    memcpy ((*vtxLNToGN)[ipart], _vtxLNToGN, _nVtx * sizeof(PDM_g_num_t));

    int *_faceVtx = (*faceVtx)[ipart];

    int *vtxEdgeIdx;
    PDM_malloc(vtxEdgeIdx,(_nVtx + 1),int);

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

    int *vtxEdge;
    PDM_malloc(vtxEdge,vtxEdgeIdx[_nVtx],int);
    int *vtxEdgeN;
    PDM_malloc(vtxEdgeN,_nVtx,int);
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
   PDM_free(vtxEdgeN);

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

   PDM_free(vtxEdge);
   PDM_free(vtxEdgeIdx);

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
 int               n_part,
 PDM_part_split_t   method,
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

    int *dCellPart;
    PDM_malloc(dCellPart,dNFace,int);
    int *dEdgeVtxIdx;
    PDM_malloc(dEdgeVtxIdx,(dNEdge+1),int);

    dEdgeVtxIdx[0] = 0;
    for (int i = 0; i < dNEdge; i++) {
      dEdgeVtxIdx[i+1] = 2 + dEdgeVtxIdx[i];
    }

    /*
     *  Split mesh i
     */

    // int ppartId;

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

    //printf("dNFace = %i | dNEdge = %i | dNVtx = %i \n", dNFace, dNEdge, dNVtx);
    PDM_part_t *ppart = PDM_part_create (pdm_mpi_comm,
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
   PDM_free(dCellPart);

    double  *elapsed = NULL;
    double  *cpu = NULL;
    double  *cpu_user = NULL;
    double  *cpu_sys = NULL;

    PDM_part_time_get (ppart,
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

    PDM_part_stat_get (ppart,
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

   PDM_free(dVtxCoord);
   PDM_free(dFaceVtxIdx);
   PDM_free(dFaceVtx);
   PDM_free(dFaceEdge);
   PDM_free(dEdgeVtxIdx);
   PDM_free(dEdgeVtx);
   PDM_free(dEdgeFace);
   PDM_free(dEdgeGroupIdx);
   PDM_free(dEdgeGroup);

    _compute_faceVtx (ppart,
                      n_part,
                      nFace,
                      faceVtxIdx,
                      faceVtx,
                      faceLNToGN,
                      nVtx,
                      vtxCoord,
                      vtxLNToGN);

    PDM_part_free (ppart);

  }
  else {
    PDM_malloc(*nFace,n_part,int);
    PDM_malloc(*faceVtxIdx,n_part,int *);
    PDM_malloc(*faceVtx,n_part,int *);
    PDM_malloc(*faceLNToGN,n_part,PDM_g_num_t *);

    PDM_malloc(*nVtx,n_part,int);
    PDM_malloc(*vtxCoord,n_part,double *);
    PDM_malloc(*vtxLNToGN,n_part,PDM_g_num_t *);

    for (int ipart = 0; ipart < n_part; ipart++) {
      int _nFace = 0;
      int _sFaceEdge = 0;
      (*nFace)[ipart] = _nFace;
      (*faceVtxIdx) PDM_malloc([ipart],(_nFace + 1),int);
      (*faceVtxIdx) [ipart][0] = 0;
      ( *faceVtx)[ipart];
      PDM_malloc(faceVtx)[ipart],_sFaceEdge,int);
      ( *faceLNToGN)[ipart];
      PDM_malloc(faceLNToGN)[ipart],_nFace,PDM_g_num_t);

      int _nVtx = 0;
      (*nVtx)[ipart] = _nVtx;
      ( *vtxCoord)[ipart];
      PDM_malloc(vtxCoord)[ipart],(3 * _nVtx),double);
      ( *vtxLNToGN)[ipart];
      PDM_malloc(vtxLNToGN)[ipart],_nVtx,PDM_g_num_t);
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

  PDM_writer_t *id_cs[2];

  id_cs[0] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CST,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "mesh1",
                                pdm_mpi_comm,
                                PDM_IO_KIND_MPI_SIMPLE,
                                1.,
                                NULL);

  id_cs[1] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CST,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "mesh2",
                                pdm_mpi_comm,
                                PDM_IO_KIND_MPI_SIMPLE,
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
                                                    PDM_WRITER_VAR_SCALAR,
                                                    PDM_WRITER_VAR_ELEMENTS,
                                                    "num_part");
    if (imesh == 0) {

      id_var_field[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                   PDM_WRITER_OFF,
                                                   PDM_WRITER_VAR_SCALAR,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "sfieldA");
    }
    else {

      id_var_field[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                   PDM_WRITER_OFF,
                                                   PDM_WRITER_VAR_SCALAR,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "rfieldB");
    }


    id_var_coo_x[imesh] = PDM_writer_var_create (id_cs[imesh],
                                       PDM_WRITER_ON,
                                       PDM_WRITER_VAR_SCALAR,
                                       PDM_WRITER_VAR_VERTICES,
                                       "coo_x");

    id_var_coo_xyz[imesh] = PDM_writer_var_create (id_cs[imesh],
                                         PDM_WRITER_ON,
                                         PDM_WRITER_VAR_VECTOR,
                                         PDM_WRITER_VAR_VERTICES,
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
                                             n_part);
    /*
     * Debut des ecritures
     */

    int *nsom_part;
    PDM_malloc(nsom_part,n_part,int);

    int *n_part_procs;
    PDM_malloc(n_part_procs,numProcs,int);

    PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
                   (void *) n_part_procs, 1, PDM_MPI_INT,
                   PDM_MPI_COMM_WORLD);

    int *debPartProcs;
    PDM_malloc(debPartProcs,(numProcs + 1),int);

    debPartProcs[0] = 0;
    for (int i = 0; i < numProcs; i++) {
      debPartProcs[i+1] = debPartProcs[i] + n_part_procs[i];
    }

   PDM_free(n_part_procs);

    PDM_writer_step_beg (id_cs[imesh], 0.);

    int **_face_nb;
    PDM_malloc(*_face_nb,n_part,int *);
    int **_face_idx;
    PDM_malloc(*_face_idx,n_part,int *);

    for (int ipart = 0; ipart < n_part; ipart++) {

      PDM_writer_geom_coord_set (id_cs[imesh],
                                 id_geom[imesh],
                                 ipart,
                                 nVtx[imesh][ipart],
                                 vtxCoord[imesh][ipart],
                                 vtxLNToGN[imesh][ipart],
                                PDM_OWNERSHIP_USER);

      PDM_malloc(_face_nb[ipart],nFace[imesh][ipart],int);
      PDM_malloc(_face_idx[ipart],nFace[imesh][ipart],int);

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
     PDM_free(_face_nb[ipart]);
     PDM_free(_face_idx[ipart]);
    }

   PDM_free(_face_nb);
   PDM_free(_face_idx);

    PDM_writer_geom_write(id_cs[imesh],
                          id_geom[imesh]);

    /* Creation des variables :
       - numero de partition
       - scalaire
       - vecteur
       - tenseur
    */

    PDM_real_t **val_num_part;
    PDM_malloc(*val_num_part,n_part,PDM_real_t *);
    PDM_real_t **val_coo_x;
    PDM_malloc(*val_coo_x,n_part,PDM_real_t *);
    PDM_real_t **val_coo_xyz;
    PDM_malloc(*val_coo_xyz,n_part,PDM_real_t *);

    for (int ipart = 0; ipart < n_part; ipart++) {

      PDM_malloc(val_num_part[ipart],nFace[imesh][ipart],PDM_real_t);
      PDM_malloc(val_coo_x[ipart],nVtx[imesh][ipart],PDM_real_t);
      PDM_malloc(val_coo_xyz[ipart],3 * nVtx[imesh][ipart],PDM_real_t);
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
     PDM_free(val_num_part[ipart]);
     PDM_free(val_coo_x[ipart]);
     PDM_free(val_coo_xyz[ipart]);
    }

   PDM_free(val_num_part);
   PDM_free(val_coo_x);
   PDM_free(val_coo_xyz);
   PDM_free(nsom_part);

    PDM_writer_step_end (id_cs[imesh]);
    PDM_writer_geom_data_free (id_cs[imesh],
                      id_geom[imesh]);

    PDM_writer_geom_free (id_cs[imesh],
                 id_geom[imesh]);
    PDM_writer_free (id_cs[imesh]);

   PDM_free(debPartProcs);

  }

}

/**
 *
 * \brief  Export overlay mesh
 *
 * \param [in]    ol    Pointer to overlay structure
 *
 */

static void
_export_ol_mesh
(
 const PDM_MPI_Comm pdm_mpi_comm,
 const PDM_ol_t *ol,
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

  PDM_writer_t *id_cs[2];

  id_cs[0] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CST,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "olmesh1",
                                pdm_mpi_comm,
                                PDM_IO_KIND_MPI_SIMPLE,
                                1.,
                                NULL);

  id_cs[1] = PDM_writer_create ("Ensight",
                                PDM_WRITER_FMT_ASCII,
                                PDM_WRITER_TOPO_CST,
                                PDM_WRITER_OFF,
                                "test_2d_surf_ens",
                                "olmesh2",
                                pdm_mpi_comm,
                                PDM_IO_KIND_MPI_SIMPLE,
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
                                                   PDM_WRITER_VAR_SCALAR,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "sOlField");
    }
    else {

      id_var_field[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                   PDM_WRITER_OFF,
                                                   PDM_WRITER_VAR_SCALAR,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "rOlField");
    }

    id_var_num_part[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                    PDM_WRITER_OFF,
                                                    PDM_WRITER_VAR_SCALAR,
                                                    PDM_WRITER_VAR_ELEMENTS,
                                                    "num_part");

    id_var_match[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                 PDM_WRITER_OFF,
                                                 PDM_WRITER_VAR_SCALAR,
                                                 PDM_WRITER_VAR_ELEMENTS,
                                                 "matching");

    id_var_cell_match[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                      PDM_WRITER_OFF,
                                                      PDM_WRITER_VAR_SCALAR,
                                                      PDM_WRITER_VAR_ELEMENTS,
                                                      "cell_matching");

    id_var_origin[imesh] = PDM_writer_var_create (id_cs[imesh],
                                                  PDM_WRITER_OFF,
                                                  PDM_WRITER_VAR_SCALAR,
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
                                             n_part);
    int *n_part_procs;
    PDM_malloc(n_part_procs,numProcs,int);

    PDM_MPI_Allgather ((void *) &n_part,      1, PDM_MPI_INT,
                   (void *) n_part_procs, 1, PDM_MPI_INT,
                   PDM_MPI_COMM_WORLD);

    int *debPartProcs;
    PDM_malloc(debPartProcs,(numProcs + 1),int);

    debPartProcs[0] = 0;
    for (int i = 0; i < numProcs; i++) {
      debPartProcs[i+1] = debPartProcs[i] + n_part_procs[i];
    }

   PDM_free(n_part_procs);

    /*
     * Debut des ecritures
     */

    PDM_g_num_t    nGOlFace;
    PDM_g_num_t    nGOlVtx;

    PDM_ol_mesh_dim_get (ol,
                         mesht,
                         &nGOlFace,
                         &nGOlVtx);

    int **_olface_nb;
    PDM_malloc(*_olface_nb,n_part,int *);
    int **_olface_idx;
    PDM_malloc(*_olface_idx,n_part,int *);
    PDM_real_t **val_num_part;
    PDM_malloc(*val_num_part,n_part,PDM_real_t *);
    PDM_real_t **val_match;
    PDM_malloc(*val_match,n_part,PDM_real_t *);
    PDM_real_t **val_cell_match;
    PDM_malloc(*val_cell_match,n_part,PDM_real_t *);
    PDM_real_t **val_origin;
    PDM_malloc(*val_origin,n_part,PDM_real_t *);
    PDM_writer_step_beg (id_cs[imesh], 0.);

    for (int ipart = 0; ipart < n_part; ipart++) {

      int           nOlFace;
      int           nOlLinkedFace;
      int           nOlVtx;
      int           sOlFaceIniVtx;
      int           sOlface_vtx;
      int           sInitToOlFace;

      PDM_ol_part_mesh_dim_get (ol,
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

      PDM_ol_mesh_entities_get (ol,
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


      PDM_malloc(val_num_part[ipart],nOlFace,PDM_real_t);
      PDM_malloc(val_match[ipart],nOlFace,PDM_real_t);
      PDM_malloc(val_cell_match[ipart],nOlFace,PDM_real_t);
      PDM_malloc(val_origin[ipart],nOlFace,PDM_real_t);

      PDM_writer_geom_coord_set (id_cs[imesh],
                                 id_geom[imesh],
                                 ipart,
                                 nOlVtx,
                                 olCoords,
                                 olvtx_ln_to_gn,
                                PDM_OWNERSHIP_USER);

      PDM_malloc(_olface_nb[ipart],nOlFace,int);
      PDM_malloc(_olface_idx[ipart],nOlFace,int);

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
     PDM_free(_olface_nb[ipart]);
     PDM_free(_olface_idx[ipart]);
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
     PDM_free(val_num_part[ipart]);
     PDM_free(val_match[ipart]);
     PDM_free(val_cell_match[ipart]);
     PDM_free(val_origin[ipart]);
    }

   PDM_free(val_num_part);
   PDM_free(val_match);
   PDM_free(val_cell_match);
   PDM_free(val_origin);

    PDM_writer_step_end (id_cs[imesh]);
    PDM_writer_geom_data_free (id_cs[imesh],
                               id_geom[imesh]);

    PDM_writer_geom_free (id_cs[imesh],
                          id_geom[imesh]);
    PDM_writer_free (id_cs[imesh]);

   PDM_free(_olface_nb);
   PDM_free(_olface_idx);
   PDM_free(debPartProcs);
  }

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
  PDM_predicate_exactinit();

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
  PDM_part_split_t method  = PDM_PART_SPLIT_HILBERT;
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
              &post,
              (int *) &method,
              &haveRandom,
              &randomTimeInit,
              &randomMeshAInit,
              &randomMeshBInit,
              &nProcData
              );

  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_rank         : %d\n", numProcs);
    PDM_printf ("  - n_vtx_segA     : %d\n", n_vtx_segA);
    PDM_printf ("  - lengthA        : %f\n", lengthA);
    PDM_printf ("  - xminA          : %f\n", xminA);
    PDM_printf ("  - yminA          : %f\n", yminA);
    PDM_printf ("  - n_partA        : %d\n", n_partA);
    PDM_printf ("  - n_vtx_segB     : %d\n", n_vtx_segB);
    PDM_printf ("  - lengthB        : %f\n", lengthB);
    PDM_printf ("  - xminB          : %f\n", xminB);
    PDM_printf ("  - yminB          : %f\n", yminB);
    PDM_printf ("  - n_partB        : %d\n", n_partB);
    PDM_printf ("  - post           : %d\n", post);
    PDM_printf ("  - method         : %d\n", method);
    PDM_printf ("  - haveRandom     : %d\n", haveRandom);
    PDM_printf ("  - randomTimeInit : %d\n", randomTimeInit);
    PDM_printf ("  - nProcData      : %d\n", nProcData);
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

    int *rankInNodes;
    PDM_malloc(rankInNodes,numProcs,int);

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

   PDM_free(rankInNodes);
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

    _create_split_mesh (activeRankMesh,
                        meshComm,
                        xmin,
                        ymin,
                        n_vtx_seg,
                        length,
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

  /*
   *  Appel des fonctions d'intersection
   */

  double projectCoeff = 0.;

  /*
   *  Creation de l'objet PDM
   */

  PDM_ol_t *ol = PDM_ol_create (n_part,
                                n_part,
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

  PDM_ol_parameter_set (ol,
                        PDM_OL_CAR_LENGTH_TOL,
                        1e-4);

  PDM_ol_parameter_set (ol,
                        PDM_OL_EXTENTS_TOL,
                        1e-4);

  /*
   *  Create field
   */

  double **sFieldA;
  PDM_malloc(*sFieldA,n_part,double *);
  double **centerA;
  PDM_malloc(*centerA,n_part,double *);
  double **rFieldB;
  PDM_malloc(*rFieldB,n_part,double *);

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_malloc(centerA[ipart],3 * nFace[0][ipart],double);
    for (int i = 0; i < 3*nFace[0][ipart]; i++) {
      centerA[ipart][i] = 0.;
    }
    PDM_malloc(sFieldA[ipart],nFace[0][ipart],double);
    PDM_malloc(rFieldB[ipart],nFace[1][ipart],double);
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

      PDM_ol_input_mesh_set (ol,
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


  /*if (post) {
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
		      }*/

  /*
   *  Calcul
   */

  PDM_ol_compute (ol);

  if (i_rank == 0){
    PDM_ol_dump_times (ol);
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

  double **sFieldOlA;
  PDM_malloc(*sFieldOlA,n_part,double *);
  double **rFieldOlB;
  PDM_malloc(*rFieldOlB,n_part,double *);
  double **surfB;
  PDM_malloc(*surfB,n_part,double *);
  double **surfOlB;
  PDM_malloc(*surfOlB,n_part,double *);

  for (int imesh = 0; imesh < 2; imesh++) {
    PDM_malloc(_nOlFace[imesh],n_part,int);
    PDM_malloc(_nOlLinkedFace[imesh],n_part,int);
    PDM_malloc(_nOlVtx[imesh],n_part,int);
    PDM_malloc(_sOlface_vtx[imesh],n_part,int);
    PDM_malloc(_sInitToOlFace[imesh],n_part,int);
    PDM_malloc(_olFaceIniVtxIdx[imesh],n_part,int*);
    PDM_malloc(_olFaceIniVtx[imesh],n_part,int*);
    PDM_malloc(_olface_vtx_idx[imesh],n_part,int*);
    PDM_malloc(_olface_vtx[imesh],n_part,int*);
    PDM_malloc(_olLinkedface_procIdx[imesh],n_part,int*);
    PDM_malloc(_olLinkedFace[imesh],n_part,int*);
    PDM_malloc(_olface_ln_to_gn[imesh],n_part,PDM_g_num_t*);
    PDM_malloc(_olCoords[imesh],n_part,double*);;
    PDM_malloc(_olvtx_ln_to_gn[imesh],n_part,PDM_g_num_t*);;
    PDM_malloc(_initToOlFaceIdx[imesh],n_part,int*);
    PDM_malloc(_initToOlFace[imesh],n_part,int*);
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

    PDM_ol_mesh_dim_get (ol,
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

      PDM_ol_part_mesh_dim_get (ol,
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

      PDM_ol_mesh_entities_get (ol,
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
        PDM_malloc(sFieldOlA[ipart],nOlFace,double);

        for (int i = 0; i < nFace[imesh][ipart]; i++) {
          for (int j = initToOlFaceIdx[i]; j < initToOlFaceIdx[i+1]; j++) {
            sFieldOlA[ipart][initToOlFace[j]-1] = sFieldA[ipart][i];
          }
        }
      }

      if (imesh == 1) {
        PDM_malloc(rFieldOlB[ipart],nOlFace,double);

        PDM_malloc(surfB[ipart],nFace[imesh][ipart],double);
        PDM_malloc(surfOlB[ipart],nOlFace,double);

        double *ol_surface_vector;
        PDM_malloc(ol_surface_vector,3 * nOlFace,double);
        double *ol_center;
        PDM_malloc(ol_center,3 * nOlFace,double);
        double *ol_characteristicLength;
        PDM_malloc(ol_characteristicLength,nOlFace,double);
        int *ol_isDegenerated;
        PDM_malloc(ol_isDegenerated,nOlFace,int);

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

       PDM_free(ol_surface_vector);
       PDM_free(ol_center);
       PDM_free(ol_characteristicLength);
       PDM_free(ol_isDegenerated);
      }
    }
  }

  int *n_send;
  PDM_malloc(n_send,numProcs,int);
  int *n_recv;
  PDM_malloc(n_recv,numProcs,int);
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

  int *i_send;
  PDM_malloc(i_send,(numProcs+1),int);
  int *i_recv;
  PDM_malloc(i_recv,(numProcs+1),int);

  i_send[0] = 0;
  i_recv[0] = 0;
  for (int iproc = 0; iproc < numProcs; iproc++) {
    i_send[iproc+1] = i_send[iproc] + 4 * n_send[iproc];
    i_recv[iproc+1] = i_recv[iproc] + 4 * n_recv[iproc];
    n_send[iproc] = 0;
    n_recv[iproc] *= 4;
  }

  int *b_send;
  PDM_malloc(b_send,i_send[numProcs],int);
  int *b_recv;
  PDM_malloc(b_recv,i_recv[numProcs],int);

  double *d_send;
  PDM_malloc(d_send,i_send[numProcs],double);
  double *d_recv;
  PDM_malloc(d_recv,i_recv[numProcs],double);

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

 PDM_free(n_send);
 PDM_free(b_send);
 PDM_free(i_send);
 PDM_free(n_recv);
 PDM_free(d_send);

  int **check_graph;
  PDM_malloc(*check_graph,n_part,int*);

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_malloc(check_graph[ipart],2 * _nOlFace[1][ipart],int);
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

 PDM_free(d_recv);

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
   PDM_free(check_graph[ipart]);
  }
 PDM_free(check_graph);

 PDM_free(i_recv);
 PDM_free(b_recv);

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
                     ol,
                     nFace,
                     sFieldOlA,
                     rFieldOlB,
                     n_part);
  }

  /*
   *  Free meshes
   */

  for (int ipart = 0; ipart < n_part; ipart++) {
   PDM_free(sFieldOlA[ipart]);
   PDM_free(sFieldA[ipart]);
   PDM_free(centerA[ipart]);
   PDM_free(rFieldOlB[ipart]);
   PDM_free(rFieldB[ipart]);
   PDM_free(surfB[ipart]);
   PDM_free(surfOlB[ipart]);
  }

 PDM_free(sFieldA);
 PDM_free(sFieldOlA);
 PDM_free(centerA);
 PDM_free(rFieldB);
 PDM_free(rFieldOlB);
 PDM_free(surfB);
 PDM_free(surfOlB);

  for (int imesh = 0; imesh < 2; imesh++) {

    for (int ipart = 0; ipart < n_part; ipart++) {
     PDM_free(faceVtxIdx[imesh][ipart]);
     PDM_free(faceVtx[imesh][ipart]);
     PDM_free(faceLNToGN[imesh][ipart]);
     PDM_free(vtxCoord[imesh][ipart]);
     PDM_free(vtxLNToGN[imesh][ipart]);
     PDM_free(_olFaceIniVtxIdx[imesh][ipart]);
     PDM_free(_olFaceIniVtx[imesh][ipart]);
     PDM_free(_olface_vtx_idx[imesh][ipart]);
     PDM_free(_olface_vtx[imesh][ipart]);
     PDM_free(_olLinkedface_procIdx[imesh][ipart]);
     PDM_free(_olLinkedFace[imesh][ipart]);
     PDM_free(_olface_ln_to_gn[imesh][ipart]);
     PDM_free(_olCoords[imesh][ipart]);
     PDM_free(_olvtx_ln_to_gn[imesh][ipart]);
     PDM_free(_initToOlFaceIdx[imesh][ipart]);
     PDM_free(_initToOlFace[imesh][ipart]);
    }

   PDM_free(faceVtxIdx[imesh]);
   PDM_free(faceVtx[imesh]);
   PDM_free(faceLNToGN[imesh]);
   PDM_free(vtxCoord[imesh]);
   PDM_free(vtxLNToGN[imesh]);
   PDM_free(nFace[imesh]);
   PDM_free(nVtx[imesh]);
   PDM_free(_olFaceIniVtxIdx[imesh]);
   PDM_free(_olFaceIniVtx[imesh]);
   PDM_free(_olface_vtx_idx[imesh]);
   PDM_free(_olface_vtx[imesh]);
   PDM_free(_olLinkedface_procIdx[imesh]);
   PDM_free(_olLinkedFace[imesh]);
   PDM_free(_olface_ln_to_gn[imesh]);
   PDM_free(_olCoords[imesh]);
   PDM_free(_olvtx_ln_to_gn[imesh]);
   PDM_free(_initToOlFaceIdx[imesh]);
   PDM_free(_initToOlFace[imesh]);

   PDM_free(_nOlFace[imesh]);
   PDM_free(_nOlLinkedFace[imesh]);
   PDM_free(_nOlVtx[imesh]);
   PDM_free(_sOlface_vtx[imesh]);
   PDM_free(_sInitToOlFace[imesh]);
  }

  /*
   *  Free Pdm
   */

  PDM_ol_del (ol);

  if (i_rank == 0) printf("-- End\n");

  PDM_MPI_Finalize ();

  return 0;

}
