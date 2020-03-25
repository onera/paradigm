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
#include "pdm_dcube_gen.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"

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
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -p      <level>  Number of points to locate.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *nVtxSeg,
           double        *length,
           int           *n_part,
	   PDM_g_num_t   *nPts,
           int           *post,
           int           *method)
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
        long _nVtxSeg = atol(argv[i]);
        *nVtxSeg = (PDM_g_num_t) _nVtxSeg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
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
        *nPts = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = PDM_PART_SPLIT_HILBERT;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}




static void
_gen_cloud_random
(
 const int      nPts,
 const double   length,
 const int      numProcs,
 const int      myRank,
 double       **coord,
 int           *nPts_l
 )
{
  *nPts_l = (int) (nPts/numProcs);
  *coord = malloc (sizeof(double) * 3 * (*nPts_l));
  double *_coord = *coord;
  double x;
  int idx = 0;
  for (PDM_g_num_t i = 0; i < numProcs*(*nPts_l); i++) {
    for (int j = 0; j < 3; j++) {
      x = length * (double) rand() / ((double) RAND_MAX);
      if (i%numProcs == myRank) {
        _coord[idx++] = x;
      }
    }
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

  PDM_g_num_t nVtxSeg = 10;
  double       length = 1.;
  int          nPart  = 1;
  int          post   = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#else
  PDM_part_split_t method  = PDM_PART_SPLIT_HILBERT;
#endif
#endif

  PDM_g_num_t nPts = 10;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nVtxSeg,
             &length,
             &nPart,
	     &nPts,
             &post,
             (int *) &method);


/*
   *  Init
   */

  struct timeval t_elaps_debut;

  int myRank;
  int numProcs;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &numProcs);

  int           dNCell;
  int           dNFace;
  int           dNVtx;
  int           nFaceGroup;
  PDM_g_num_t *dFaceCell = NULL;
  int          *dFaceVtxIdx = NULL;
  PDM_g_num_t *dFaceVtx = NULL;
  double       *dVtxCoord = NULL;
  int          *dFaceGroupIdx = NULL;
  PDM_g_num_t *dFaceGroup = NULL;
  int           dFaceVtxL;
  int           dFaceGroupL;

  /*
   *  Create distributed cube
   */

  int          id;

  const double xmin = 0;
  const double ymin = 0;
  const double zmin = 0;

  const double xmax = xmin + length;
  const double ymax = ymin + length;
  const double zmax = zmin + length;

  if (myRank == 0) {
    printf("-- Build cube\n");
    fflush(stdout);
  }

  PDM_dcube_gen_init(&id,
                     PDM_MPI_COMM_WORLD,
                     nVtxSeg,
                     length,
                     xmin,
                     ymin,
                     zmin);

  PDM_dcube_gen_dim_get(id,
                      &nFaceGroup,
                      &dNCell,
                      &dNFace,
                      &dNVtx,
                      &dFaceVtxL,
                      &dFaceGroupL);

  PDM_dcube_gen_data_get(id,
                       &dFaceCell,
                       &dFaceVtxIdx,
                       &dFaceVtx,
                       &dVtxCoord,
                       &dFaceGroupIdx,
                       &dFaceGroup);
  int ppartId = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dCellPart = 0;

  int *dCellPart = (int *) malloc(dNCell*sizeof(int));

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int nPropertyCell = 0;
  int nPropertyFace = 0;

  if (myRank == 0) {
    printf("-- Part\n");
    fflush(stdout);
  }

  PDM_part_create(&ppartId,
                  PDM_MPI_COMM_WORLD,
                  method,
                  "PDM_PART_RENUM_CELL_NONE",
                  "PDM_PART_RENUM_FACE_NONE",
                  nPropertyCell,
                  renum_properties_cell,
                  nPropertyFace,
                  renum_properties_face,
                  nPart,
                  dNCell,
                  dNFace,
                  dNVtx,
                  nFaceGroup,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  have_dCellPart,
                  dCellPart,
                  dFaceCell,
                  dFaceVtxIdx,
                  dFaceVtx,
                  NULL,
                  dVtxCoord,
                  NULL,
                  dFaceGroupIdx,
                  dFaceGroup);

  free(dCellPart);



  /************************
   *
   * Point cloud definition
   *
   ************************/
  int _nPts_l;
  double *coords = NULL;
  _gen_cloud_random (nPts,
		     length,
		     numProcs,
		     myRank,
		     &coords,
		     &_nPts_l);

  int id_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD);

  double *char_length = malloc(sizeof(double) * _nPts_l);

  for (int i = 0; i < _nPts_l; i++) {
    char_length[i] = length * 1.e-6;
  }

  PDM_gnum_set_from_coords (id_gnum, 0, _nPts_l, coords, char_length);

  PDM_gnum_compute (id_gnum);

  PDM_g_num_t *gnum = PDM_gnum_get(id_gnum, 0);

  PDM_gnum_free (id_gnum, 1);



  
  /************************
   *
   * Mesh location struct initializaiton
   *
   ************************/

  int id_loc = PDM_mesh_location_create (PDM_MESH_NATURE_SURFACE_MESH,//???
					 1,//const int n_point_cloud,
					 PDM_MPI_COMM_WORLD);

  /* Set point cloud(s) */
  PDM_mesh_location_n_part_cloud_set (id_loc,
				      0,//i_point_cloud,
				      1);//n_part
  
  PDM_mesh_location_cloud_set (id_loc,
			       0,//i_point_cloud,
			       0,//i_part,
			       _nPts_l,
			       coords,
			       gnum);

  PDM_mesh_location_mesh_global_data_set (id_loc,
					  nPart);
     
  /* Set mesh */
  for (int ipart = 0; ipart < nPart; ipart++) {
    
    int nCell;
    int nFace;
    int nFacePartBound;
    int nVtx;
    int nProc;
    int nTPart;
    int sCellFace;
    int sFaceVtx;
    int sFaceGroup;
    int nEdgeGroup2;

    PDM_part_part_dim_get(ppartId,
                       ipart,
                       &nCell,
                       &nFace,
                       &nFacePartBound,
                       &nVtx,
                       &nProc,
                       &nTPart,
                       &sCellFace,
                       &sFaceVtx,
                       &sFaceGroup,
                       &nEdgeGroup2);

    int         *cellTag;
    int         *cellFaceIdx;
    int         *cellFace;
    PDM_g_num_t *cellLNToGN;
    int         *faceTag;
    int         *faceCell;
    int         *faceVtxIdx;
    int         *faceVtx;
    PDM_g_num_t *faceLNToGN;
    int         *facePartBoundProcIdx;
    int         *facePartBoundPartIdx;
    int         *facePartBound;
    int         *vtxTag;
    double      *vtx;
    PDM_g_num_t *vtxLNToGN;
    int         *faceGroupIdx;
    int         *faceGroup;
    PDM_g_num_t *faceGroupLNToGN;

    PDM_part_part_val_get (ppartId,
			   ipart,
			   &cellTag,
			   &cellFaceIdx,
			   &cellFace,
			   &cellLNToGN,
			   &faceTag,
			   &faceCell,
			   &faceVtxIdx,
			   &faceVtx,
			   &faceLNToGN,
			   &facePartBoundProcIdx,
			   &facePartBoundPartIdx,
			   &facePartBound,
			   &vtxTag,
			   &vtx,
			   &vtxLNToGN,
			   &faceGroupIdx,
			   &faceGroup,
			   &faceGroupLNToGN);

    PDM_mesh_location_part_set (id_loc,
				ipart,
				nCell,
				cellFaceIdx,
				cellFace,
				cellLNToGN,
				nFace,
				faceVtxIdx,
				faceVtx,
				faceLNToGN,
				nVtx,
				vtx,
				vtxLNToGN);
  }















  const double tolerance = 1.e-6;
  PDM_mesh_location_tolerance_set (id_loc,
				   tolerance);
  
  PDM_mesh_location_method_set (id_loc,
				PDM_MESH_LOCATION_OCTREE);


  

  PDM_mesh_location_compute (id_loc);

  PDM_g_num_t *location_elt_gnum = NULL;
  PDM_mesh_location_get (id_loc,
			 0,//i_point_cloud,
			 0,//i_part,
			 &location_elt_gnum);

#if 1
  /* Check results */
  if (myRank == 0) {
    printf("-- Check\n");
    fflush(stdout);
  }

  PDM_g_num_t nFaceSeg = nVtxSeg - 1;
  double cell_side = length / ((double) nFaceSeg);
  
  for (int ipt = 0; ipt < _nPts_l; ipt++) {
    double *pt_coord = coords + 3*ipt;

    int i = (int) floor (pt_coord[0] / cell_side);
    int j = (int) floor (pt_coord[1] / cell_side);
    int k = (int) floor (pt_coord[2] / cell_side);

    PDM_g_num_t box_gnum = 1 + i + nFaceSeg*(j + nFaceSeg*k);

    //printf("%d: (%ld) | (%ld)\n", ipt, location_elt_gnum[ipt], box_gnum);
    assert (location_elt_gnum[ipt] == box_gnum);
  }
#endif


  PDM_mesh_location_free (id_loc,
			  0);

  
  
  PDM_MPI_Finalize();

  if (myRank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }
   
  return 0;
}
