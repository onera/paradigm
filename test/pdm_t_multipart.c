#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_multipart.h"
#include "pdm_dcube_gen.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_dmesh.h"

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
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
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
           PDM_g_num_t  *nVtxSeg,
           double        *length,
           int           *n_part,
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
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

/**
 *
 * \Write data in a file in order to rebuild a CGNS
 *
 */
static void _dumpJsonData
(
 const char  *filename,
 const int    nVtx,
 const int    nFace,
 const int    nCell,
 const int    nFaceGroup,
 double      *VtxCoord,
 int         *FaceVtxIdx,
 PDM_g_num_t *FaceVtx,
 PDM_g_num_t *FaceCell,
 int         *cellFaceIdx,
 int         *cellFace,
 int         *FaceGroupIdx,
 PDM_g_num_t *FaceGroup
 )
{
  FILE * fp;
  fp = fopen(filename, "w");
  char endcomma = ' ';

  fprintf(fp, "{\n");
  fprintf(fp, "  \"nbVtx\" : %d,\n", nVtx);
  fprintf(fp, "  \"nbFace\" : %d,\n", nFace);
  fprintf(fp, "  \"nbCell\" : %d,\n", nCell);
  fprintf(fp, "  \"nbFaceGroup\" : %d,\n", nFaceGroup);

  fprintf(fp, "  \"VtxCoord\" :\n  [");
  for (int k=0; k<nVtx; k++)
  {
    fprintf(fp, "%c\n    [%12.5e,  %12.5e,  %12.5e]", endcomma, VtxCoord[3*k], VtxCoord[3*k+1], VtxCoord[3*k+2]);
    endcomma = ',';
  }
  fprintf(fp, "\n  ],\n");

  endcomma = ' ';
  fprintf(fp, "  \"NGonConnectivity\" :\n  [");
  for (int k=0; k<nFace; k++)
  {
    fprintf(fp, "%c\n    [", endcomma);
    for (int i=FaceVtxIdx[k]; i<FaceVtxIdx[k+1]-1; i++)
      fprintf(fp, "%d,  ", FaceVtx[i]);
    fprintf(fp, "%d]", FaceVtx[FaceVtxIdx[k+1]-1]);
    endcomma = ',';
  }
  fprintf(fp, "\n  ],\n");

  endcomma = ' ';
  fprintf(fp, "  \"NGonParentElement\" :\n  [");
  for (int k=0; k<nFace; k++)
  {
    fprintf(fp, "%c\n    [%d, %d]", endcomma, FaceCell[2*k], FaceCell[2*k+1]);
    endcomma = ',';
  }
  fprintf(fp, "\n  ],\n");

  if (cellFace != NULL)
  {
    endcomma = ' ';
    fprintf(fp, "  \"NFaceConnectivity\" :\n  [");
    for (int k=0; k<nCell; k++)
    {
      fprintf(fp, "%c\n    [", endcomma);
      for (int i=cellFaceIdx[k]; i<cellFaceIdx[k+1]-1; i++)
        fprintf(fp, "%d,  ", cellFace[i]);
      fprintf(fp, "%d]", cellFace[cellFaceIdx[k+1]-1]);
      endcomma = ',';
    }
    fprintf(fp, "\n  ],\n");
  }

  endcomma = ' ';
  fprintf(fp, "  \"FaceGroups\" :\n  [");
  for (int k=0; k<nFaceGroup; k++)
  {
    fprintf(fp, "%c\n    [", endcomma);
    for (int i=FaceGroupIdx[k]; i<FaceGroupIdx[k+1]-1; i++)
      fprintf(fp, "%d,  ", FaceGroup[i]);
    // Cas particulier ou le groupe n'a aucun élément
    if (FaceGroupIdx[k+1] > FaceGroupIdx[k])
      fprintf(fp, "%d]", FaceGroup[FaceGroupIdx[k+1]-1]);
    else
      fprintf(fp, "]");
    endcomma = ',';
  }
  fprintf(fp, "\n  ]\n");

  fprintf(fp, "}");
  fclose(fp);
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

  PDM_g_num_t        nVtxSeg = 10;
  double             length  = 1.;
  int                nPart   = 1;
  int                post    = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nVtxSeg,
             &length,
             &nPart,
             &post,
             (int *) &method);

  /*
   *  Init
   */

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

  int              cubeid;
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

  //Alloue les pointeurs de la structure interne dFaceCell, dFaceVtxIdx, dFaceVtx, dVtxCoord, dFaceGroupIdx, dFaceGroup
  PDM_dcube_gen_init(&cubeid,
                      comm,
                      nVtxSeg,
                      length,
            		      0.,
		                  0.,
		                  0.);

  PDM_dcube_gen_dim_get(cubeid,
                         &nFaceGroup,
                         &dNCell,
                         &dNFace,
                         &dNVtx,
                         &dFaceVtxL,
                         &dFaceGroupL);

  PDM_dcube_gen_data_get(cubeid,
                          &dFaceCell,
                          &dFaceVtxIdx,
                          &dFaceVtx,
                          &dVtxCoord,
                          &dFaceGroupIdx,
                          &dFaceGroup);

  int nblocks = 1;
  int mpartId;
  mpartId = PDM_multipart_create(nblocks, nPart, PDM_FALSE, method, comm);
  PDM_printf("From exe : created a multipart object, id is %i \n", mpartId);

  int *dmeshIds = (int *) malloc(nblocks * sizeof(int));
  for (int iblock = 0; iblock<nblocks; iblock++)
  {
    int dmeshId = 0;
    dmeshId = PDM_dmesh_create(dNCell, dNFace, dNVtx, nFaceGroup);
    PDM_dmesh_set(dmeshId, dVtxCoord, dFaceVtxIdx, dFaceVtx, dFaceCell, dFaceGroupIdx, dFaceGroup);
    PDM_multipart_register_block(mpartId, iblock, dmeshId);
    dmeshIds[iblock] = dmeshId;
  }


  char filename[25];
  snprintf(filename,sizeof(filename),"distributed_P%d.pdm",myRank);
  _dumpJsonData(filename, dNVtx, dNFace, dNCell, nFaceGroup,
            dVtxCoord, dFaceVtxIdx, dFaceVtx, dFaceCell, NULL, NULL, dFaceGroupIdx, dFaceGroup);

  PDM_multipart_run_ppart(mpartId);

  double  *elapsed = NULL;
  double  *cpu = NULL;
  double  *cpu_user = NULL;
  double  *cpu_sys = NULL;

  PDM_part_time_get(mpartId,
                 &elapsed,
                 &cpu,
                 &cpu_user,
                 &cpu_sys);

  if (myRank == 0)
  {
    PDM_printf("[%i]   - elapsed total                    : %12.5e\n", myRank, elapsed[0]);
    PDM_printf("[%i]   - elapsed building graph           : %12.5e\n", myRank, elapsed[1]);
    PDM_printf("[%i]   - elapsed splitting graph          : %12.5e\n", myRank, elapsed[2]);
    PDM_printf("[%i]   - elapsed building mesh partitions : %12.5e\n", myRank, elapsed[3]);

    PDM_printf("[%i]   - cpu total                        : %12.5e\n", myRank, cpu[0]);
    PDM_printf("[%i]   - cpu building graph               : %12.5e\n", myRank, cpu[1]);
    PDM_printf("[%i]   - cpu splitting graph              : %12.5e\n", myRank, cpu[2]);
    PDM_printf("[%i]   - cpu building mesh partitions     : %12.5e\n", myRank, cpu[3]);

    PDM_printf("[%i]   - cpu_user total                   : %12.5e\n", myRank, cpu_user[0]);
    PDM_printf("[%i]   - cpu_user building graph          : %12.5e\n", myRank, cpu_user[1]);
    PDM_printf("[%i]   - cpu_user splitting graph         : %12.5e\n", myRank, cpu_user[2]);
    PDM_printf("[%i]   - cpu_user building mesh partitions: %12.5e\n", myRank, cpu_user[3]);

    PDM_printf("[%i]   - cpu_sys total                    : %12.5e\n", myRank, cpu_sys[0]);
    PDM_printf("[%i]   - cpu_sys building graph           : %12.5e\n", myRank, cpu_sys[1]);
    PDM_printf("[%i]   - cpu_sys splitting graph          : %12.5e\n", myRank, cpu_sys[2]);
    PDM_printf("[%i]   - cpu_sys building mesh partitions : %12.5e\n", myRank, cpu_sys[3]);
  }

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
    int nFaceGroup2;

    PDM_multipart_part_dim_get(mpartId,
                       0,
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
                       &nFaceGroup2);

    int          *cellTag;
    int          *cellFaceIdx;
    int          *cellFace;
    PDM_g_num_t *cellLNToGN;
    int          *faceTag;
    int          *faceCell;
    int          *faceVtxIdx;
    int          *faceVtx;
    PDM_g_num_t *faceLNToGN;
    int          *facePartBoundProcIdx;
    int          *facePartBoundPartIdx;
    int          *facePartBound;
    int          *vtxTag;
    double       *vtx;
    PDM_g_num_t *vtxLNToGN;
    int          *faceGroupIdx;
    int          *faceGroup;
    PDM_g_num_t *faceGroupLNToGN;


    PDM_multipart_part_val_get(mpartId,
                       0,
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

    snprintf(filename,sizeof(filename),"partitionned_P%dN%d.pdm",myRank, ipart);
    _dumpJsonData(filename, nVtx, nFace, nCell, nFaceGroup2, vtx,
                  faceVtxIdx, faceVtx, faceCell, cellFaceIdx, cellFace,faceGroupIdx, faceGroup);

  }

  PDM_multipart_free(mpartId);

  for (int iblock = 0; iblock<nblocks; iblock++)
    PDM_dmesh_free(dmeshIds[iblock]);
  free(dmeshIds);

  //Desalloue les pointeurs de la structure interne dFaceCell, dFaceVtxIdx, dFaceVtx, dVtxCoord, dFaceGroupIdx, dFaceGroup
  PDM_dcube_gen_free(cubeid);

  PDM_MPI_Finalize();

  return 0;
}
