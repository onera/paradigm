#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include</stck/jcoulet/dev/dev-Test/json-c/install/include/json-c/json.h>

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
  #define MAXBUFLEN 10000000 // Taille max (en char --> wc file.dat) du ficher json
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
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scotch       Call PT-Scotch.\n\n"
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
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *nVtxSeg,
           double        *length,
           int           *n_part,
	         int           *method,
           char         **distri_dir)
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
    else if (strcmp(argv[i], "-distri_dir") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *distri_dir = (argv[i]);
      }
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
 const int    iblock,
 const int    nVtx,
 const int    nFace,
 const int    nCell,
 const int    nFaceGroup,
 const int    nTpart,
 double      *VtxCoord,
 int         *FaceVtxIdx,
 PDM_g_num_t *FaceVtx,
 PDM_g_num_t *FaceCell,
 int         *cellFaceIdx,
 int         *cellFace,
 int         *FaceGroupIdx,
 PDM_g_num_t *FaceGroup,
 int         *facePartBoundProcIdx,
 int         *facePartBoundPartIdx,
 int         *facePartBound

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
  fprintf(fp, "\n  ],\n");

  endcomma = ' ';
  fprintf(fp, "  \"PartBoundGroups\" :\n  {");
  for (int iMatch=0; iMatch < nTpart; iMatch++) //Local number of parts
  {
    if (facePartBoundPartIdx[iMatch+1] != facePartBoundPartIdx[iMatch])
    {
      int lastindex = -1;
      fprintf(fp, "%c\n    \"Match%i\" :\n    {", endcomma, iMatch);

      fprintf(fp, "\n      \"PointList\" : [");
      for (int k=facePartBoundPartIdx[iMatch]; k < facePartBoundPartIdx[iMatch+1]-1; k++)
        fprintf(fp, "%i, ", facePartBound[4*k] );
      lastindex = 4*(facePartBoundPartIdx[iMatch+1]-1);
      fprintf(fp, "%i]", facePartBound[lastindex]);

      fprintf(fp, ",\n      \"PointListDonor\" : [");
          for (int k=facePartBoundPartIdx[iMatch]; k < facePartBoundPartIdx[iMatch+1]-1; k++)
        fprintf(fp, "%i, ", facePartBound[4*k+3] );
      lastindex = 4*(facePartBoundPartIdx[iMatch+1]-1)+3;
      fprintf(fp, "%i]", facePartBound[lastindex]);

      char targetname[25];
      int procOpp = facePartBound[4*facePartBoundPartIdx[iMatch] + 1];
      int partOpp = facePartBound[4*facePartBoundPartIdx[iMatch] + 2];
      snprintf(targetname,sizeof(targetname),"P%dB%dN%d",procOpp, iblock, partOpp-1);
      fprintf(fp, ",\n      \"targetPart\" : \"%s\"", targetname);
      fprintf(fp, "\n    }");
    endcomma = ',';
    }
  }
  fprintf(fp, "\n  }\n");

  fprintf(fp, "}");
  fclose(fp);
}

static int _readJsonNumberOfBlocks
(
 const char *filename
)
{
  FILE *fp;
  char buffer[MAXBUFLEN];
  fp = fopen(filename, "r");
  fread(buffer, sizeof(char), MAXBUFLEN, fp);
  fclose(fp);

  struct json_object *parsed_json;
  struct json_object *intdata;
  int nbZones = -1;

  parsed_json = json_tokener_parse(buffer);
  json_object_object_get_ex(parsed_json, "nbZones", &intdata);
  nbZones = json_object_get_int(intdata);
  return nbZones;
}

static void _readJsonBlock
(
  const char  *filename,
  const int    blockId,
  int         *meshIds,
  int         *zoneIds
)
{
  PDM_printf("Parsing json data from %s\n", filename);
  FILE *fp;
  char buffer[MAXBUFLEN];
  fp = fopen(filename, "r");
  fread(buffer, sizeof(char), MAXBUFLEN, fp);
  fclose(fp);

  struct json_object *parsed_json;
  struct json_object *zonedata;
  struct json_object *intdata;
  struct json_object *array;
  struct json_object *arraydata;

  parsed_json = json_tokener_parse(buffer);
  // json_object_object_get_ex(parsed_json, "nbZones", &intdata);
  // PDM_printf("Nb of zones in json is %d\n\n", json_object_get_int(intdata));
  json_object_object_get_ex(parsed_json, "Zones", &zonedata);
  zonedata = json_object_array_get_idx(zonedata, blockId);

  int nbVtx, nbFace, nbCell, nbBound, nbJoin;
  int zoneGId = -1;
  json_object_object_get_ex(zonedata, "nbVtx", &intdata);
  nbVtx = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nbFace", &intdata);
  nbFace = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nbCell", &intdata);
  nbCell = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nbBoundary", &intdata);
  nbBound = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "nbInterface", &intdata);
  nbJoin = json_object_get_int(intdata);
  json_object_object_get_ex(zonedata, "ZoneGId", &intdata);
  zoneGId = json_object_get_int(intdata);

  // PDM_printf("Found in json : nbVtx = %d, nbFace = %d, nbCell = %d, nbBound = %d, \
  nbJoin = %d for zoneGId %d\n", nbVtx, nbFace, nbCell, nbBound, nbJoin, zoneGId);

  json_object_object_get_ex(zonedata, "VtxCoord", &array);
  assert(json_object_array_length(array) == nbVtx);
  double * vtxCoord = (double *) malloc(3*nbVtx * sizeof(double));
  for (int i = 0; i < nbVtx; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    vtxCoord[3*i] = json_object_get_double(json_object_array_get_idx(arraydata, 0));
    vtxCoord[3*i+1] = json_object_get_double(json_object_array_get_idx(arraydata, 1));
    vtxCoord[3*i+2] = json_object_get_double(json_object_array_get_idx(arraydata, 2));
  }

  json_object_object_get_ex(zonedata, "NGonConnectivity", &array);
  assert(json_object_array_length(array) == nbFace);
  int * faceVtxIdx = (int *) malloc((nbFace+1) * sizeof(int));
  faceVtxIdx[0] = 0;
  //First gets size of faceVtx and fills faceVtxIdx
  for (int i = 0; i < nbFace; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    faceVtxIdx[i+1] = json_object_array_length(arraydata) + faceVtxIdx[i];
  }
  PDM_g_num_t * faceVtx = (PDM_g_num_t *) malloc(faceVtxIdx[nbFace] * sizeof(PDM_g_num_t));
  //Now fill faceVtx
  int iface = 0;
  for (int i = 0; i < nbFace; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    for (int k = 0; k < json_object_array_length(arraydata); k++)
    {
      faceVtx[iface] = json_object_get_int(json_object_array_get_idx(arraydata, k));
      iface += 1;
    }
  }

  json_object_object_get_ex(zonedata, "NGonParentElement", &array);
  assert(json_object_array_length(array) == nbFace);
  PDM_g_num_t * faceCell = (PDM_g_num_t *) malloc(2*nbFace * sizeof(PDM_g_num_t));
  for (int i = 0; i < nbFace; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    faceCell[2*i] = json_object_get_int(json_object_array_get_idx(arraydata, 0));
    faceCell[2*i+1] = json_object_get_int(json_object_array_get_idx(arraydata, 1));
  }

  //First get size of faceBound and fills faceBoundIdx
  int * faceBoundIdx = (int *) malloc((nbBound + 1) * sizeof(int));
  faceBoundIdx[0] = 0;
  json_object_object_get_ex(zonedata, "BoundGroups", &array);
  assert(json_object_array_length(array) == nbBound);
  for (int i = 0; i < nbBound; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    faceBoundIdx[i+1] = json_object_array_length(arraydata) + faceBoundIdx[i];
  }
  //Now fill faceBound
  PDM_g_num_t * faceBound = (PDM_g_num_t *) malloc(faceBoundIdx[nbBound] * sizeof(PDM_g_num_t));
  iface = 0;
  for (int i = 0; i < nbBound; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    for (int k = 0; k < json_object_array_length(arraydata); k++)
    {
      faceBound[iface] = json_object_get_int(json_object_array_get_idx(arraydata, k));
      iface += 1;
    }
  }

  //First get size of faceJoin and fills faceJoinIdx
  int * faceJoinIdx = (int *) malloc((nbJoin + 1) * sizeof(int));
  faceJoinIdx[0] = 0;
  json_object_object_get_ex(zonedata, "InterfaceGroups", &array);
  assert(json_object_array_length(array) == nbJoin);
  for (int i = 0; i < nbJoin; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    json_object_object_get_ex(arraydata, "PointList", &arraydata);
    faceJoinIdx[i+1] = json_object_array_length(arraydata) + faceJoinIdx[i];
  }
  //Now fill faceJoin
  PDM_g_num_t * faceJoin = (PDM_g_num_t *) malloc(faceJoinIdx[nbJoin] * sizeof(PDM_g_num_t));
  int * joinZoneOpp = (int *) malloc(nbJoin * sizeof(int));
  iface = 0;
  for (int i = 0; i < nbJoin; i++)
  {
    arraydata = json_object_array_get_idx(array, i);
    json_object_object_get_ex(arraydata, "targetGId", &intdata);
    joinZoneOpp[i] = json_object_get_int(intdata);
    json_object_object_get_ex(arraydata, "PointList", &arraydata);
    for (int k = 0; k < json_object_array_length(arraydata); k++)
    {
      faceJoin[iface] = json_object_get_int(json_object_array_get_idx(arraydata, k));
      iface += 1;
    }
  }

  int dmeshId = PDM_dmesh_create(nbCell, nbFace, nbVtx, nbBound, nbJoin);
  PDM_dmesh_set(dmeshId, vtxCoord, faceVtxIdx, faceVtx, faceCell,
                faceBoundIdx, faceBound, joinZoneOpp, faceJoinIdx, faceJoin, NULL);
  meshIds[blockId] = dmeshId;
  zoneIds[blockId] = zoneGId;
  //Fuite mémoire -> les données allouées ici sont perdues
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
  char               *distri_dir = strdup("\0");

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
             (int *) &method,
             &distri_dir);

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

  int nbzone = 1;
  int mpartId = -1;
  int multipartSize = 1;
  int *dmeshIds;             // For a given block, store the identifier of the corresponding blockdata (size = nzone)
  int *dblockIds;            // For a given block, store the globalId of the parent zone

  if (strcmp(distri_dir, "\0") == 0) // No distributed data provided -> generate cube (only 1 zone)
  {
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

    dmeshIds  = (int *) malloc(nbzone * sizeof(int));
    dblockIds = (int *) malloc(nbzone * sizeof(int));
    int dmeshId = -1;
    dmeshId = PDM_dmesh_create(dNCell, dNFace, dNVtx, nFaceGroup, 0);
    PDM_dmesh_set(dmeshId, dVtxCoord, dFaceVtxIdx, dFaceVtx, dFaceCell, dFaceGroupIdx, dFaceGroup, NULL, NULL, NULL, NULL);
    dmeshIds[0] = dmeshId;
    dblockIds[0] = 1;

    // Dump arrays
    char filename[25];
    snprintf(filename,sizeof(filename),"distributed_P%d.json",myRank);
    _dumpJsonData(filename, 0, dNVtx, dNFace, dNCell, nFaceGroup, 0,
                  dVtxCoord, dFaceVtxIdx, dFaceVtx, dFaceCell, NULL, NULL, dFaceGroupIdx, dFaceGroup,
                  NULL, NULL, NULL);
  }
  else  // -> parse Json data
  {
    char filename[25];
    snprintf(filename,sizeof(filename),"distributed_P%d.json",myRank);
    char path[50];
    strcpy(path, distri_dir);
    strcat(path, filename);
    // 1. Get number of zones known by each zone
    nbzone = _readJsonNumberOfBlocks(path);
    PDM_printf("[%i] Found %d blocks in distributed data\n", myRank, nbzone);
    dmeshIds  = (int *) malloc(nbzone * sizeof(int));
    dblockIds = (int *) malloc(nbzone * sizeof(int));

    // 2. Read the data for each zone and store it in arrays
    for (int iblock=0; iblock < nbzone; iblock++)
    {
      _readJsonBlock(path, iblock, dmeshIds, dblockIds);
      // PDM_dmesh_dims_get(dmeshId, &dNCell, &dNFace, &dNVtx, &nFaceGroup);
      // PDM_dmesh_data_get(dmeshId, &dVtxCoord, &dFaceVtxIdx, &dFaceVtx, &dFaceCell, &dFaceGroupIdx, &dFaceGroup);
    }
    // 3. Compute the multizone size (total number of different zones)
    int maxZoneGId = -1;
    for (int k = 0; k < nbzone; k++)
      maxZoneGId = (dblockIds[k] > maxZoneGId) ? dblockIds[k] : maxZoneGId;
    PDM_MPI_Allreduce(&maxZoneGId, &multipartSize, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  }

  int * nPartArray = (int *) malloc((multipartSize) * sizeof(int));
  for (int k=0; k<multipartSize; k++)
    nPartArray[k] = nPart;
   mpartId = PDM_multipart_create(multipartSize, nPartArray, PDM_FALSE, method, comm);
   PDM_printf("From exe : created a multipart object, id is %i \n", mpartId);
   for (int iblock=0; iblock < nbzone; iblock++)
   {
    // PDM_printf("[%i] check -- iblock %i : gid %i meshid %i\n", myRank, iblock, dblockIds[iblock], dmeshIds[iblock]);
    PDM_multipart_register_block(mpartId, dblockIds[iblock]-1, dmeshIds[iblock]);
   }
   PDM_MPI_Barrier(comm);



  PDM_multipart_run_ppart(mpartId);

  double  *elapsed = NULL;
  double  *cpu = NULL;
  double  *cpu_user = NULL;
  double  *cpu_sys = NULL;

  PDM_multipart_time_get(mpartId,
                         multipartSize-1,
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

  for (int iblock = 0; iblock < multipartSize; iblock++) {
    for (int ipart = 0; ipart < nPartArray[iblock]; ipart++) {

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
                                 iblock,
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
                                 iblock,
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

      char filename[25];
      snprintf(filename,sizeof(filename),"partitionned_P%dB%dN%d.json",myRank, iblock, ipart);
      _dumpJsonData(filename, iblock, nVtx, nFace, nCell, nFaceGroup2, nTPart, vtx,
                    faceVtxIdx, faceVtx, faceCell, cellFaceIdx, cellFace,faceGroupIdx, faceGroup,
                    facePartBoundProcIdx, facePartBoundPartIdx, facePartBound);

    }
  }

  PDM_multipart_free(mpartId);

  for (int iblock = 0; iblock<nbzone; iblock++)
    PDM_dmesh_free(dmeshIds[iblock]);
  free(dmeshIds);
  free(dblockIds);
  free(nPartArray);

  if (strcmp(distri_dir, "\0") == 0) // No distributed data provided
  {
  //Desalloue les pointeurs de la structure interne dFaceCell, dFaceVtxIdx, dFaceVtx, dVtxCoord, dFaceGroupIdx, dFaceGroup
    PDM_dcube_gen_free(cubeid);
  }

  PDM_MPI_Finalize();

  return 0;
}
